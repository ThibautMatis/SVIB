'use strict';
/* SVIB V7: independent, reference-guided two-component trace fitting.
 * Exploratory research implementation, NOT Tracy, NOT validated variant calling.
 * All AB1 processing is local. No known variant or expected ABI position is embedded.
 */
(()=>{
const $=id=>document.getElementById(id), DNA='ACGT', COMP={A:'T',C:'G',G:'C',T:'A'};
const results=[];
const median=a=>{if(!a.length)return NaN;const b=[...a].sort((x,y)=>x-y);return b[Math.floor(b.length/2)]};
function profile(read,k){if(k<0||k>=read.positions.length)return null;const p=read.positions[k],v=Array.from(DNA,b=>Math.max(0,read.traces[b][p]||0)),s=v.reduce((a,b)=>a+b,0);return s<50?null:v.map(x=>x/s)}
function refAt(pos,dir){const b=state.ref.seq[pos];return b?(dir===1?b:COMP[b]||'N'):'N'}
function err(v,a,b,f){if(!v||!DNA.includes(a)||!DNA.includes(b))return NaN;let e=0;for(let i=0;i<4;i++){const x=.015+.94*((a===DNA[i]?1-f:0)+(b===DNA[i]?f:0));e+=(v[i]-x)**2}return e}
// Anchor using ONLY high-quality, unambiguous matches upstream of a proposed transition.
function anchor(read,k){const al=read.alignment;if(!al)return null;const dir=al.strand==='+'?1:-1, offsets=[];
 for(const c of al.cols){if(c.q===null||c.r===null||c.base!==c.ref)continue;const q=dir===1?c.q:read.bases.length-1-c.q;
 if(q>=k-62&&q<=k-10&&DNA.includes(read.bases[q]))offsets.push(c.r-dir*q)}
 if(offsets.length<12)return null;const off=median(offsets);let matches=0,compared=0;
 for(let q=Math.max(0,k-38);q<k-7;q++){const b=refAt(off+dir*q,dir);if(DNA.includes(b)){compared++;if(read.bases[q]===b)matches++}}
 return compared>=22&&matches/compared>=.68?{dir,off,matches,compared}:null;
}
function fit(read,k,a){const {dir,off}=a;let pre=0,np=0;for(let q=k-23;q<k-5;q++){const b=refAt(off+dir*q,dir),e=err(profile(read,q),b,b,0);if(Number.isFinite(e)){pre+=e;np++}}
 if(np<13||pre/np>.25)return null;
 const post=[];for(let q=k+4;q<Math.min(k+44,read.bases.length-3);q++){const b=refAt(off+dir*q,dir),v=profile(read,q);if(v&&DNA.includes(b))post.push({q,b,v})}
 if(post.length<23)return null;
 const nullErr=post.reduce((s,x)=>s+err(x.v,x.b,x.b,0),0)/post.length;
 let best=null;
 // Shift is in the direction of ABI acquisition; sign alone does NOT determine insertion/deletion.
 for(let shift=-20;shift<=20;shift++){if(!shift)continue;
  for(let frac=.25;frac<=.751;frac+=.05){let total=0,n=0;
   for(const x of post){if(x.q-k<=Math.abs(shift)+2)continue;const alt=refAt(off+dir*(x.q+shift),dir),e=err(x.v,x.b,alt,frac);if(Number.isFinite(e)){total+=e;n++}}
   if(n<19)continue;const loss=total/n;
   // Penalize shifts whose compared positions contain little informative sequence.
   const gain=nullErr-loss;
   if(!best||gain>best.gain)best={shift,frac,loss,gain,n};
  }
 }
 if(!best||best.gain<.045||best.loss>.24||best.gain/nullErr<.18)return null;
 return {...best,read,k,dir,off,refPos:off+dir*k+1,preLoss:pre/np,nullErr,anchorMatches:a.matches};
}
function scan(read){const all=[];for(let k=45;k<read.bases.length-51;k+=1){const a=anchor(read,k);if(!a)continue;const c=fit(read,k,a);if(c)all.push(c)}
 // Collapse neighboring windows, preserving the best fitting transition per region.
 all.sort((a,b)=>b.gain-a.gain);const kept=[];
 for(const c of all){if(kept.every(x=>Math.abs(x.k-c.k)>23))kept.push(c);if(kept.length>=8)break}
 return kept;
}
function seq(read,c){const {k,off,dir,shift,frac}=c,one=[],two=[],bases=[];
 for(let q=Math.max(0,k-16);q<Math.min(read.bases.length,k+52);q++){
  const b=refAt(off+dir*q,dir),v=profile(read,q),altRef=refAt(off+dir*(q+shift),dir);
  one.push(b);
  if(q<k){two.push(b);continue}
  if(!v||!DNA.includes(b)){two.push('N');continue}
  // Reference-guided residual estimate: infer secondary base from measured signal.
  const residual=DNA.map((x,i)=>({base:x,score:v[i]-(x===b?.94*(1-frac):0)})).sort((a,b)=>b.score-a.score);
  const picked=residual[0];two.push(picked.score>.09?picked.base:'N');
  if(q>=k+Math.abs(shift)+3&&DNA.includes(altRef))bases.push(picked.base===altRef?1:0);
 }
 return {wt:one.join(''),alt:two.join(''),support:bases.length?bases.reduce((a,b)=>a+b,0)/bases.length:NaN};
}
function show(c){$('readselect').value=state.reads.indexOf(c.read);viewTo(Math.max(1,c.k-23),66);const s=seq(c.read,c);
 $('v7detail').textContent=`Lecture : ${c.read.name}\nTransition estimée : base ABI ${c.k+1} ; référence ≈ ${c.refPos} (non normalisée)\nOrientation : ${c.dir===1?'+':'−'} ; décalage dans le sens d’acquisition : ${c.shift>0?'+':''}${c.shift} nt\nFraction de mélange modélisée : ${c.frac.toFixed(2)} ; gain non calibré : ${c.gain.toFixed(3)}\nConcordance des bases secondaires avec le décalage : ${Number.isFinite(s.support)?(100*s.support).toFixed(0)+' %':'indéterminée'}\n\nAllèle de référence (fenêtre) :\n${s.wt}\nAllèle secondaire estimé (fenêtre, N = incertain) :\n${s.alt}\n\nATTENTION : estimation locale guidée par la référence, pas un alignement complet des deux allèles. La taille, le type d’indel et la nomenclature HGVS ne sont PAS déterminés. Les résultats ne sont pas validés pour le diagnostic.`;
}
function run(){results.length=0;$('v7rows').replaceChildren();$('v7detail').textContent='';const status=$('v7status');if(!state.ref||!state.reads.length||state.reads.some(r=>!r.alignment)){status.textContent='Chargez un transcrit, importez les ABI puis cliquez sur « Lire et analyser localement ».';return}
 status.textContent='Analyse indépendante en cours…';setTimeout(()=>{try{for(const read of state.reads)results.push(...scan(read));results.sort((a,b)=>b.gain-a.gain);for(const c of results){const tr=document.createElement('tr');for(const v of [c.read.name,c.k+1,c.refPos,c.shift,c.gain.toFixed(3),c.frac.toFixed(2)]){const td=document.createElement('td');td.textContent=v;tr.append(td)}tr.onclick=()=>show(c);$('v7rows').append(tr)}status.textContent=results.length?`${results.length} régions exploratoires détectées. Chaque orientation a été analysée indépendamment. Aucun variant HGVS n’est confirmé.`:'Aucune région ne passe les critères exploratoires. Un résultat vide n’exclut pas un variant.'}catch(e){status.textContent='Erreur V7 : '+e.message;console.error(e)}},30);
}
$('v7scan').addEventListener('click',run);
$('v7export').addEventListener('click',()=>{const rows=[['lecture','base_ABI','position_NM_approximative','decalage_acquisition','gain_non_calibre','fraction_modele'],...results.map(c=>[c.read.name,c.k+1,c.refPos,c.shift,c.gain,c.frac])];const csv=rows.map(row=>row.map(v=>'"'+String(v).replace(/"/g,'""')+'"').join(',')).join('\r\n');const url=URL.createObjectURL(new Blob(['\ufeff'+csv],{type:'text/csv;charset=utf-8'}));const a=document.createElement('a');a.href=url;a.download='SVIB_V7_exploration.csv';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000)});
})();
