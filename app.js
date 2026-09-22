'use strict';
// ABIF v1 parser: no network, no external dependencies. Untrusted inputs are bounds-checked.
const $=id=>document.getElementById(id);
const state={ref:null,reads:[],rows:[],view:null,drag:null};
const rc=s=>s.split('').reverse().map(x=>({A:'T',T:'A',C:'G',G:'C',R:'Y',Y:'R',M:'K',K:'M',S:'S',W:'W',N:'N'}[x]||'N')).join('');
const valid=s=>/^[ACGTNRYKMSWBDHV]+$/.test(s);
function parseABI(buffer){
 const d=new DataView(buffer),u8=new Uint8Array(buffer),dec=new TextDecoder('ascii');
 function check(p,n){if(!Number.isSafeInteger(p)||!Number.isSafeInteger(n)||p<0||n<0||p+n>buffer.byteLength)throw Error('Fichier ABI tronqué ou invalide');}
 function u16(p){check(p,2);return d.getUint16(p,false)} function u32(p){check(p,4);return d.getUint32(p,false)}
 check(0,34);if(dec.decode(u8.slice(0,4))!=='ABIF')throw Error('Signature ABIF absente');
 const count=u32(18),offset=u32(26);if(count>10000)throw Error('Répertoire ABIF trop grand');check(offset,count*28);
 const entries=new Map();for(let i=0;i<count;i++){let p=offset+i*28;let name=dec.decode(u8.slice(p,p+4)),num=u32(p+4),type=u16(p+8),size=u16(p+10),n=u32(p+12),bytes=u32(p+16),at=bytes<=4?p+20:u32(p+20);check(at,bytes);entries.set(name+num,{type,size,n,bytes,at})}
 function tag(name,num){const e=entries.get(name+num);if(!e)throw Error('Tag ABI absent : '+name+num);return e}
 function chars(name,num){let e=tag(name,num);return dec.decode(u8.slice(e.at,e.at+e.bytes)).replace(/\0/g,'')}
 function numbers(name,num){let e=tag(name,num);if(e.n>2e7||e.size<1||e.n*e.size>e.bytes)throw Error('Tag numérique ABI invalide');let out=new Array(e.n);for(let i=0;i<e.n;i++){let p=e.at+i*e.size;out[i]=e.size===2?u16(p):e.size===4?u32(p):u8[p]}return out}
 let bases=chars('PBAS',2);if(!bases||!valid(bases))bases=chars('PBAS',1);let positions=numbers('PLOC',2);if(positions.length!==bases.length)positions=numbers('PLOC',1);
 let order=chars('FWO_',1).slice(0,4);if(!/^[ACGT]{4}$/.test(order)||new Set(order).size!==4)throw Error('Ordre des canaux FWO_ invalide');
 let traces={};for(let i=0;i<4;i++)traces[order[i]]=numbers('DATA',9+i);
 let n=traces.A.length;if(!n||Object.values(traces).some(t=>t.length!==n))throw Error('Longueurs des canaux incohérentes');
 let quality=[];try{quality=numbers('PCON',2)}catch(e){};
 return {bases,positions,traces,quality,order};
}
function fasta(text){let lines=text.trim().split(/\r?\n/);let header=lines.find(x=>x.startsWith('>'))||'';let seq=lines.filter(x=>!x.startsWith('>')).join('').replace(/\s/g,'').toUpperCase().replace(/U/g,'T');if(!seq||!/^[ACGTN]+$/.test(seq))throw Error('FASTA invalide (séquence ADN attendue)');if(!header.includes('NM_000059.4'))throw Error('En-tête FASTA : NM_000059.4 attendu. Vérifiez la version de référence.');return {seq,header}}
$('loadref').onclick=async()=>{try{let file=$('fasta').files[0];let text=file?await file.text():$('refpaste').value;state.ref=fasta(text);$('refstatus').textContent=`Référence chargée : ${state.ref.header}\n${state.ref.seq.length} nucléotides.`}catch(e){state.ref=null;$('refstatus').textContent='Erreur : '+e.message}};
// Référence statique servie par GitHub Pages : seuls les octets publics du FASTA sont téléchargés.
// Aucun chromatogramme n'est transmis au serveur.
(async function autoReference(){try{
 const response=await fetch('./BRCA2_NM_000059.4.fasta',{cache:'no-store'});
 if(!response.ok)throw Error('Fichier de référence absent du dépôt (HTTP '+response.status+')');
 const loaded=fasta(await response.text());
 if(loaded.seq.length!==11954)throw Error('Longueur inattendue : '+loaded.seq.length+' (11954 attendus)');
 state.ref=loaded; $('refstatus').textContent='Référence BRCA2 NM_000059.4 intégrée : '+loaded.seq.length+' nucléotides. Aucune donnée patient envoyée.';
 }catch(e){$('refstatus').textContent='Référence intégrée indisponible : '+e.message+' — chargez votre FASTA manuellement.'}
})();
// Banded-free Smith-Waterman local alignment with traceback. Limit read/reference size to bound runtime.
function align(read,ref){const n=read.length,m=ref.length;if(n>3000||m>25000)throw Error('Référence/lecture trop longue pour cet alignement de démonstration');
 let w=m+1,H=new Int32Array((n+1)*w),T=new Uint8Array(H.length),best=0,bi=0,bj=0;
 for(let i=1;i<=n;i++){let row=i*w,prev=(i-1)*w;for(let j=1;j<=m;j++){let match=(read[i-1]===ref[j-1]&&'ACGT'.includes(read[i-1]))?3:-3;let a=H[prev+j-1]+match,b=H[prev+j]-4,c=H[row+j-1]-4;let v=Math.max(0,a,b,c);H[row+j]=v;T[row+j]=v===0?0:v===a?1:v===b?2:3;if(v>best){best=v;bi=i;bj=j}}}
 let i=bi,j=bj,cols=[];while(i>0&&j>0&&H[i*w+j]>0){let t=T[i*w+j];if(t===1){cols.push({q:i-1,r:j-1,base:read[i-1],ref:ref[j-1]});i--;j--}else if(t===2){cols.push({q:i-1,r:null,base:read[i-1],ref:'-'});i--}else if(t===3){cols.push({q:null,r:j-1,base:'-',ref:ref[j-1]});j--}else break}cols.reverse();return {score:best,cols,start:i,end:bi,refStart:j,refEnd:bj}}
function chooseAlignment(read,ref){let a=align(read,ref),b=align(rc(read),ref);return b.score>a.score?{...b,strand:'-',oriented:rc(read),length:read.length}:{...a,strand:'+',oriented:read,length:read.length}}
function cpos(refZero,cds){let x=refZero+1;if(!cds||cds<1)return 'non défini';let d=x-cds+1;return d>0?'c.'+d:'c.'+d}
function callRows(read){let a=read.alignment,cds=Number($('cdsstart').value),rows=[];for(let k=0;k<a.cols.length;k++){let col=a.cols[k];if(col.q===null||col.r===null||col.base===col.ref||!/[ACGT]/.test(col.base)||!/[ACGT]/.test(col.ref))continue;
 let qOriginal=a.strand==='+'?col.q:a.length-1-col.q;let label=Number.isInteger(cds)&&cds>0?`NM_000059.4:${cpos(col.r,cds)}${col.ref}>${col.base}`:'CDS non renseignée';
 rows.push({name:read.name,q:qOriginal+1,r:col.r+1,ref:col.ref,alt:col.base,hgvs:label,read})}return rows}
$('analyze').onclick=async()=>{state.reads=[];state.rows=[];$('variants').replaceChildren();$('alignments').replaceChildren();$('readselect').replaceChildren();try{
 for(let [id,name] of [['forward','Forward'],['reverse','Reverse']]){let f=$(id).files[0];if(!f)continue;let ab=parseABI(await f.arrayBuffer());let read={...ab,name:name+' · '+f.name};if(state.ref){read.alignment=chooseAlignment(ab.bases,state.ref.seq);state.rows.push(...callRows(read))}state.reads.push(read)}if(!state.reads.length)throw Error('Importez au moins un fichier ABI');
 for(let [i,r] of state.reads.entries()){let opt=document.createElement('option');opt.value=i;opt.textContent=r.name;$('readselect').append(opt);let p=document.createElement('p');p.textContent=`${r.name} — ${r.bases.length} bases${r.alignment?`, alignement ${r.alignment.strand} ; référence ${r.alignment.refStart+1}–${r.alignment.refEnd} ; score ${r.alignment.score}`:''}`;$('alignments').append(p)}
 for(let row of state.rows.slice(0,1000)){let tr=document.createElement('tr');for(let value of [row.name,row.q,row.r,row.ref,row.alt,row.hgvs]){let td=document.createElement('td');td.textContent=value;tr.append(td)}tr.onclick=()=>{let idx=state.reads.indexOf(row.read);$('readselect').value=idx;$('start').value=Math.max(1,row.q-20);draw()};$('variants').append(tr)}
 $('status').textContent=`${state.reads.length} chromatogramme(s) lu(s) ; ${state.rows.length} discordance(s) simples indicatives. ${state.ref?'Référence présente.':'Aucune référence : visualisation uniquement.'}`;draw();
 }catch(e){$('status').textContent='Erreur : '+e.message}}
function draw(){let r=state.reads[Number($('readselect').value)||0];if(!r)return;let start=Math.max(0,Math.min(r.bases.length-5,(Number($('start').value)||1)-1)),count=Math.min(r.bases.length,Math.max(8,Number($('count').value)||45));let end=Math.min(r.bases.length,start+count);if(end<=start)return;
 let canvas=$('trace'),ctx=canvas.getContext('2d'),W=canvas.width,H=canvas.height;ctx.clearRect(0,0,W,H);let x0=Math.max(0,r.positions[start]-12),x1=Math.min(r.traces.A.length-1,r.positions[end-1]+12),span=Math.max(1,x1-x0);let ymax=1;for(let b of 'ACGT'){let t=r.traces[b];for(let i=x0;i<=x1;i++)if(t[i]>ymax)ymax=t[i]};let color={A:'#16854e',C:'#2166cc',G:'#252b32',T:'#d33a37'};ctx.lineWidth=1.5;
 for(let b of 'ACGT'){ctx.beginPath();ctx.strokeStyle=color[b];let t=r.traces[b];for(let i=x0;i<=x1;i++){let x=(i-x0)/span*W,y=H-45-t[i]/ymax*(H-80);if(i===x0)ctx.moveTo(x,y);else ctx.lineTo(x,y)}ctx.stroke()}
 ctx.font='13px Segoe UI';ctx.textAlign='center';for(let k=start;k<end;k++){let x=(r.positions[k]-x0)/span*W;ctx.fillStyle=color[r.bases[k]]||'#8c60a7';ctx.fillText(r.bases[k],x,H-22);if(k%5===0){ctx.fillStyle='#647587';ctx.font='10px Segoe UI';ctx.fillText(String(k+1),x,H-6);ctx.font='13px Segoe UI'}}state.view={r,x0,x1,start,end};$('baseinfo').textContent=`${r.name} · bases ${start+1}–${end} · cliquer sur un pic pour afficher sa position.`}
function viewTo(start,count){let r=state.reads[Number($('readselect').value)||0];if(!r)return;count=Math.max(8,Math.min(r.bases.length,Math.round(count)));start=Math.max(1,Math.min(r.bases.length-count+1,Math.round(start)));$('start').value=start;$('count').value=count;draw()}
$('zoomout').onclick=()=>viewTo(Number($('start').value)-Number($('count').value)*.25,Number($('count').value)*1.5);
$('resetview').onclick=()=>{let r=state.reads[Number($('readselect').value)||0];if(r)viewTo(1,r.bases.length)};
const trace=$('trace');
trace.addEventListener('pointerdown',e=>{if(!state.view)return;state.drag={x:e.clientX,start:Number($('start').value),moved:false};trace.setPointerCapture(e.pointerId);trace.style.cursor='grabbing'});
trace.addEventListener('pointermove',e=>{if(!state.drag)return;let dx=e.clientX-state.drag.x;if(Math.abs(dx)>3)state.drag.moved=true;if(!state.drag.moved)return;let count=Number($('count').value),rect=trace.getBoundingClientRect();viewTo(state.drag.start-Math.round(dx/rect.width*count),count)});
trace.addEventListener('pointerup',()=>{trace.style.cursor='grab';if(state.drag){state.suppressClick=state.drag.moved;state.drag=null;setTimeout(()=>state.suppressClick=false,80)}});
trace.addEventListener('dblclick',e=>{e.preventDefault();let v=state.view;if(!v)return;let rect=trace.getBoundingClientRect(),fraction=(e.clientX-rect.left)/rect.width,count=Number($('count').value),next=Math.max(8,Math.round(count/2));viewTo(Number($('start').value)+Math.round(fraction*(count-next)),next)});
trace.addEventListener('wheel',e=>{e.preventDefault();let v=state.view;if(!v)return;let rect=trace.getBoundingClientRect(),fraction=(e.clientX-rect.left)/rect.width,count=Number($('count').value),next=Math.round(count*(e.deltaY<0?.8:1.25));viewTo(Number($('start').value)+Math.round(fraction*(count-next)),next)},{passive:false});
$('draw').onclick=draw;$('readselect').onchange=draw;$('trace').onclick=e=>{if(state.suppressClick)return;let v=state.view;if(!v)return;let rect=$('trace').getBoundingClientRect(),sample=v.x0+(e.clientX-rect.left)/rect.width*(v.x1-v.x0),k=v.start;for(let i=v.start;i<v.end;i++)if(Math.abs(v.r.positions[i]-sample)<Math.abs(v.r.positions[k]-sample))k=i;let a=v.r.alignment,ref='non alignée';if(a){let q=a.strand==='+'?k:a.length-1-k;let c=a.cols.find(c=>c.q===q&&c.r!==null);if(c)ref=`référence ${c.r+1} ; ${cpos(c.r,Number($('cdsstart').value))}`} $('baseinfo').textContent=`${v.r.name} · base ABI ${k+1} : ${v.r.bases[k]} · ${ref}`};
$('export').onclick=()=>{let header=['lecture','position_ABI','position_reference_1based','reference','observe','HGVS_indicatif'];let esc=v=>'"'+String(v).replace(/"/g,'""')+'"';let data=[header,...state.rows.map(r=>[r.name,r.q,r.r,r.ref,r.alt,r.hgvs])].map(row=>row.map(esc).join(',')).join('\r\n');let blob=new Blob(['\ufeff'+data],{type:'text/csv;charset=utf-8'}),url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;a.download='sangervariant_discordances_indicatives.csv';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000)};

// Score heuristique : recherche d'une augmentation des signaux secondaires compatibles
// avec une superposition de la référence WT et d'une référence décalée de 4 bases.
// Le score ne constitue pas une probabilité, un génotype ou une validation d'indel.
function candidateScore(read,breakQ,shift){
 const a=read.alignment,ref=state.ref.seq,points=[];
 const qmap=new Map(a.cols.filter(c=>c.q!==null&&c.r!==null).map(c=>[c.q,c.r]));
 for(let q=breakQ+2;q<Math.min(breakQ+28,a.oriented.length-2);q++){
   const r=qmap.get(q);if(r===undefined||r+shift<0||r+shift>=ref.length)continue;
   const original=a.strand==='+'?q:a.length-1-q;
   const pos=read.positions[original];if(!Number.isInteger(pos)||pos<0)continue;
   const observed={};let total=0;
   for(let base of 'ACGT'){
     const v=read.traces[base][pos]||0;observed[base]=v;total+=v;
   }
   if(total<=0)continue;
   const wt=ref[r],del=ref[r+shift];if(!'ACGT'.includes(wt)||!'ACGT'.includes(del)||wt===del)continue;
   const expected=a.strand==='+'?del:rc(del);
   const normal=a.strand==='+'?wt:rc(wt);
   const mixed=observed[expected]/total;
   const primary=observed[normal]/total;
   points.push({mixed,primary});
 }
 if(points.length<8)return null;
 return {score:points.reduce((s,x)=>s+x.mixed,0)/points.length,
   primary:points.reduce((s,x)=>s+x.primary,0)/points.length,n:points.length};
}
$('scanindel').onclick=()=>{
 const out=$('indelrows');out.replaceChildren();
 if(!state.ref){$('indelstatus').textContent='Chargez la référence BRCA2 avant la recherche.';return}
 let all=[];
 for(const read of state.reads){
  const a=read.alignment;if(!a)continue;
  const candidates=[];
  for(let q=Math.max(a.start+12,15);q<Math.min(a.end-32,a.length-32);q+=2){
   const r=a.cols.find(c=>c.q===q&&c.r!==null)?.r;
   if(r===undefined)continue;
   // Les deux sens sont évalués : un indel hétérozygote peut être ambigu
   // en présence de pics mixtes et d'un alignement déjà dégradé.
   for(const shift of [-4,4]){
    const result=candidateScore(read,q,shift);
    if(result)candidates.push({read,q,r,shift,...result});
   }
  }
  candidates.sort((x,y)=>y.score-x.score);
  // Regrouper les maxima proches : ne pas afficher dix positions contiguës.
  let chosen=[];
  for(const item of candidates){if(chosen.every(x=>Math.abs(x.q-item.q)>14)){chosen.push(item);if(chosen.length>=5)break}}
  all.push(...chosen);
 }
 all.sort((x,y)=>y.score-x.score);
 for(const item of all){
  const tr=document.createElement('tr');
  const original=item.read.alignment.strand==='+'?item.q+1:item.read.bases.length-item.q;
  for(const val of [item.read.name,original,item.r+1,'Décalage '+(item.shift>0?'+':'')+item.shift+' bases (hypothèse)',item.score.toFixed(3)+' (n='+item.n+')']){
   const td=document.createElement('td');td.textContent=val;tr.append(td)
  }
  tr.onclick=()=>{ $('readselect').value=state.reads.indexOf(item.read);viewTo(original-25,55)};
  out.append(tr)
 }
 $('indelstatus').textContent=all.length?`${all.length} régions exploratoires classées par score de signal secondaire. Les scores ne sont pas calibrés et ne permettent PAS d'identifier ni de confirmer une délétion de 4 bases.`:'Aucune région évaluable : couverture ou alignement insuffisant.';
};
