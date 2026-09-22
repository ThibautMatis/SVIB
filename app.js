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
// Recherche publique NCBI : aucune séquence de patient, chromatogramme ou variant n'entre dans ces URL.
const API='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
let refRequest=0;
async function ncbi(endpoint,params){
 const url=new URL(API+endpoint);for(const [k,v] of Object.entries(params))url.searchParams.set(k,String(v));
 const response=await fetch(url,{mode:'cors',cache:'no-store',credentials:'omit',referrerPolicy:'no-referrer'});
 if(!response.ok)throw Error('NCBI HTTP '+response.status);
 const body=await response.text();if(body.length>8e6)throw Error('Réponse NCBI trop volumineuse');
 if(body.trim().startsWith('<!DOCTYPE html'))throw Error('NCBI a renvoyé une page HTML');return body;
}
function xml(text){const doc=new DOMParser().parseFromString(text,'application/xml');if(doc.querySelector('parsererror'))throw Error('XML NCBI invalide');return doc}
function tagText(node,name){return node?.getElementsByTagName(name)[0]?.textContent?.trim()||''}
function fasta(text){
 const lines=text.replace(/^\uFEFF/,'').trim().split(/\r?\n/);const header=lines[0]||'';
 if(!header.startsWith('>'))throw Error('En-tête FASTA absent');
 const seq=lines.slice(1).join('').replace(/\s/g,'').toUpperCase();
 if(!seq||!/^[ACGTN]+$/.test(seq))throw Error('FASTA invalide : nucléotides A/C/G/T/N attendus');
 const accession=header.match(/^>(NM_\d+\.\d+)\b/)?.[1];
 if(!accession)throw Error('Accession NM_ versionnée absente de l’en-tête FASTA');
 return {seq,header,accession,cdsStart:null,cdsEnd:null};
}
function setReference(ref){
 if(!ref?.seq||!/^NM_\d+\.\d+$/.test(ref.accession))throw Error('Référence NM_ invalide');
 state.ref=ref;state.reads=[];state.rows=[];v5results=[];$('readselect').replaceChildren();$('alignments').replaceChildren();$('variants').replaceChildren();$('v5rows').replaceChildren();$('v5status').textContent='Référence modifiée : relancez « Lire et analyser localement » avant la comparaison.';$('cdsstart').value=ref.cdsStart||'';
 $('refstatus').textContent=`Référence chargée : ${ref.accession} · ${ref.seq.length} nt · CDS ${ref.cdsStart&&ref.cdsEnd?ref.cdsStart+'–'+ref.cdsEnd:'non documentée'} · source ${ref.source||'FASTA manuel'}. Relancez l’analyse ABI pour mettre à jour les alignements.`;
}
$('loadref').onclick=async()=>{try{const file=$('fasta').files[0],text=file?await file.text():$('refpaste').value;setReference(fasta(text))}catch(e){$('refstatus').textContent='Erreur : '+e.message}};
$('searchgene').onclick=async()=>{
 const gene=$('gene').value.trim().toUpperCase(),select=$('transcripts');
 if(!/^[A-Z0-9][A-Z0-9.-]{0,29}$/.test(gene)){$('refstatus').textContent='Symbole de gène invalide';return}
 select.replaceChildren();$('refstatus').textContent='Recherche des transcrits RefSeq de '+gene+'…';
 try{
 const term=`${gene}[Gene] AND Homo sapiens[Organism] AND biomol_mrna[PROP] AND refseq[filter]`;
 const search=JSON.parse(await ncbi('esearch.fcgi',{db:'nuccore',term,retmode:'json',retmax:200})).esearchresult;
 const ids=search?.idlist||[];if(!ids.length)throw Error('Aucun résultat NCBI pour ce gène');
 const summary=JSON.parse(await ncbi('esummary.fcgi',{db:'nuccore',id:ids.join(','),retmode:'json'})).result;
 const found=new Map();for(const id of summary.uids||[]){const item=summary[id],acc=item?.accessionversion||'';
 if(/^NM_\d+\.\d+$/.test(acc)&&!found.has(acc))found.set(acc,item.title||'Transcrit RefSeq');}
 for(const [acc,title] of found){const option=document.createElement('option');option.value=acc;option.textContent=acc+' — '+title.slice(0,100);select.append(option)}
 if(!found.size)throw Error('Aucun NM_ dans les résultats (essayez l’accession directe)');
 $('refstatus').textContent=found.size+' transcrit(s) NM_ trouvé(s) pour '+gene+'. Sélectionnez une version puis cliquez sur Charger.'+(Number(search.count)>200?' Liste limitée aux 200 premiers résultats.':'');
 }catch(e){$('refstatus').textContent='Recherche impossible : '+e.message+'. Essayez une accession directe ou un FASTA manuel.'}
};
async function fetchTranscript(accession){
 if(!/^NM_\d+\.\d+$/.test(accession))throw Error('Saisissez un NM_ avec son numéro de version (ex. NM_000059.4)');
 const request=++refRequest;$('refstatus').textContent='Téléchargement de '+accession+' et de ses annotations CDS…';
 const document=xml(await ncbi('efetch.fcgi',{db:'nuccore',id:accession,rettype:'gb',retmode:'xml'}));
 const record=document.getElementsByTagName('GBSeq')[0];if(!record)throw Error('Enregistrement GenBank absent');
 const returned=tagText(record,'GBSeq_accession-version');if(returned!==accession)throw Error('NCBI a retourné '+returned+' au lieu de '+accession);
 const seq=tagText(record,'GBSeq_sequence').toUpperCase();if(!/^[ACGTN]+$/.test(seq))throw Error('Séquence NCBI invalide');
 let cdsStart=null,cdsEnd=null;
 for(const feature of record.getElementsByTagName('GBFeature')){
 if(tagText(feature,'GBFeature_key')!=='CDS')continue;
 const intervals=Array.from(feature.getElementsByTagName('GBInterval'));
 const starts=intervals.map(i=>Number(tagText(i,'GBInterval_from')||tagText(i,'GBInterval_point')));
 const ends=intervals.map(i=>Number(tagText(i,'GBInterval_to')||tagText(i,'GBInterval_point')));
 if(starts.length&&starts.every(Number.isInteger)&&ends.every(Number.isInteger)){
  const start=Math.min(...starts),end=Math.max(...ends);
  if(start>=1&&end<=seq.length&&seq.slice(start-1,start+2)==='ATG'){cdsStart=start;cdsEnd=end;break}
 }
 }
 if(request!==refRequest)return;
 setReference({seq,header:'>'+returned,accession:returned,cdsStart,cdsEnd,source:'NCBI RefSeq'});
 if(!cdsStart)$('refstatus').textContent+=' Annotation CDS non résolue : aucun HGVS c. calculé.';
}
$('loadtranscript').onclick=()=>fetchTranscript($('transcripts').value).catch(e=>$('refstatus').textContent='Chargement impossible : '+e.message);
$('loadaccession').onclick=()=>fetchTranscript($('accession').value.trim().toUpperCase()).catch(e=>$('refstatus').textContent='Chargement impossible : '+e.message);
// Banded-free Smith-Waterman local alignment with traceback. Limit read/reference size to bound runtime.
function align(read,ref){const n=read.length,m=ref.length;if(n>3000||m>25000)throw Error('Référence/lecture trop longue pour cet alignement de démonstration');
 let w=m+1,H=new Int32Array((n+1)*w),T=new Uint8Array(H.length),best=0,bi=0,bj=0;
 for(let i=1;i<=n;i++){let row=i*w,prev=(i-1)*w;for(let j=1;j<=m;j++){let match=(read[i-1]===ref[j-1]&&'ACGT'.includes(read[i-1]))?3:-3;let a=H[prev+j-1]+match,b=H[prev+j]-4,c=H[row+j-1]-4;let v=Math.max(0,a,b,c);H[row+j]=v;T[row+j]=v===0?0:v===a?1:v===b?2:3;if(v>best){best=v;bi=i;bj=j}}}
 let i=bi,j=bj,cols=[];while(i>0&&j>0&&H[i*w+j]>0){let t=T[i*w+j];if(t===1){cols.push({q:i-1,r:j-1,base:read[i-1],ref:ref[j-1]});i--;j--}else if(t===2){cols.push({q:i-1,r:null,base:read[i-1],ref:'-'});i--}else if(t===3){cols.push({q:null,r:j-1,base:'-',ref:ref[j-1]});j--}else break}cols.reverse();return {score:best,cols,start:i,end:bi,refStart:j,refEnd:bj}}
function chooseAlignment(read,ref){let a=align(read,ref),b=align(rc(read),ref);return b.score>a.score?{...b,strand:'-',oriented:rc(read),length:read.length}:{...a,strand:'+',oriented:read,length:read.length}}
function cpos(refZero,cds){let x=refZero+1;if(!cds||cds<1)return 'non défini';let d=x-cds+1;return d>0?'c.'+d:'c.'+d}
function callRows(read){let a=read.alignment,cds=Number($('cdsstart').value),rows=[];for(let k=0;k<a.cols.length;k++){let col=a.cols[k];if(col.q===null||col.r===null||col.base===col.ref||!/[ACGT]/.test(col.base)||!/[ACGT]/.test(col.ref))continue;
 let qOriginal=a.strand==='+'?col.q:a.length-1-col.q;let label=Number.isInteger(cds)&&cds>0&&state.ref?.cdsStart===cds&&col.r+1>=cds&&col.r+1<=state.ref.cdsEnd?`${state.ref.accession}:${cpos(col.r,cds)}${col.ref}>${col.base} (indicatif)`:'CDS non vérifiée / hors CDS';
 rows.push({name:read.name,q:qOriginal+1,r:col.r+1,ref:col.ref,alt:col.base,hgvs:label,read})}return rows}
$('analyze').onclick=async()=>{state.reads=[];state.rows=[];v5results=[];$('v5rows').replaceChildren();$('v5status').textContent='Analyse ABI en cours…';$('variants').replaceChildren();$('alignments').replaceChildren();$('readselect').replaceChildren();try{
 for(let [id,name] of [['forward','Forward'],['reverse','Reverse']]){let f=$(id).files[0];if(!f)continue;let ab=parseABI(await f.arrayBuffer());let read={...ab,name:name+' · '+f.name};if(state.ref){read.alignment=chooseAlignment(ab.bases,state.ref.seq);state.rows.push(...callRows(read))}state.reads.push(read)}if(!state.reads.length)throw Error('Importez au moins un fichier ABI');
 for(let [i,r] of state.reads.entries()){let opt=document.createElement('option');opt.value=i;opt.textContent=r.name;$('readselect').append(opt);let p=document.createElement('p');p.textContent=`${r.name} — ${r.bases.length} bases${r.alignment?`, alignement ${r.alignment.strand} ; référence ${r.alignment.refStart+1}–${r.alignment.refEnd} ; score ${r.alignment.score}`:''}`;$('alignments').append(p)}
 for(let row of state.rows.slice(0,1000)){let tr=document.createElement('tr');for(let value of [row.name,row.q,row.r,row.ref,row.alt,row.hgvs]){let td=document.createElement('td');td.textContent=value;tr.append(td)}tr.onclick=()=>{let idx=state.reads.indexOf(row.read);$('readselect').value=idx;$('start').value=Math.max(1,row.q-20);draw()};$('variants').append(tr)}
 $('v5status').textContent='ABI chargés : vous pouvez lancer la comparaison exploratoire.'; $('status').textContent=`${state.reads.length} chromatogramme(s) lu(s) ; ${state.rows.length} discordance(s) simples indicatives. ${state.ref?'Référence présente.':'Aucune référence : visualisation uniquement.'}`;draw();
 }catch(e){$('status').textContent='Erreur : '+e.message;$('v5status').textContent='Comparaison indisponible : '+e.message;console.error('SVIB analyse ABI',e)}}
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



// V5: reference-guided two-allele signal fitting. Original implementation, no Tracy code.
// Candidate hypotheses are NOT diagnostic calls or complete de novo haplotype reconstruction.
const DNA='ACGT';let v5results=[];
function intensity(read,k){
 if(k<0||k>=read.positions.length)return null;
 const p=read.positions[k],v=DNA.split('').map(b=>Math.max(0,read.traces[b][p]||0));
 const sum=v.reduce((a,b)=>a+b,0);return sum<25?null:v.map(x=>x/sum);
}
function fit(obs,a,b,f){
 const ia=DNA.indexOf(a),ib=DNA.indexOf(b);if(ia<0||ib<0)return null;
 const pred=[.025,.025,.025,.025];pred[ia]+=.9*(1-f);pred[ib]+=.9*f;
 return pred.reduce((x,v,i)=>x+(v-obs[i])**2,0);
}
// Alignment anchors in original acquisition order, with strand-specific reference bases.
function anchors(read){
 const al=read.alignment,reverse=al.strand==='-';let map=new Map();
 for(const col of al.cols)if(col.q!==null&&col.r!==null)map.set(reverse?read.bases.length-1-col.q:col.q,col.r);
 return map;
}
function orientedBase(ref,r,strand){return r>=0&&r<ref.length?(strand===1?ref[r]:rc(ref[r])):'N'}
function evaluate(read,k,r,kind,len,window=23){
 const ref=state.ref.seq,strand=read.alignment.strand==='+'?1:-1;
 // A positive shift in acquisition direction models deletion; negative models insertion.
 // For an insertion the inserted bases are inferred from post-transition peaks, not assumed from reference.
 const delta=kind==='del'?len:-len;
 let pre0=0,preN=0,post0=0,postN=0,valid=0,best=null;
 for(let d=-12;d<=-3;d++){
  const obs=intensity(read,k+d),wt=orientedBase(ref,r+strand*d,strand);
  if(!obs||wt==='N')continue;pre0+=fit(obs,wt,wt,0);preN++;
 }
 if(preN<7)return null;
 const points=[];
 for(let d=3;d<=window;d++){
  const obs=intensity(read,k+d),wt=orientedBase(ref,r+strand*d,strand);
  let alt=orientedBase(ref,r+strand*(d+delta),strand);
  if(!obs||wt==='N'||alt==='N')continue;
  // Insertions: first inserted segment cannot be reconstructed from reference; omit it here.
  if(kind==='ins'&&d<=len+2)continue;
  points.push({obs,wt,alt});post0+=fit(obs,wt,wt,0);postN++;
 }
 if(postN<9)return null;
 for(let f=.2;f<=.801;f+=.05){
  let loss=0;for(const point of points)loss+=fit(point.obs,point.wt,point.alt,f);
  if(!best||loss<best.loss)best={loss,f};
 }
 const improvement=(post0-best.loss)/postN,preError=pre0/preN,postError=best.loss/postN;
 // Reconstruct a local reference-guided sequence; unknown inserted bases are explicitly N.
 const before=Array.from({length:12},(_,i)=>orientedBase(ref,r+strand*(i-12),strand)).join('');
 const after=Array.from({length:25},(_,i)=>orientedBase(ref,r+strand*(i+delta),strand)).join('');
 const allele1=before+Array.from({length:25},(_,i)=>orientedBase(ref,r+strand*i,strand)).join('');
 const allele2=before+(kind==='ins'?'N'.repeat(len):'')+after;
 return {read,k,r,kind,len,improvement,preError,postError,f:best.f,postN,allele1,allele2};
}
function addV5row(c){
 const tr=document.createElement('tr');
 for(const value of [c.read.name,c.k+1,c.r+1,c.kind==='del'?`Délétion ${c.len} nt`:`Insertion ${c.len} nt`,c.improvement.toFixed(3),c.f.toFixed(2)]){
  const td=document.createElement('td');td.textContent=value;tr.append(td);
 }
 tr.onclick=()=>{ $('readselect').value=state.reads.indexOf(c.read);viewTo(c.k-22,55);
  $('v5alleles').textContent=`Hypothèse ${c.kind} ${c.len} nt · lecture ${c.read.name} · orientation d’acquisition\nAllèle 1 (guidé par référence) : ${c.allele1}\nAllèle 2 (modèle) :            ${c.allele2}\n${c.kind==='ins'?'N = bases insérées inconnues : le modèle ne les détermine pas.\n':''}Coordonnées approximatives ; reconstruction partielle, non confirmée.`;
 };$('v5rows').append(tr);
}
$('v5scan').onclick=()=>{
 const status=$('v5status');$('v5rows').replaceChildren();$('v5alleles').textContent='';v5results=[];
 if(!state.ref||!state.reads.length||state.reads.some(x=>!x.alignment)){status.textContent='Chargez le NM_, puis cliquez sur « Lire et analyser localement ».';return}
 status.textContent='Évaluation des modèles bi-alléliques en cours…';
 // Yield to the UI; all patient-derived data remain in this tab.
 setTimeout(()=>{try{
 for(const read of state.reads){
  const map=anchors(read),strand=read.alignment.strand==='+'?1:-1,ranked=[];
  for(let k=17;k<read.bases.length-35;k+=2){
   const r=map.get(k);if(r===undefined)continue;
   let matched=0;for(let d=-9;d<=-2;d++)if(map.get(k+d)===r+strand*d)matched++;
   if(matched<5)continue;
   for(let len=1;len<=12;len++)for(const kind of ['del','ins']){
    const candidate=evaluate(read,k,r,kind,len);
    if(candidate&&candidate.improvement>0&&candidate.preError<.38)ranked.push(candidate);
   }
  }
  ranked.sort((a,b)=>b.improvement-a.improvement);
  const selected=[];for(const c of ranked){if(selected.every(x=>Math.abs(x.k-c.k)>14))selected.push(c);if(selected.length>=8)break}
  v5results.push(...selected);
 }
 v5results.sort((a,b)=>b.improvement-a.improvement);v5results.forEach(addV5row);
 status.textContent=v5results.length?`${v5results.length} hypothèse(s) de mélange retenue(s) (1–12 nt). Cliquez sur une ligne pour afficher les séquences modélisées. Le score n’est pas une probabilité, et le type/les coordonnées ne sont PAS confirmés.`:'Aucune hypothèse ne présente de gain positif avec les filtres de qualité actuels. Cela n’exclut pas un variant.';
 }catch(e){status.textContent='Erreur du moteur V5 : '+e.message;console.error(e)}},20);
};
$('v5export').onclick=()=>{
 const rows=[['lecture','base_ABI','reference_approximative','hypothese','taille_nt','gain_non_calibre','fraction_modele','allele_reference_modele','allele_alternatif_modele'],...v5results.map(c=>[c.read.name,c.k+1,c.r+1,c.kind,c.len,c.improvement.toFixed(4),c.f.toFixed(2),c.allele1,c.allele2])];
 const csv=rows.map(row=>row.map(v=>'"'+String(v).replace(/"/g,'""')+'"').join(',')).join('\r\n');
 const url=URL.createObjectURL(new Blob(['\ufeff'+csv],{type:'text/csv;charset=utf-8'})),link=document.createElement('a');link.href=url;link.download='SVIB_V5_hypotheses_non_validees.csv';link.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
};
