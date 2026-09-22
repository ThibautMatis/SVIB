'use strict';
// ABIF v1 parser: no network, no external dependencies. Untrusted inputs are bounds-checked.
const $=id=>document.getElementById(id);
const state={ref:null,reads:[],rows:[],view:null};
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
function draw(){let r=state.reads[Number($('readselect').value)||0];if(!r)return;let start=Math.max(0,(Number($('start').value)||1)-1),count=Math.min(120,Math.max(10,Number($('count').value)||45));let end=Math.min(r.bases.length,start+count);if(end<=start)return;
 let canvas=$('trace'),ctx=canvas.getContext('2d'),W=canvas.width,H=canvas.height;ctx.clearRect(0,0,W,H);let x0=Math.max(0,r.positions[start]-12),x1=Math.min(r.traces.A.length-1,r.positions[end-1]+12),span=Math.max(1,x1-x0);let ymax=1;for(let b of 'ACGT'){let t=r.traces[b];for(let i=x0;i<=x1;i++)if(t[i]>ymax)ymax=t[i]};let color={A:'#16854e',C:'#2166cc',G:'#252b32',T:'#d33a37'};ctx.lineWidth=1.5;
 for(let b of 'ACGT'){ctx.beginPath();ctx.strokeStyle=color[b];let t=r.traces[b];for(let i=x0;i<=x1;i++){let x=(i-x0)/span*W,y=H-45-t[i]/ymax*(H-80);if(i===x0)ctx.moveTo(x,y);else ctx.lineTo(x,y)}ctx.stroke()}
 ctx.font='13px Segoe UI';ctx.textAlign='center';for(let k=start;k<end;k++){let x=(r.positions[k]-x0)/span*W;ctx.fillStyle=color[r.bases[k]]||'#8c60a7';ctx.fillText(r.bases[k],x,H-22);if(k%5===0){ctx.fillStyle='#647587';ctx.font='10px Segoe UI';ctx.fillText(String(k+1),x,H-6);ctx.font='13px Segoe UI'}}state.view={r,x0,x1,start,end};$('baseinfo').textContent=`${r.name} · bases ${start+1}–${end} · cliquer sur un pic pour afficher sa position.`}
$('draw').onclick=draw;$('readselect').onchange=draw;$('trace').onclick=e=>{let v=state.view;if(!v)return;let rect=$('trace').getBoundingClientRect(),sample=v.x0+(e.clientX-rect.left)/rect.width*(v.x1-v.x0),k=v.start;for(let i=v.start;i<v.end;i++)if(Math.abs(v.r.positions[i]-sample)<Math.abs(v.r.positions[k]-sample))k=i;let a=v.r.alignment,ref='non alignée';if(a){let q=a.strand==='+'?k:a.length-1-k;let c=a.cols.find(c=>c.q===q&&c.r!==null);if(c)ref=`référence ${c.r+1} ; ${cpos(c.r,Number($('cdsstart').value))}`} $('baseinfo').textContent=`${v.r.name} · base ABI ${k+1} : ${v.r.bases[k]} · ${ref}`};
$('export').onclick=()=>{let header=['lecture','position_ABI','position_reference_1based','reference','observe','HGVS_indicatif'];let esc=v=>'"'+String(v).replace(/"/g,'""')+'"';let data=[header,...state.rows.map(r=>[r.name,r.q,r.r,r.ref,r.alt,r.hgvs])].map(row=>row.map(esc).join(',')).join('\r\n');let blob=new Blob(['\ufeff'+data],{type:'text/csv;charset=utf-8'}),url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;a.download='sangervariant_discordances_indicatives.csv';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000)};
