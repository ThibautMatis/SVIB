/* SVIB Tracy WASM adapter. NO remote uploads; fail closed if WASM is missing. */
(()=>{'use strict';
const $=id=>document.getElementById(id), status=$('tracystatus'), run=$('tracyrun');
let createTracy=null;
const names=['tracy/tracy.js','tracy/tracy.wasm'];
async function probe(){
 run.disabled=true;
 try{
  for(const name of names){const r=await fetch(name,{cache:'no-store'});if(!r.ok)throw Error(name+' absent (HTTP '+r.status+')');
   const type=r.headers.get('content-type')||'';if(name.endsWith('.wasm')&& !/wasm|octet-stream/i.test(type)) status.textContent='Attention : type MIME WASM inattendu : '+type;
  }
  if(!createTracy){const mod=await import('./tracy/tracy.js');createTracy=mod.default;
    if(typeof createTracy!=='function')throw Error('tracy.js doit exporter une factory Emscripten (MODULARIZE/EXPORT_ES6)');}
  status.textContent='Fichiers WASM détectés. Initialisation de Tracy à vérifier à l’exécution.';
  run.disabled=false;
 }catch(e){status.textContent='Tracy non disponible : '+e.message+' — aucune analyse Tracy effectuée.';}
}
$('tracycheck').addEventListener('click',probe);
run.addEventListener('click',async()=>{
 run.disabled=true;$('tracyresults').textContent='';
 try{
  if(!state.ref?.seq)throw Error('Chargez une référence NM_ avant de lancer Tracy.');
  const files=[['Forward',$('forward').files[0]],['Reverse',$('reverse').files[0]]].filter(x=>x[1]);
  if(!files.length)throw Error('Importez au moins un ABI.');
  // A fresh module per run avoids leftover virtual files from a previous patient.
  const Module=await createTracy({locateFile:(file)=>'tracy/'+file,print:()=>{},printErr:(s)=>{status.textContent='Tracy : '+s}});
  if(!Module.FS||!Module.callMain)throw Error('Compilation incompatible : FS et callMain requis (FORCE_FILESYSTEM, INVOKE_RUN=0, EXPORTED_RUNTIME_METHODS=FS,callMain).');
  Module.FS.writeFile('/reference.fa','>reference\n'+state.ref.seq+'\n');
  const report=[];
  for(const [label,file] of files){
    const input='/'+label+'.ab1',prefix='/'+label.toLowerCase();
    Module.FS.writeFile(input,new Uint8Array(await file.arrayBuffer()));
    let exitCode;
    try{exitCode=Module.callMain(['decompose','-r','/reference.fa','-o',prefix,input]);}
    catch(e){throw Error(label+' : échec de Tracy : '+e.message)}
    const found=Module.FS.readdir('/').filter(n=>n.startsWith(label.toLowerCase()+'.'));
    report.push(label+' : code '+String(exitCode??'non communiqué')+' ; sorties : '+found.join(', '));
    for(const suffix of ['.align1','.align2']){
      const path=prefix+suffix;
      try{report.push('\n'+label+' '+suffix+'\n'+Module.FS.readFile(path,{encoding:'utf8'}));}
      catch{report.push(label+' '+suffix+' absent');}
    }
  }
  $('tracyresults').textContent=report.join('\n');
  status.textContent='Tracy exécuté localement. Alignements bruts à vérifier ; aucun HGVS ni concordance F/R automatique.';
 }catch(e){status.textContent='Analyse Tracy impossible : '+e.message;}
 finally{run.disabled=false;}
});
})();
