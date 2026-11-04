"use strict";

const PROFILES=[{id:"cytb_l15601_h15748",name:"CYTB L15601/H15748",marker:"MT-CYB",forward:"TACGCAATCCTACGATCAATTCC",reverse:"GGTTGTCCTCCAATTCATGTTAG",minInsert:125,maxInsert:175,maxErrors:3}];
const REFERENCE_URL="https://raw.githubusercontent.com/NU3Gen/metabarcoding/main/references/releases/nu3gen_animals_v1.fasta";
const COMPLEMENT={A:"T",C:"G",G:"C",T:"A",R:"Y",Y:"R",M:"K",K:"M",B:"V",V:"B",D:"H",H:"D",N:"N"};
let lastResult=null;

const $=selector=>document.querySelector(selector);
const esc=value=>String(value).replace(/[&<>"']/g,c=>({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[c]));
const reverseComplement=seq=>[...seq.toUpperCase()].reverse().map(x=>COMPLEMENT[x]||"N").join("");
const kmerSet=(sequence,k=7)=>{const set=new Set();for(let i=0;i<=sequence.length-k;i++)set.add(sequence.slice(i,i+k));return set;};

function levenshtein(a,b,limit=Infinity){
  if(a===b)return 0;if(!a)return b.length;if(!b)return a.length;
  if(a.length>b.length)[a,b]=[b,a];
  let prev=Array.from({length:a.length+1},(_,i)=>i);
  for(let i=1;i<=b.length;i++){
    const cur=[i];let rowMin=i;
    for(let j=1;j<=a.length;j++){
      cur[j]=Math.min(cur[j-1]+1,prev[j]+1,prev[j-1]+(a[j-1]===b[i-1]?0:1));rowMin=Math.min(rowMin,cur[j]);
    }
    if(rowMin>limit)return limit+1;prev=cur;
  }
  return prev[a.length];
}

function approximateIndex(seq,pattern,maxErrors,start=0,end=seq.length){
  const exact=seq.indexOf(pattern,start);if(exact>=0&&exact<end)return {start:exact,end:exact+pattern.length,errors:0};
  const candidates=new Set(),parts=maxErrors+1,seedLength=Math.floor(pattern.length/parts);
  for(let part=0;part<parts;part++){
    const offset=part*seedLength,seed=pattern.slice(offset,part===parts-1?pattern.length:offset+seedLength);let found=seq.indexOf(seed,start);
    while(found>=0&&found<end){for(let shift=-maxErrors;shift<=maxErrors;shift++)candidates.add(found-offset+shift);found=seq.indexOf(seed,found+1);}
  }
  let best=null;
  for(const i of candidates){if(i<start||i>=end)continue;for(let delta=-maxErrors;delta<=maxErrors;delta++){
    const length=pattern.length+delta;if(length<1||i+length>seq.length)continue;const errors=levenshtein(pattern,seq.slice(i,i+length),maxErrors);
    if(errors<=maxErrors&&(!best||errors<best.errors))best={start:i,end:i+length,errors};
  }
  }
  return best;
}

function findAmplicon(sequence,profile){
  const binding=reverseComplement(profile.reverse);
  for(const oriented of [sequence.toUpperCase(),reverseComplement(sequence)]){
    const left=approximateIndex(oriented,profile.forward,profile.maxErrors,0,Math.min(oriented.length,500));
    if(!left)continue;
    const from=left.end+profile.minInsert,to=Math.min(oriented.length,left.end+profile.maxInsert+binding.length+profile.maxErrors);
    const right=approximateIndex(oriented,binding,profile.maxErrors,from,to+1);
    if(right){const insert=oriented.slice(left.end,right.start);if(insert.length>=profile.minInsert&&insert.length<=profile.maxInsert)return insert;}
  }
  return null;
}

function parseFasta(text){
  const records=[];let header="",parts=[];
  const emit=()=>{if(!header)return;const fields=Object.fromEntries(header.split("|").filter(x=>x.includes("=")).map(x=>x.split(/=(.*)/s).slice(0,2))),sequence=parts.join("").toUpperCase();records.push({accession:fields.accession||header.split("|")[0],taxon:fields.taxon||"unclassified",source:fields.source||"unknown",audit:fields.audit||"unreviewed",sequence,kmers:kmerSet(sequence)});};
  for(const raw of text.split(/\r?\n/)){const line=raw.trim();if(!line)continue;if(line.startsWith(">")){emit();header=line.slice(1);parts=[];}else parts.push(line);}emit();return records;
}

async function readFile(file){
  if(file.name.toLowerCase().endsWith(".gz")){
    if(!window.DecompressionStream)throw new Error("Tento prehliadač nepodporuje rozbalenie gzip. Použite Chrome/Edge alebo nahrajte FASTQ.");
    return new Response(file.stream().pipeThrough(new DecompressionStream("gzip"))).text();
  }
  return file.text();
}

function classify(sequence,references){
  const queryKmers=kmerSet(sequence),shortlist=references.map(ref=>{let shared=0;for(const kmer of queryKmers)if(ref.kmers.has(kmer))shared++;return {ref,shared};}).sort((a,b)=>b.shared-a.shared).slice(0,6).map(x=>x.ref);
  const scores=shortlist.map(ref=>{const maxLength=Math.max(sequence.length,ref.sequence.length),distance=sequence.length===ref.sequence.length?[...sequence].reduce((n,c,i)=>n+(c!==ref.sequence[i]),0):levenshtein(sequence,ref.sequence,Math.ceil(maxLength*.20));return {ref,identity:1-distance/maxLength};}).sort((a,b)=>b.identity-a.identity);
  const best=scores[0],second=scores[1];if(!best||best.identity<.90)return {category:"unassigned",taxon:"Nepriradené",identity:best?.identity||0};
  const margin=second?best.identity-second.identity:1;
  let category="low-confidence";
  if(margin<.01)category="ambiguous";else if(best.identity>=.97)category="assigned";else if(best.identity>=.95)category="probable";
  return {category,taxon:best.ref.taxon,accession:best.ref.accession,source:best.ref.source,audit:best.ref.audit,identity:best.identity,sequence};
}

function buildKmerIndex(references){
  const index=new Map();
  references.forEach((ref,refIndex)=>{for(const kmer of ref.kmers){const matches=index.get(kmer);if(matches)matches.push(refIndex);else index.set(kmer,[refIndex]);}});
  return index;
}

function recoverLowAbundance(sequence,references,index){
  const queryKmers=kmerSet(sequence);if(!queryKmers.size)return {category:"unassigned"};
  const shared=new Uint16Array(references.length);
  for(const kmer of queryKmers){const matches=index.get(kmer);if(matches)for(const refIndex of matches)shared[refIndex]++;}
  const byTaxon=new Map();
  references.forEach((ref,i)=>{const denominator=Math.min(queryKmers.size,ref.kmers.size),score=denominator?shared[i]/denominator:0,old=byTaxon.get(ref.taxon);if(!old||score>old.score)byTaxon.set(ref.taxon,{ref,score});});
  const ranked=[...byTaxon.values()].sort((a,b)=>b.score-a.score),best=ranked[0],second=ranked[1];
  if(!best||best.score<.67)return {category:"unassigned"};
  if(second&&best.score-second.score<.03)return {category:"ambiguous"};
  return {category:"low-confidence",taxon:best.ref.taxon,ref:best.ref,score:best.score};
}

function parseFastq(text){
  const lines=text.split(/\r?\n/),records=[];
  for(let i=0;i<lines.length;){if(!lines[i]){i++;continue;}if(i+3>=lines.length||!lines[i].startsWith("@")||!lines[i+2].startsWith("+"))throw new Error(`Neplatný alebo neúplný FASTQ záznam pri riadku ${i+1}.`);const sequence=lines[i+1].trim().toUpperCase(),quality=lines[i+3].trim();if(!/^[ACGTRYSWKMBDHVN]+$/i.test(sequence))throw new Error(`Sekvencia obsahuje neplatné znaky pri riadku ${i+2}.`);if(sequence.length!==quality.length)throw new Error(`Sekvencia a kvalita majú inú dĺžku pri riadku ${i+1}.`);const meanQ=quality.length?[...quality].reduce((n,c)=>n+c.charCodeAt(0)-33,0)/quality.length:0;records.push({sequence,meanQ});i+=4;}
  return records;
}

async function detectProfile(records){
  const sample=records.length<=1000?records:records.filter((_,i)=>i%Math.ceil(records.length/1000)===0).slice(0,1000),scores=[];
  for(const profile of PROFILES){let hits=0;for(let i=0;i<sample.length;i++){if(findAmplicon(sample[i].sequence,profile))hits++;if(i%100===0)await new Promise(r=>setTimeout(r,0));}scores.push({profile,hits,sampled:sample.length,support:sample.length?hits/sample.length:0});}
  return scores.sort((a,b)=>b.hits-a.hits)[0];
}

async function analyzeTexts(sample,texts){
  setStatus("Načítavam referenčný panel…",10);
  const refResponse=await fetch(REFERENCE_URL);if(!refResponse.ok)throw new Error("Referenčný panel sa nepodarilo načítať.");
  const references=parseFasta(await refResponse.text());
  const records=texts.flatMap(parseFastq);if(!records.length)throw new Error("FASTQ neobsahuje žiadne úplné záznamy.");
  const meanQuality=records.reduce((n,r)=>n+r.meanQ,0)/records.length,meanLength=records.reduce((n,r)=>n+r.sequence.length,0)/records.length;
  setStatus(`Kontrolujem primerový profil na vzorke čítaní…`,15);const detection=await detectProfile(records),profile=detection?.profile;
  const counts={assigned:0,probable:0,"low-confidence":0,ambiguous:0,unassigned:0,"insufficient length":0},variants=new Map();
  if(!detection||detection.hits===0){counts.unassigned=records.length;const result={schema_version:"prototype-0.3",sample,created_at:new Date().toISOString(),marker:null,primer_profile:null,primer_detection:{status:"unsupported",sampled:detection?.sampled||0,hits:0,support:0},total_reads:records.length,assigned_reads:0,mean_read_length:meanLength,mean_read_quality:meanQuality,categories:counts,mixed_sample:false,taxa:[],dominant_consensus:"",warnings:["FASTQ je platný, ale nenašiel sa žiadny podporovaný primerový pár.","Taxonomická identifikácia sa nevykonala, aby nevznikol falošný výsledok."]};setStatus("FASTQ bol načítaný, marker však zatiaľ nie je podporovaný.",100);return result;}
  setStatus(`Rozpoznaný ${profile.name}; analyzujem ${records.length.toLocaleString("sk-SK")} čítaní…`,25);
  for(let i=0;i<records.length;i++){
    const seq=records[i].sequence;
    if(seq.length<profile.forward.length+profile.reverse.length+profile.minInsert){counts["insufficient length"]++;continue;}
    const insert=findAmplicon(seq,profile);if(!insert){counts.unassigned++;continue;}
    variants.set(insert,(variants.get(insert)||0)+1);
    if(i%250===0){setStatus(`Analyzujem čítanie ${i+1} z ${records.length}…`,25+55*i/Math.max(1,records.length));await new Promise(r=>setTimeout(r,0));}
  }
  const hits=[],rankedVariants=[...variants].sort((a,b)=>b[1]-a[1]),selectedVariants=rankedVariants.slice(0,250),deferredVariants=rankedVariants.slice(250),deferredReads=deferredVariants.reduce((n,row)=>n+row[1],0);
  for(const [sequence,reads] of selectedVariants){const hit=classify(sequence,references);counts[hit.category]+=reads;if(hit.category!=="unassigned")hits.push({...hit,reads});}
  const grouped=new Map();for(const hit of hits){const key=`${hit.taxon}|${hit.accession}`;const old=grouped.get(key);if(!old||hit.identity>old.identity)grouped.set(key,{...hit,reads:(old?.reads||0)+hit.reads});else old.reads+=hit.reads;}
  setStatus(`Dopĺňam nízko početné varianty…`,82);
  const kmerIndex=buildKmerIndex(references),recovered=new Map();
  for(let i=0;i<deferredVariants.length;i++){
    const [sequence,reads]=deferredVariants[i],match=recoverLowAbundance(sequence,references,kmerIndex);counts[match.category]+=reads;
    if(match.taxon){const old=recovered.get(match.taxon);if(!old)recovered.set(match.taxon,{reads,representative:sequence,representativeReads:reads,ref:match.ref});else{old.reads+=reads;if(reads>old.representativeReads){old.representative=sequence;old.representativeReads=reads;old.ref=match.ref;}}}
    if(i%1000===0){setStatus(`Dopĺňam variant ${i+1} z ${deferredVariants.length}…`,82+13*i/Math.max(1,deferredVariants.length));await new Promise(r=>setTimeout(r,0));}
  }
  for(const [taxon,recovery] of recovered){
    const anchor=classify(recovery.representative,references);if(anchor.taxon!==taxon)continue;
    const key=`${anchor.taxon}|${anchor.accession}`,old=grouped.get(key);if(old)old.reads+=recovery.reads;else grouped.set(key,{...anchor,category:"low-confidence",reads:recovery.reads});
  }
  const taxa=[...grouped.values()].sort((a,b)=>b.reads-a.reads),assignedReads=counts.assigned+counts.probable+counts["low-confidence"]+counts.ambiguous;
  taxa.forEach(x=>x.fraction=assignedReads?x.reads/assignedReads:0);
  const mixed=taxa.filter(x=>x.fraction>=.05).length>1;
  const result={schema_version:"prototype-0.4",sample,created_at:new Date().toISOString(),marker:profile.marker,primer_profile:profile.id,primer_detection:{status:"detected",sampled:detection.sampled,hits:detection.hits,support:detection.support},total_reads:records.length,assigned_reads:assignedReads,mean_read_length:meanLength,mean_read_quality:meanQuality,unique_inserts:variants.size,classified_variants:selectedVariants.length,recovered_low_abundance_reads:[...recovered.values()].reduce((n,row)=>n+row.reads,0),deferred_low_abundance_reads:deferredReads,categories:counts,mixed_sample:mixed,taxa,dominant_consensus:taxa[0]?.sequence||"",warnings:["Výsledok je platný iba v rozsahu lokálneho referenčného panelu.","Najpočetnejšie varianty sa klasifikujú presnou podobnosťou; nízko početné varianty sa dopĺňajú konzervatívnym k-mer porovnaním a uvádzajú sa ako low-confidence.","Toto je validačný prototyp, nie akreditovaná diagnostická metóda."]};
  setStatus("Analýza je dokončená.",100);return result;
}

function render(result){
  lastResult=result;$("#results").hidden=false;$("#result-title").textContent=result.sample;$("#mixed-badge").textContent=result.mixed_sample?"Zmiešaná vzorka":"Jednodruhová / bez preukázanej zmesi";
  const metrics=[[result.total_reads,"Celkové čítania"],[result.assigned_reads,"Priradené čítania"],[result.mean_read_quality?.toFixed(1)||"—","Priemerná kvalita Q"],[result.taxa[0]?`${(100*result.taxa[0].identity).toFixed(1)} %`:"—","Najvyššia identita"]];
  $("#metrics").innerHTML=metrics.map(([v,l])=>`<div class="metric"><strong>${esc(v)}</strong><span>${esc(l)}</span></div>`).join("");
  const labels={assigned:"Assigned",probable:"Probable","low-confidence":"Low confidence",ambiguous:"Ambiguous",unassigned:"Unassigned","insufficient length":"Insufficient length"};
  $("#category-chart").innerHTML=Object.entries(result.categories).map(([k,v])=>`<div class="bar-row"><span>${labels[k]}</span><div class="bar-track"><div class="bar" style="width:${100*v/Math.max(1,result.total_reads)}%"></div></div><b>${v}</b></div>`).join("");
  $("#composition").innerHTML=result.taxa.length?result.taxa.map(x=>`<div class="species-row"><div class="species-label"><i>${esc(x.taxon)}</i><b>${(100*x.fraction).toFixed(1)} %</b></div><div class="bar-track"><div class="bar" style="width:${100*x.fraction}%"></div></div></div>`).join(""):"<p class=empty>Nenašla sa klasifikovateľná zhoda.</p>";
  $("#taxa").innerHTML=result.taxa.map(x=>`<tr><td><i>${esc(x.taxon)}</i></td><td>${x.reads}</td><td>${(100*x.fraction).toFixed(1)} %</td><td>${esc(x.accession)}</td><td>${(100*x.identity).toFixed(2)} %</td><td>${esc(x.source)}</td><td>${esc(x.audit)}</td></tr>`).join("");
  const dominant=result.taxa[0];$("#summary").textContent=dominant?`Automaticky rozpoznaný marker ${result.marker} (${(100*result.primer_detection.support).toFixed(1)} % podpory v kontrolnej vzorke readov). Vo vzorke dominuje ${dominant.taxon} (${(100*dominant.fraction).toFixed(1)} % priradených čítaní; zhoda s ${dominant.accession}: ${(100*dominant.identity).toFixed(2)} %). ${result.mixed_sample?"Vzorka spĺňa prototypové kritérium zmesi.":"Významná minoritná zložka nebola preukázaná."}`:result.primer_detection?.status==="unsupported"?"Súbor je platný FASTQ, ale neobsahuje zatiaľ podporovaný primerový pár. Taxonomická identifikácia nebola vykonaná.":"V dostupnom referenčnom paneli sa nenašla dostatočne podobná sekvencia.";
  $("#results").scrollIntoView({behavior:"smooth"});
}

function download(name,content,type){const url=URL.createObjectURL(new Blob([content],{type}));const a=document.createElement("a");a.href=url;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);}
function exportResult(kind){if(!lastResult)return;const safe=lastResult.sample.replace(/[^a-z0-9_-]+/gi,"_");if(kind==="json")download(`${safe}.json`,JSON.stringify(lastResult,null,2),"application/json");if(kind==="fasta")download(`${safe}_consensus.fasta`,lastResult.taxa.map((x,i)=>`>${safe}|cluster_${i+1}|taxon=${x.taxon}|reads=${x.reads}\n${x.sequence}`).join("\n"),"text/plain");if(kind==="csv"){const head="sample,total_reads,assigned_reads,assigned,probable,low_confidence,ambiguous,unassigned,insufficient_length,mixed_sample,taxon,reads,fraction,accession,identity,source,audit";const rows=(lastResult.taxa.length?lastResult.taxa:[{}]).map(x=>[lastResult.sample,lastResult.total_reads,lastResult.assigned_reads,lastResult.categories.assigned,lastResult.categories.probable,lastResult.categories["low-confidence"],lastResult.categories.ambiguous,lastResult.categories.unassigned,lastResult.categories["insufficient length"],lastResult.mixed_sample,x.taxon||"",x.reads||0,x.fraction||0,x.accession||"",x.identity||0,x.source||"",x.audit||""].map(v=>`"${String(v).replaceAll('"','""')}"`).join(","));download(`${safe}.csv`,[head,...rows].join("\n"),"text/csv");}}
function setStatus(text,value){$("#status").textContent=text;$("#progress").hidden=false;$("#progress").value=value;}

$("#files").addEventListener("change",e=>{$("#analyze").disabled=!e.target.files.length;$("#file-summary").textContent=e.target.files.length?`${e.target.files.length} súbor(ov) · ${(Array.from(e.target.files).reduce((n,f)=>n+f.size,0)/1048576).toFixed(1)} MB`:"Jeden súbor alebo viac súborov naraz · bez ZIP";});
$("#analyze").addEventListener("click",async()=>{try{$("#analyze").disabled=true;const texts=[];for(const file of $("#files").files)texts.push(await readFile(file));render(await analyzeTexts($("#sample").value.trim()||"sample",texts));}catch(e){setStatus(`Chyba: ${e.message}`,0);}finally{$("#analyze").disabled=false;}});
$("#demo").addEventListener("click",async()=>{const profile=PROFILES[0],insert="TAACAAACTAGGAGGAGTACTAGCCCTAATCTCATCAATCCTAATCTTGATCCTTATACCCCTCCTCCATACATCTAAACAACGCAGTATAATGTTCCGGCCATTTAGTCAATGCCTATTCTGAATCCTAGTAGCTGACCTATTAACA";const read=profile.forward+insert+reverseComplement(profile.reverse);const text=Array.from({length:24},(_,i)=>`@demo_${i+1}\n${read}\n+\n${"I".repeat(read.length)}`).join("\n")+"\n";try{render(await analyzeTexts("demo_roedeer",[text]));}catch(e){setStatus(`Chyba: ${e.message}`,0);}});
document.querySelectorAll("[data-download]").forEach(b=>b.addEventListener("click",()=>exportResult(b.dataset.download)));
$("#blast").addEventListener("click",()=>{if(!lastResult?.dominant_consensus)return;const form=document.createElement("form");form.method="get";form.action="https://blast.ncbi.nlm.nih.gov/Blast.cgi";form.target="_blank";for(const [name,value] of Object.entries({PROGRAM:"blastn",DATABASE:"nt",PAGE_TYPE:"BlastSearch",QUERY:lastResult.dominant_consensus})){const input=document.createElement("input");input.type="hidden";input.name=name;input.value=value;form.append(input);}document.body.append(form);form.submit();form.remove();});
