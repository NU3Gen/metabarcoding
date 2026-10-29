# Nanopore Metabarcoding

Webová a príkazová pipeline na druhovú identifikáciu živočíchov z nezarovnaných
Oxford Nanopore `FASTQ`/`FASTQ.gz` dát. Používateľ nemusí súbory spájať ani
vytvárať ZIP: môže nahrať jeden súbor, viac súborov alebo celý priečinok
`fastq_pass`.

## Aktuálny stav (MVP)

- streamované čítanie FASTQ a FASTQ.gz,
- automatická detekcia primerového páru v oboch orientáciách,
- primer-aware extrakcia insertu,
- dĺžkový a kvalitatívny filter,
- dereplikácia a zhlukovanie podobných insertov,
- porovnanie reprezentantov zhlukov s lokálnou FASTA databázou,
- konzervatívne priradenie druhu/rodu/skupiny,
- JSON, CSV a HTML report,
- konsenzuálne sekvencie vo FASTA a priame nezávislé overenie v NCBI BLAST nt,
- accession, zdroj a auditná úroveň každej referenčnej zhody,
- webový upload viacerých súborov a priečinka,
- validačný profil L15601/H15748 pre približne 148 bp fragment MT-CYB.

Webová aplikácia používa verzovaný panel
`references/releases/nu3gen_animals_v1.fasta`. Súbor
`references/demo_cytb.fasta` zostáva iba validačným minipanelom pre AMP25.
Výsledok z neúplného panelu je `provisional`, nie definitívna identifikácia.

Prvá produkčná databázová vrstva je definovaná v `references/`. Základ tvorí
NCBI RefSeq; MIDORI2 a verejné BOLD záznamy možno pridať iba so zachovaným
accessionom, taxonómiou, markerom, zdrojom a auditnou úrovňou. Builder a zoznam
cieľových hospodárskych a voľne žijúcich druhov sú opísané v
`references/README.md`.

## Najjednoduchší vstup

Používateľ vyberie priamo jeden alebo viac súborov:

```text
*.fastq
*.fastq.gz
```

Pri barcodovanom behu sa súbory každého barcode analyzujú ako samostatná
vzorka. V MVP možno jeden upload označiť jedným názvom vzorky; automatické
rozdelenie podľa barcode je pripravené ako ďalší krok.

## Lokálne spustenie

```bash
python -m venv .venv
. .venv/bin/activate
pip install -e .
uvicorn metabarcoder.webapp:app --host 0.0.0.0 --port 8000
```

Otvorte <http://localhost:8000>.

## CLI

```bash
metabarcoder analyze \
  --input /cesta/k/fastq_pass \
  --reference references/releases/nu3gen_animals_v1.fasta \
  --output results/AMP25 \
  --sample AMP25
```

CLI prijíma súbor, viac `--input` argumentov alebo priečinok. V priečinku
rekurzívne nájde `.fastq` a `.fastq.gz`.

## Testy

Testy jadra nevyžadujú webové závislosti:

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
```

## Vstupný kontrakt

- Sekvencie musia byť nezarovnané reads vo FASTQ/FASTQ.gz.
- FASTQ môže byť priamo z MinKNOW (HAC alebo SUP).
- `fastq_fail` možno pridať, ale predvolene sa odporúča `fastq_pass`.
- POD5 bude podporovaný neskôr ako voliteľný GPU/SUP režim.
- BAM/BAI/BCF/VCF nie sú primárny vstup.

## Primerové profily

Profily sú vo verzovanom súbore
`src/metabarcoder/data/primer_profiles.json`. Každý záznam obsahuje publikovanú
sekvenciu primerov, marker, rozsah insertu, DOI a taxonomický rozsah. Nové
profily sa pridávajú až po kontrole primárneho článku.

## Produkčný plán

1. vytvoriť a validačne uzamknúť prvé vydanie kurátorovanej CYTB/COI/12S/16S databázy,
2. zaradiť vsearch/amplicon_sorter a robustný konsenzus,
3. detekovať chiméry a NUMT, pri kódujúcich markeroch kontrolovať preklad,
4. automaticky rozdeľovať barcode vzorky,
5. pridať negatívne kontroly a pravidlá kontaminácie,
6. spracovanie úloh presunúť do worker queue,
7. voliteľný Dorado SUP basecalling z POD5.

## Licencia

MIT. Biologické výsledky musia byť interpretované v kontexte použitých
primerov, referenčnej databázy, kontrol a validačných prahov; softvér sám osebe
nenahrádza validáciu laboratórnej metódy.

