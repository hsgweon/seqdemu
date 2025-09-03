# 🚀 seqdemu: A No-Nonsense Nanopore Demultiplexer

Tired of sifting through your nanopore sequences manually? Meet **seqdemu**, your no-nonsense, super-efficient demultiplexer for sequences with custom barcodes! 🎉  

---

## ✨ Features

✅ Demultiplex sequences with custom barcodes—on both ends!

✅ Got barcodes in the middle of the sequence? No problem!

✅ Allows for mismatches—because we all make mistakes. 😏

✅ Detects barcode issues (like multiple forward/reverse primers) and automatically excludes problematic sequences.

✅ Supports multiple CPUs for faster processing.

## ⚡ Installation

### Clone the Repository & Set Up the Environment

```bash
# Clone the repo
git clone https://github.com/hsgweon/seqdemu.git
cd seqdemu

# Create the seqdemu environment (ensure that you have conda installed)
mamba create -n seqdemu_env -y -c conda-forge -c bioconda conda-forge::biopython progressbar2
```

## 🎬 Running seqdemu

### 1️⃣ What You Need:
1. A gzipped basecalled FASTQ file (your sequence data).
2. A barcode file (check out the example format below).

Barcode File Format
Your barcodes should be listed like this:

| Forward Primer  | Reverse Primer  | Sample ID |
|----------------|----------------|-----------|
| GCCTCAGGCTTA  | CGCATAAGGCAA  | S01       |
| CTTATGCAATGC  | ATGGACTCGCAA  | S02       |
| TCCACCAGAGGT  | CAATGCGTCCAA  | S03       |
| GGAGAAGAAGAA  | TGTGGTTGCTAA  | S04       |
| AAGCGGAGAGAA  | CCTCATGGAAGA  | S05       |
| CTCTCTCCAGAA  | TGGAGACCAAGA  | S06       |
| TCTCTCCAGGAA  | TACGAGAAGAGA  | S07       |
| GCTACGTTACAA  | GTAGTCCGCAGA  | S08       |
| ATGTCGGTTAGA  | TAGCCTGCGTCA  | S09       |
| CATGTATCTGGA  | CTCCGGTATCTA  | S10       |
| AGACCTAACCGA  | CGAGCGCTGTTA  | S11       |
| GTTCGAATGTGA  | GTCACCGAGAAG  | S12       |
| AGCCTCTTAACA  | AGTTAGCGGAAG  | S13       |
| GCATGTTAGACA  | CCGGTTAACAAG  | S14       |
| ACAATGGCCGCA  | TCACCATGTAAG  | S15       |
| GAGAACCTTGCA  | TAAGTGCCACAG  | S16       |

### 2️⃣ Activate the seqdemu environment

```bash
mamba activate seqdemu_env
```

### 3️⃣ Run seqdemu
Make sure seqdemu is in your PATH, then run:

```bash
seqdemu.py -i GZIPPED_DORADO_FASTQ_FILE -b BARCODE_FILE -m NUMBER_OF_MISMATCH -o OUTPUT_FILENAME -t NUMBER_OF_CPUs -f both
```

## 🧪 Running a Test Example
```bash
cd test
../seqdemu.py -i testseqs.fastq.gz -b barcodes.csv -m 0 -o testseqs -t 10 -f both
```


# 📂 Understanding the Output Files

You'll get different output files, each with a specific purpose:

🟢 _barc → Sequences with barcodes attached.

🟡 _full → Sequences without trimming.

🔵 _noba → Barcodes removed—use this for downstream analysis!

---

That's it! 🚀 Now go forth and demultiplex with confidence! 🎉

