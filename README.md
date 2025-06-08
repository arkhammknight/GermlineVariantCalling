# Variant Calling Pipeline for *Mycobacterium ulcerans* CSURP7741

This repository contains a pipeline for performing variant calling on whole-genome sequencing (WGS) data of *Mycobacterium ulcerans* CSURP7741, a clinical isolate from French Guiana. The pipeline processes paired-end FASTQ files, aligns them to a reference genome, and identifies single nucleotide polymorphisms (SNPs) and insertions/deletions (indels) using `bcftools`. The pipeline uses `conda` for environment management and tools like `bwa`, `samtools`, and `sickle`.

## Prerequisites

- **Software**:
  - `conda` (with `conda-forge` and `bioconda` channels)
  - Tools installed in a `conda` environment: `bwa`, `samtools`, `sickle`, `python`, `biopython`
- **Input Data**:
  - Paired-end FASTQ files: `P7741_R1.fastq.gz` and `P7741_R2.fastq.gz` (downloaded from [ENA: ERR3335404](https://trace.ncbi.nlm.nih.gov/Traces/index.html?view=run_browser&acc=ERR3335404&display=data-access))
  - Reference genome: `sequence.fasta` (*Mycobacterium ulcerans* Agy99, [NCBI: CP000325.1](https://www.ncbi.nlm.nih.gov/nuccore/CP000325.1))
- **Hardware**:
  - Linux-based system with sufficient storage and memory (e.g., ~300 MB for SAM files, ~90 MB for BAM files)
  - Multi-core CPU for parallel processing (e.g., 8 threads used in this pipeline)

## Setup

1. **Configure `conda` channels**:
   ```bash
   conda config --add channels conda-forge
   conda config --add channels bioconda
   ```

2. **Create and activate a `conda` environment**:
   ```bash
   conda create -n variants bwa bcftools samtools sickle python biopython
   source /path/to/miniconda3/etc/profile.d/conda.sh
   conda activate variants
   ```

3. **Download input data**:
   ```bash
   wget http://ftp.sra.ebi.ac.uk/vol1/run/ERR333/ERR3335404/P7741_R1.fastq.gz
   wget http://ftp.sra.ebi.ac.uk/vol1/run/ERR333/ERR3335404/P7741_R2.fastq.gz
   ```

4. **Prepare reference genome**:
   - Download the reference genome (`sequence.fasta`) from [NCBI: CP000325.1](https://www.ncbi.nlm.nih.gov/nuccore/CP000325.1).
   - Create a `ref` directory and move the reference genome there:
     ```bash
     mkdir ref
     mv sequence.fasta ref/
     ```

## Pipeline Steps

### 1. Quality Trimming with `sickle`
Trim low-quality reads from the paired-end FASTQ files to improve alignment accuracy.

```bash
sickle pe -f P7741_R1.fastq.gz -r P7741_R2.fastq.gz -t sanger -q 20 -l 20 -g -o trimmed_R1.fastq.gz -p trimmed_R2.fastq.gz -s trimmed_S.fast.qz
```

- **Explanation**:
  - `-f`, `-r`: Input paired-end FASTQ files (forward and reverse reads).
  - `-t sanger`: Specifies Sanger quality encoding for FASTQ files.
  - `-q 20`: Minimum quality score threshold (Phred score ≥ 20).
  - `-l 20`: Minimum read length after trimming (≥ 20 bp).
  - `-g`: Output compressed (gzip) FASTQ files.
  - `-o`, `-p`: Output trimmed paired-end files (`trimmed_R1.fastq.gz`, `trimmed_R2.fastq.gz`).
  - `-s`: Output singleton reads (unpaired reads after trimming) to `trimmed_S.fast.qz`.
- **Output**:
  - Total input: 546,142 reads (273,071 pairs).
  - Kept: 538,788 paired reads (269,394 pairs), 3,188 singletons.
  - Discarded: 978 paired reads (489 pairs), 3,188 singletons.
  - Files: `trimmed_R1.fastq.gz`, `trimmed_R2.fastq.gz`, `trimmed_S.fast.qz`.

### 2. Index the Reference Genome
Create index files for the reference genome using `bwa`.

```bash
bwa index ref/sequence.fasta
```

- **Explanation**:
  - Generates index files (`.amb`, `.ann`, `.bwt`, `.pac`, `.sa`) required for alignment.
- **Output**:
  - Files in `ref/`: `sequence.fasta.amb`, `sequence.fasta.ann`, `sequence.fasta.bwt`, `sequence.fasta.pac`, `sequence.fasta.sa`.

### 3. Align Reads to Reference
Align trimmed FASTQ files to the reference genome using `bwa mem` and sort the output using `samtools`.

```bash
bwa mem -t 8 ref/sequence.fasta trimmed_R1.fastq.gz trimmed_R2.fastq.gz | samtools sort -o output.sorted.bam -
```

- **Explanation**:
  - `bwa mem`: Aligns paired-end reads using the BWA-MEM algorithm.
    - `-t 8`: Uses 8 threads for parallel processing.
    - Inputs: Reference genome (`sequence.fasta`), trimmed FASTQ files.
  - `| samtools sort`: Pipes the SAM output to `samtools sort` to produce a sorted BAM file.
  - **Note**: An initial attempt failed due to a typo (`trimmed_R1.fastq.qz` instead of `.gz`). Corrected by using the proper file names.
- **Output**:
  - Sorted BAM file: `output.sorted.bam` (~91 MB).
  - Alignment stats (from `samtools flagstat`):
    - Total reads: 566,095 (538,788 primary, 27,307 supplementary).
    - Mapped reads: 502,161 (88.71%), with 474,854 primary mapped (88.13%).
    - Properly paired: 444,342 (82.47%).
    - Singletons: 4,136 (0.77%).
  - Stats saved to `mapping.stats.txt`.

### 4. Variant Calling with `bcftools`
Generate a pileup and call variants using `bcftools`.

```bash
bcftools mpileup -O b -o raw.bcf -f ref/sequence.fasta --threads 8 -q 20 -Q 30 output.sorted.bam
bcftools call --ploidy 1 -m -v -o variants.raw.vcf raw.bcf
```

- **Explanation**:
  - `bcftools mpileup`:
    - `-O b`: Output in BCF format (`raw.bcf`).
    - `-f`: Reference genome FASTA file.
    - `--threads 8`: Uses 8 threads.
    - `-q 20`: Minimum mapping quality (Phred score ≥ 20).
    - `-Q 30`: Minimum base quality (Phred score ≥ 30).
    - `-d 250`: Maximum read depth per input file (default setting observed).
  - `bcftools call`:
    - `--ploidy 1`: Assumes haploid genome (appropriate for *Mycobacterium ulcerans*).
    - `-m`: Multiallelic caller mode.
    - `-v`: Outputs only variant sites.
    - Outputs VCF file: `variants.raw.vcf`.
- **Output**:
  - BCF file: `raw.bcf`.
  - VCF file: `variants.raw.vcf`.

### 5. Analyze Variants
Extract and count variants (SNPs and indels) from the VCF file.

```bash
# Total variants (excluding header lines)
grep -v -c '^#' variants.raw.vcf
# Result: 32,665 variants

# SNPs only
bcftools view -v snps variants.raw.vcf | grep -v -c '^#'
# Result: 31,018 SNPs

# Indels only
bcftools view -v indels variants.raw.vcf | grep -v -c '^#'
# Result: 1,647 indels

# Extract variant positions
bcftools query -f '%POS\n' variants.raw.vcf > pos.txt
```

- **Explanation**:
  - `grep -v -c '^#'`: Counts non-header lines (variants) in the VCF file.
  - `bcftools view -v snps`: Filters for SNPs.
  - `bcftools view -v indels`: Filters for indels.
  - `bcftools query -f '%POS\n'`: Extracts genomic positions of variants to `pos.txt`.
- **Output**:
  - Total variants: 32,665.
  - SNPs: 31,018.
  - Indels: 1,647.
  - Variant positions: `pos.txt` (e.g., positions 78, 432, 478, etc.).

## Directory Structure

After running the pipeline, the directory contains:

```bash
P7741_R1.fastq.gz      # Input forward reads
P7741_R2.fastq.gz      # Input reverse reads
ref/                   # Reference genome directory
  sequence.fasta       # Reference genome
  sequence.fasta.*     # BWA index files
trimmed_R1.fastq.gz    # Trimmed forward reads
trimmed_R2.fastq.gz    # Trimmed reverse reads
trimmed_S.fast.qz      # Trimmed singleton reads
output.sorted.bam      # Sorted BAM file
mapping.stats.txt      # Alignment statistics
raw.bcf                # Pileup BCF file
variants.raw.vcf       # Raw variant calls
pos.txt                # Variant positions
```

## Notes

- The pipeline assumes a haploid genome (`--ploidy 1`) for *Mycobacterium ulcerans*.
- The initial `bwa mem` command had a typo in the file extension (`.qz` instead of `.gz`), which was corrected.
- Intermediate files (`output.bam`, `output.sam`) were removed to save space, and the pipeline was optimized by piping `bwa mem` directly to `samtools sort`.
- The reference genome is from *Mycobacterium ulcerans* Agy99 ([CP000325.1](https://www.ncbi.nlm.nih.gov/nuccore/CP000325.1)), as described in the study: [Whole-Genome Sequence of Mycobacterium ulcerans CSURP7741](https://pmc.ncbi.nlm.nih.gov/articles/PMC6639603/).
- Ensure read permissions for input files (e.g., `chmod +r *.fastq.gz`) to avoid errors.

## References

- FASTQ data: [ENA: ERR3335404](https://trace.ncbi.nlm.nih.gov/Traces/index.html?view=run_browser&acc=ERR3335404&display=data-access)
- Reference genome: [NCBI: CP000325.1](https://www.ncbi.nlm.nih.gov/nuccore/CP000325.1)
- Study: [Whole-Genome Sequence of Mycobacterium ulcerans CSURP7741](https://pmc.ncbi.nlm.nih.gov/articles/PMC6639603/)
