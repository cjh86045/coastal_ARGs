# ARGs Analysis Pipeline

A **four-stage** metagenomic pipeline for identifying, quantifying, and host-associating **Antibiotic Resistance Genes (ARGs)** from paired-end sequencing data:

1. **Preprocessing** (`preprocessing.py`) — Quality control, assembly, ORF prediction, CARD annotation, and ARG-carrying contig (ACC) extraction.
2. **Abundance profiling** (`arg_abundance.py`) — Estimate cell counts (ncells) and compute ARG abundance normalized as **copies/cell**.
3. **Binning** (`binning.py`) — Reconstruct metagenome-assembled genomes (MAGs/bins) and assess their quality.
4. **Bin annotation** (`bin_annotation.py`) — Filter high-quality bins, assign taxonomy with GTDB-tk, and link bins to ARGs.

## Overview

### Stage 1 — Preprocessing

1. **fastp** — Raw reads quality control
2. **MEGAHIT** — Metagenome assembly
3. **Prodigal** — ORF prediction
4. **DIAMOND** — Alignment against the CARD database
5. **ARG-like ORFs extraction** — Extract ORFs that align to CARD
6. **ARG-carrying contigs (ACCs) extraction** — Extract contigs carrying ARGs
7. **ACC ORF prediction** — Predict ORFs on ACCs

### Stage 2 — Abundance profiling

1. **count_cells** — Align clean reads to KO30 single-copy marker genes (DIAMOND BLASTx) to estimate cell counts (ncells)
2. **run_diamond** — Align clean reads to per-sample ARG-like ORFs (DIAMOND BLASTx)
3. **calculate_abundance** — Compute ARG abundance normalized by ncells (gene-level + antibiotic-type-level)
4. **merge_results** — Merge all sample results into a single summary

### Stage 3 — Binning

1. **filter_needless_sample** — Remove samples lacking contigs longer than 1000 bp
2. **generate_bam_and_depth** — Align clean reads to contigs (bowtie2), sort BAM, compute per-contig coverage (jgi_summarize_bam_contig_depths)
3. **run_metabat** — Bin contigs with MetaBAT2
4. **run_maxbin** — Bin contigs with MaxBin2
5. **run_concoct** — Bin contigs with CONCOCT (cut contigs → coverage → cluster → extract bins)
6. **run_dastool** — Integrate the three bin sets with DAS Tool to obtain optimized bins
7. **run_checkm** — Assess bin completeness / contamination with CheckM

### Stage 4 — Bin annotation

1. **get_quality_bins** — Select high-quality bins (Completeness > 50%, Contamination < 10%) from CheckM results
2. **extract_quality_bins** — Copy high-quality bin FASTA files into per-sample folders
3. **run_gtdbtk** — Merge all high-quality bins (with sample prefix) and run GTDB-tk **once** for taxonomic classification
4. **merge_results** — Combine quality metrics + GTDB-tk taxonomy + ARG-host linkage

## Requirements

- Python 3.7+
- [fastp](https://github.com/OpenGene/fastp)
- [MEGAHIT](https://github.com/voutcn/megahit)
- [Prodigal](https://github.com/hyattpd/Prodigal)
- [DIAMOND](https://github.com/bbuchfink/diamond)
- [Bowtie2](https://github.com/BenLangmead/bowtie2)
- [samtools](https://github.com/samtools/samtools)
- [jgi_summarize_bam_contig_depths](https://bitbucket.org/berkeleylab/metabat/src/master/)
- [MetaBAT2](https://bitbucket.org/berkeleylab/metabat/)
- [MaxBin2](https://sourceforge.net/projects/maxbin2/)
- [CONCOCT](https://github.com/BinPro/CONCOCT)
- [DAS_Tool](https://github.com/cmks/DAS_Tool)
- [CheckM](https://github.com/checkm-genome/checkm1)
- [GTDB-tk](https://github.com/Ecogenomics/GTDBTk)
- [Biopython](https://biopython.org/)
- [pandas](https://pandas.pydata.org/)
- [tqdm](https://github.com/tqdm/tqdm)

## Configuration

Before running the pipeline, configure the following in `configure.py`:

- **CARD_DB** — Path to the CARD database (DIAMOND format)
- **CARD_STRUCTURE** — Path to the CARD metadata table (TSV)
- **KO30_DB** — Path to the KO30 single-copy marker gene database (DIAMOND format)
- **KO30_STRUCTURE** — Path to the KO30 annotation file (TSV: `sseqid` \t `ko30`)
- **THREAD_CONFIG** — Thread settings for each tool (e.g., `fastp`, `megahit`, `prodigal`, `diamond`)
- **THREAD_CONFIG** — Pool size for parallel sample processing (e.g., `fastp_pool`, `megahit_pool`, `diamond_pool`)
- **THREAD_CONFIG (binning)** — Binning tool threads & pools: `bowtie2`, `samtools`, `metabat`, `maxbin`, `concoct`, `dastool`, `checkm`, plus `*_pool` sizes
- **THREAD_CONFIG (bin annotation)** — `gtdbtk` thread count

## Input Data Format

Place raw paired-end reads in the `data/` folder under the project directory. Files must follow the naming convention:

```
{sample_name}_1.fastq.gz   # Forward reads
{sample_name}_2.fastq.gz   # Reverse reads
```

## Usage

### Stage 1 — Preprocessing

```python
from preprocessing import Preprocessing

preprocessing = Preprocessing(input_folder="/path/to/your/project")
preprocessing.run()
```

Or run directly:

```bash
python preprocessing.py
```

### Stage 2 — Abundance profiling

```python
from arg_abundance import ARGsAbundanceProfile

abundance = ARGsAbundanceProfile(input_folder="/path/to/your/project")
abundance.run()
```

Or run directly:

```bash
python arg_abundance.py
```

> **Note:** Update the `input_folder` path in the `__main__` section before running. Stage 2 depends on the outputs of Stage 1.

### Stage 3 — Binning

```python
from binning import Binning

binning = Binning(input_folder="/path/to/your/project")
binning.run()
```

Or run directly:

```bash
python binning.py
```

> Depends on Stage 1 (assembly, ACC) and Stage 2 (ARG-like ORFs) outputs.

### Stage 4 — Bin annotation

```python
from bin_annotation import BinAnnotation

annotator = BinAnnotation(input_folder="/path/to/your/project")
annotator.run()
```

Or run directly:

```bash
python bin_annotation.py
```

> Depends on Stage 3 outputs (`dastool_out/`, `checkm_out/`) and Stage 1 ACC outputs.

---

## Output Directory Structure

After running the full pipeline, the following directory structure is generated:

```
project/
├── data/                                    # Input data directory
│   ├── sample1_1.fastq.gz                  # Raw forward reads
│   ├── sample1_2.fastq.gz                  # Raw reverse reads
│   ├── sample2_1.fastq.gz
│   └── sample2_2.fastq.gz
│
├── preprocessing/                           # Stage 1 outputs
│   ├── preprocessing.log                    # Pipeline log
│   │
│   ├── fastp_output/                        # Step 1: fastp QC
│   │   ├── clean_sample1_1.fastq.gz        # Clean forward reads
│   │   ├── clean_sample1_2.fastq.gz        # Clean reverse reads
│   │   ├── sample1_fastp.json              # QC report (JSON)
│   │   ├── sample1_fastp.html              # QC report (HTML)
│   │   ├── clean_sample2_1.fastq.gz
│   │   ├── clean_sample2_2.fastq.gz
│   │   ├── sample2_fastp.json
│   │   └── sample2_fastp.html
│   │
│   ├── megahit_output/                      # Step 2: MEGAHIT assembly
│   │   ├── sample1_final.contigs.fa        # Final assembled contigs
│   │   └── sample2_final.contigs.fa
│   │
│   ├── prodigal_output/                     # Step 3: Prodigal ORF prediction
│   │   ├── sample1_prodigal.gkb            # Gene coordinates (GBK format)
│   │   ├── sample1_prodigal.faa            # Protein sequences
│   │   ├── sample2_prodigal.gkb
│   │   └── sample2_prodigal.faa
│   │
│   ├── ARGs_blast_result_output/            # Step 4: DIAMOND vs CARD
│   │   ├── sample1_card_results.tsv        # Alignment results
│   │   └── sample2_card_results.tsv
│   │
│   ├── arg_like_orfs_output/                # Step 5: ARG-like ORFs
│   │   ├── sample1_arg_like_orfs.faa       # ORFs that aligned to CARD
│   │   └── sample2_arg_like_orfs.faa
│   │
│   ├── acc_output/                          # Step 6: ARG-carrying contigs
│   │   ├── sample1_acc_contigs.fa          # Contigs carrying ARGs
│   │   └── sample2_acc_contigs.fa
│   │
│   └── acc_orfs/                            # Step 7: ORFs on ACCs
│       ├── sample1_acc_orfs.faa            # Protein sequences on ACCs
│       ├── sample1_acc_orfs.gbk            # Gene coordinates on ACCs
│       ├── sample2_acc_orfs.faa
│       └── sample2_acc_orfs.gbk
│
├── args_abundance/                          # Stage 2 outputs
│   ├── args_abundance_profile.log           # Abundance pipeline log
│   │
│   ├── count_cells_output/                  # Phase 1: KO30 BLASTx results
│   │   ├── sample1_count_cells_results.tsv # Clean reads vs KO30 markers
│   │   └── sample2_count_cells_results.tsv
│   │
│   ├── arg_like_orfs_diamond_output/        # Phase 2: ARG-like ORFs BLASTx results
│   │   ├── sample1_diamond_results.tsv     # Clean reads vs ARG-like ORFs
│   │   └── sample2_diamond_results.tsv
│   │
│   ├── abundance_output/                    # Phase 3: Gene-level abundance
│   │   ├── sample1_abundance.tsv           # Per-ARG abundance (copies/cell)
│   │   └── sample2_abundance.tsv
│   │
│   ├── type_abundance_output/               # Phase 3: Antibiotic-type abundance
│   │   ├── sample1_type_abundance.tsv      # Per-type abundance (copies/cell)
│   │   └── sample2_type_abundance.tsv
│   │
│   └── merge_results/                       # Phase 4: Merged summary
│       └── merged_abundance_results.tsv    # Total abundance per sample
│
├── binning/                                 # Stage 3 outputs
│   ├── binning.log                          # Binning pipeline log
│   │
│   ├── get_depth_index_output/              # bowtie2 indexes (per sample)
│   ├── get_depth_bam_output/                # Sorted BAM + BAI (per sample)
│   ├── get_depth_output/                    # Per-contig coverage (*_depth.txt)
│   │
│   ├── metabat_output/                      # MetaBAT2 bins (per sample)
│   │   └── {sample}_metabat/
│   │       └── bin.1.fa, bin.2.fa, ...
│   │
│   ├── maxbin_output/                       # MaxBin2 bins (per sample)
│   │   └── {sample}_maxbin/
│   │       └── bin.001.fasta, ...
│   │
│   ├── concoct_output/                      # CONCOCT bins (per sample)
│   │   └── {sample}_concoct/
│   │       ├── contigs_10K.fa               # 10kbp-cut contigs
│   │       ├── coverage_table.tsv           # Coverage table
│   │       ├── clustering_merged.csv        # Merged clustering
│   │       └── fasta_bins/                  # Final bins (*.fa)
│   │
│   ├── dastool_out/                         # DAS Tool integrated bins
│   │   └── {sample}_dastool/
│   │       ├── output_DASTool_bins/         # Optimized bins
│   │       ├── output_DASTool_contig2bin.tsv# contig → bin mapping
│   │       └── DASTool_summary.tsv
│   │
│   └── checkm_out/                          # CheckM quality assessment
│       └── {sample}_checkm/
│           └── storage/
│               └── bin_stats_ext.tsv        # Completeness / contamination
│
└── bin_annotation/                          # Stage 4 outputs
    ├── bin_annotation.log                   # Bin annotation log
    │
    ├── quality_bins_results/                # High-quality bin list (per sample)
    │   └── {sample}_quality_bins_results.tsv
    ├── bin_count_sample_results/            # Quality bin count (per sample)
    │   └── bin_count_sample_results.tsv
    ├── quality_bins_fasta/                  # High-quality bin FASTA (per sample)
    │   └── {sample}_quality_bins/
    │       └── *.fa
    ├── merged_quality_bins/                 # ⭐ All quality bins merged (sample-prefixed)
    │   └── {sample}__{bin}.fa
    ├── gtdb_results/                        # GTDB-tk results (single run)
    │   └── gtdbtk_all/
    │       └── gtdb_.bac120.summary.tsv     # Taxonomic classification
    ├── final_results_with_species/          # Final per-sample results
    │   └── {sample}_final_results_with_species.tsv
    └── sample_merge_results/                # ARG-host summary across samples
        └── sample_merge_results.tsv
```

## Output Files Description

### Stage 1 — Preprocessing

| Step | Folder | File Pattern | Description |
|------|--------|--------------|-------------|
| 1 | `fastp_output/` | `clean_{sample}_*.fastq.gz` | Quality-controlled reads |
| 1 | `fastp_output/` | `{sample}_fastp.json/html` | QC reports |
| 2 | `megahit_output/` | `{sample}_final.contigs.fa` | Assembled contigs |
| 3 | `prodigal_output/` | `{sample}_prodigal.gkb` | ORF coordinates |
| 3 | `prodigal_output/` | `{sample}_prodigal.faa` | ORF protein sequences |
| 4 | `ARGs_blast_result_output/` | `{sample}_card_results.tsv` | DIAMOND vs CARD results |
| 5 | `arg_like_orfs_output/` | `{sample}_arg_like_orfs.faa` | ARG-like ORF sequences |
| 6 | `acc_output/` | `{sample}_acc_contigs.fa` | ARG-carrying contigs |
| 7 | `acc_orfs/` | `{sample}_acc_orfs.faa` | ORFs on ACCs (proteins) |
| 7 | `acc_orfs/` | `{sample}_acc_orfs.gbk` | ORFs on ACCs (coordinates) |

### Stage 2 — Abundance profiling

| Phase | Folder | File Pattern | Description |
|-------|--------|--------------|-------------|
| 1 | `count_cells_output/` | `{sample}_count_cells_results.tsv` | KO30 marker alignment results |
| 2 | `arg_like_orfs_diamond_output/` | `{sample}_diamond_results.tsv` | ARG-like ORFs alignment results |
| 3 | `abundance_output/` | `{sample}_abundance.tsv` | Gene-level abundance (copies/cell) |
| 3 | `type_abundance_output/` | `{sample}_type_abundance.tsv` | Antibiotic-type-level abundance |
| 4 | `merge_results/` | `merged_abundance_results.tsv` | Per-sample total abundance summary |

### Stage 3 — Binning

| Step | Folder | File Pattern | Description |
|------|--------|--------------|-------------|
| depth | `get_depth_output/` | `{sample}_depth.txt` | Per-contig coverage |
| depth | `get_depth_bam_output/` | `{sample}.sorted.bam[.bai]` | Sorted BAM + index |
| 3 | `metabat_output/` | `{sample}_metabat/bin.*.fa` | MetaBAT2 bins |
| 4 | `maxbin_output/` | `{sample}_maxbin/bin.*.fasta` | MaxBin2 bins |
| 5 | `concoct_output/` | `{sample}_concoct/fasta_bins/*.fa` | CONCOCT bins |
| 6 | `dastool_out/` | `{sample}_dastool/output_DASTool_bins/` | DAS Tool optimized bins |
| 6 | `dastool_out/` | `{sample}_dastool/output_DASTool_contig2bin.tsv` | contig → bin mapping |
| 7 | `checkm_out/` | `{sample}_checkm/storage/bin_stats_ext.tsv` | CheckM completeness / contamination |

### Stage 4 — Bin annotation

| Step | Folder | File Pattern | Description |
|------|--------|--------------|-------------|
| 1 | `quality_bins_results/` | `{sample}_quality_bins_results.tsv` | High-quality bins + quality metrics |
| 2 | `quality_bins_fasta/` | `{sample}_quality_bins/*.fa` | High-quality bin FASTA |
| 3 | `merged_quality_bins/` | `{sample}__{bin}.fa` | Merged bins (sample-prefixed) |
| 3 | `gtdb_results/gtdbtk_all/` | `gtdb_.bac120.summary.tsv` | GTDB-tk taxonomic classification |
| 4 | `final_results_with_species/` | `{sample}_final_results_with_species.tsv` | Final per-sample result (with `contains_ARG`) |
| 4 | `sample_merge_results/` | `sample_merge_results.tsv` | ARG-host summary across all samples |

### Abundance table columns (`{sample}_abundance.tsv`)

| Column | Description |
|--------|-------------|
| sseqid | ARG gene ID (from CARD) |
| Type | Antibiotic type / class |
| Resistance_Mechanism_Type | Resistance mechanism category |
| abundance(copies/cell) | ARG abundance normalized by cell count |
| sample_name | Sample name |

### Merged summary columns (`merged_abundance_results.tsv`)

| Column | Description |
|--------|-------------|
| sample_name | Sample name |
| total_abundance(copies/cell) | Sum of all ARG abundances in the sample |

### Final bin table columns (`{sample}_final_results_with_species.tsv`)

| Column | Description |
|--------|-------------|
| bin_name | Bin name |
| Completeness | Bin completeness (%) |
| Contamination | Bin contamination (%) |
| Genome size | Genome size (bp) |
| classification | GTDB-tk taxonomic classification |
| contains_ARG | Whether the bin carries any ARG (True/False) |

## DIAMOND Output Format

The DIAMOND BLASTx alignment output uses `--outfmt 6` with 8 columns:

```
qseqid sseqid pident length qlen slen evalue bitscore
```

| Column | Description |
|--------|-------------|
| qseqid | Query sequence ID (read ID) |
| sseqid | Subject sequence ID (marker / ARG gene ID) |
| pident | Percentage of identical matches |
| length | Alignment length |
| qlen | Query (read) length |
| slen | Subject (reference) length |
| evalue | Expectation value |
| bitscore | Bit score |

> **Note:** Adjust the `names=[...]` list in `pd.read_table` to match your actual `--outfmt` column set; the default BLAST tabular format has 12 columns, but this pipeline uses an 8-column subset.

## Abundance Calculation

### Cell count estimation (ncells)

Clean reads are aligned to 30 universal single-copy marker genes (KO30). For each marker gene, the total coverage is summed; dividing by 30 gives the estimated cell count:

```
ncells = Σ_over_30_markers ( Σ_reads(length / slen) ) / 30
```

- `length` — alignment length of the read on the marker
- `slen` — full length of the marker gene
- `length / slen` — coverage contributed by one read
- `/ 30` — average across the 30 single-copy marker genes ≈ number of cells

### ARG abundance (copies/cell)

Clean reads are aligned to the per-sample ARG-like ORFs. After filtering, the abundance of each ARG is normalized by ncells:

```
abundance(copies/cell) = (qlen / slen) / ncells
```

- `qlen` — read length
- `slen` — ARG reference sequence length
- `qlen / slen` — copies of the ARG contributed by one read
- `/ ncells` — normalized to per-cell abundance

#### Filtering thresholds

| Criterion | Threshold |
|-----------|-----------|
| pident | ≥ 80% |
| evalue | ≤ 1e-7 |
| length | ≥ 8 aa |
| qcov (length / qlen) | ≥ 0.25 |

#### Aggregation

- **Gene level** — Sum `abundance(copies/cell)` grouped by ARG gene (`sseqid`)
- **Type level** — Sum `abundance(copies/cell)` grouped by antibiotic type (`Type`)
- **Sample level** — Sum all gene-level abundances → `total_abundance(copies/cell)`

> **Note:** Risk-index weighting has been removed; abundance is reported as raw copies/cell.

## Binning & Annotation Details

### High-quality bin selection

Bins are filtered by CheckM metrics:

```
Completeness > 50%  AND  Contamination < 10%
```

### Bin merging & naming (Stage 4, key point)

All per-sample high-quality bins are copied into a **single directory** (`merged_quality_bins/`) and run through GTDB-tk **once**. To avoid name collisions and preserve sample-of-origin, every bin is renamed with a sample prefix:

```
Naming rule:  {sample_name}__{original_bin_name}.fa
Example:      sampleA__bin.1.fa,  sampleB__bin.1.fa
```

GTDB-tk's `user_genome` column uses the same naming, so `merge_results` splits it back:

```
user_genome = "sampleA__bin.1"
  →  sample_name = "sampleA"
  →  bin_name    = "bin.1"
```

> **Separator note:** The double underscore `__` is the delimiter. If your original bin names themselves contain `__`, change the separator (e.g., to `||`) in both `_prepare_merged_bins` and `merge_results` to avoid ambiguous splitting.

### ARG-host linkage

For each sample, ACC contig IDs (from `acc_output/`) are intersected with the `contig2bin` mapping (`output_DASTool_contig2bin.tsv`) to identify which bins carry ARGs:

```
arg_bins = contig2bin[ contig2bin.contig_id ∈ ACC_contig_ids ].bin_name.unique()
```

The resulting `contains_ARG` column marks bins that host at least one ARG-carrying contig.

## Pipeline Features

- **Parallel processing** — Sample-level parallelism via a configurable process pool (Stage 1 & 3); serial processing with progress bar where computation is cheap (Stage 2 abundance & extraction)
- **Single-run GTDB-tk** — All high-quality bins are merged with sample-prefixed names and classified in one GTDB-tk invocation, avoiding redundant taxonomic placement
- **Resume support** — MEGAHIT supports `--continue`; every step skips samples whose non-empty output already exists
- **Progress tracking** — Real-time `tqdm` progress bars with success counts
- **Robust error handling** — Failed samples are logged but do not stop the pipeline
- **Direct DIAMOND database** — `.faa` files are supplied directly to `--db` (no separate `makedb` step required)

## Known Considerations

- **DIAMOND `--query` with paired files** — Both `count_cells` and `run_diamond` pass forward and reverse reads as `--query fwd rev`. This works with versions that accept multiple query files; if your DIAMOND build rejects it, concatenate the paired reads or run two separate alignments and merge the results.
- **Empty abundance outputs** — Samples with no valid ARG alignments after filtering produce empty `*_abundance.tsv` files (marked as processed). `merge_results` will record these with `total_abundance = 0`.
- **GTDB-tk summary filename** — This pipeline uses the `--prefix gtdb_` output (`gtdb_.bac120.summary.tsv`). Some GTDB-tk ≥2.x versions produce `gtdbtk.bac120.summary.tsv` instead; adjust the filename in `merge_results` accordingly.
- **Bin name separator** — The `__` delimiter used when merging bins for GTDB-tk must not appear inside original bin names; change it on both the merge and split sides if needed.
- **process_samples_in_parallel return value** — This helper returns a list of `(sample_name, success)` tuples; call sites must receive the list first, then iterate, rather than using `for result in process_samples_in_parallel(...)`.

## License

This project is for research purposes. Please ensure you comply with the licenses of the respective tools used in this pipeline.
