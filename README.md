
### Benchmarking DeepVariant and Conventional Variant Callers Using an Assembly-Derived Truth Set for the *Chelonia mydas* Genome

---

##  Overview

##  This repository contains all scripts and workflow files used to reproduce the analyses presented in the manuscript:
“Benchmarking DeepVariant and Conventional Variant Callers Using an Assembly-Derived Truth Set for the Chelonia mydas Genome.”

The study evaluates the performance of five short-variant callers on *Chelonia mydas* whole-genome sequencing data using a polished, assembly-derived SNP truth set. All pipelines for assembly generation, polishing, alignment, variant calling, filtering, and benchmarking are included.

---

## Repository Structure

```
Cmydas-Assembly-Derived-TruthSet/
│
├── script/                      # Scripts and pipelines used in the paper
│   ├── commands_alignment.txt
│   ├── commands_variant_calling.txt
│   ├── commands_trimming_all_tools.txt
│   ├── commands_consensus_truthset.txt
│   ├── make_snps.sh
│   ├── filter_all_tools.sh
│   ├── normal_snp_jaccard.sh            # (Not used in paper)
│   ├── normal_snp_overlap.sh            # (Not used in paper)
│   ├── polca.txt
│   └── other scripts
│
├── methods.md                 # Detailed methodology used in the manuscript
├── data_summary.md            # Final BUSCO, QUAST, ANI, and alignment results
├── LICENSE
└── README.md
```


## 🧬 Pipeline Summary (Methods Used in the Paper)

### 1. Read Preprocessing

* fastp filtering, trimming, and quality correction
* FastQC + MultiQC summary

### 2. Genome Assembly

* ABySS v2.3.6
* Best assembly selected using BUSCO (sauropsida_odb10), QUAST, and ANI

### 3. Assembly Polishing

* POLCA (MaSuRCA v4.1.4) using original Illumina reads
* Polished contigs used to derive the truth set

### 4. Reference-Based Alignment

* BWA-MEM v0.7.17 (primary aligner)
* Bowtie2 v2.5.2 (comparison)
* minimap2 v2.28 (comparison)
* BWA-MEM selected for variant calling based on mapping metrics

### 5. Variant Calling

Five short-variant callers evaluated:

* DeepVariant v1.9.0
* GATK HaplotypeCaller v4.6.2.0
* FreeBayes v1.3.10
* VarScan v2.4.6
* BCFtools v1.22

### 6. Unified Variant Filtering

Criteria used:

* FILTER = PASS or “.”
* Alternate allele depth (AD) > 3
* Variant allele frequency (VAF) > 0.02
* SNPs extracted using bcftools

### 7. Assembly-Derived Truth Set Construction

* Polished ABySS assembly aligned to reference with minimap2 (asm5)
* Variants called with BCFtools mpileup + call
* High-confidence SNPs used as the truth set

### 8. Benchmarking

* RTG Tools vcfeval used to compute

  * True positives (TP)
  * False positives (FP)
  * False negatives (FN)
  * Precision
  * Recall
  * F1-score
* Evaluations performed in:

  * Unfiltered mode (threshold = none)
  * RTG-optimised QUAL threshold mode

> 📌 No Venn diagrams, Jaccard similarity, or overlap analyses were included in this manuscript.

---

## 📥 Data Availability

Scripts and pipelines for all analyses are publicly available at:

🔗 [https://github.com/FahadAslam988/Cmydas-Assembly-Derived-TruthSet/](https://github.com/FahadAslam988/Cmydas-Assembly-Derived-TruthSet/)

Raw FASTQ files, ABySS assemblies, POLCA-corrected assemblies, and BAM alignment files are available upon reasonable request.

---

## 👤 Author

Fahad Aslam / Wasiq Aslam
GitHub: [https://github.com/FahadAslam988](https://github.com/FahadAslam988)

---

## 📄 License

Licensed under the MIT License.


Just tell me!
