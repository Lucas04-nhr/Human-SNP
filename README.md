# Human SNP

This is a repository for the Human SNP project.

See notes at WPS cloud document: [Human SNP](https://kdocs.cn/l/ctgvmxKPKfYD).

## Mirrored Repositories

- [GitHub](https://github.com/Lucas04-nhr/Human-SNP)
- [Gitee](https://gitee.com/lucas04/Human-SNP)

---

## Dataset Location

### Raw Data

``` bash
data
├── 00_bwa_index                  # Index files and reference genome, known sites for the pipeline
│   ├── chm13v2                   # Reference genome and known sites for CHM13v2, deprecated
│   │   └── known_sites
│   └── GRCh38                    # Reference genome and known sites for GRCh38
│       ├── known-sites           # Known sites for GRCh38, using when performing SNP calling
│       │   ├── 1000g             # hg38_v0_1000G_phase1.snps.high_confidence.modified.hg38.vcf(.idx)
│       │   ├── dbsnp138          # hg38_v0_Homo_sapiens_assembly38.dbsnp138.modified.vcf(.idx)
│       │   ├── hapmap            # hg38_v0_hapmap_3.3.hg38.modified.vcf(.idx)
│       │   └── omni              # hg38_v0_1000G_omni2.5.hg38.modified.vcf(.idx)
│       └── ref                   # 1000G reference genome for GRCh38, using when performing genome phasing
│           └── 1000G
├── 01_align                      # Alignment results, end with *.sam
│   ├── Beijing
│   └── Guangzhou
├── 02_filter                     # Alignment results, filtered out unmapped reads, end with *.bam
│   ├── Beijing
│   └── Guangzhou
├── 03_sort                       # Alignment results, sorted by read name, end with *.bam, index files end with *.bam.bai
│   ├── Beijing
│   └── Guangzhou
├── 05_BaseRecalibrator           # Base recalibration results, end with *.recal_data.table
│   ├── Beijing
│   └── Guangzhou
├── 06_ApplyBQSR                  # BQSR results, end with *.recalibrated.bam
│   ├── Beijing
│   └── Guangzhou
├── 07_HaplotypeCaller            # HaplotypeCaller results, end with *.called.vcf.gz
│   ├── Beijing
│   └── Guangzhou
├── skin_metagenome               # Skin metagenome original data, end with *.fastq.gz
├── skin_microbiome_2022          # Compressed skin metagenome original data, end with *.fastq.gz
└── test                          # Simlink files to sample BJ001, used for testing the pipeline
```

### Pipeline Scripts

``` bash
Human-SNP
├── filter                        # Filter out unmapped reads
│   ├── filter
│   └── mapq_analysis
├── gatk                          # GATK pipeline, used for SNP calling
│   ├── combine                   # Combine VCF files to GVCF files
│   ├── single                    # Single calling
│   ├── joint                     # Joint calling
│   └── prepare                   # Prepare scripts for GATK pipeline
├── manta                         # Manta pipeline, used for SV calling, deprecated
│   ├── 00_convert_format
│   ├── 01_configurate
│   └── 01_execute
├── old                           # Old scripts, deprecated
│   ├── BQSR
│   ├── bwa
│   ├── env
│   ├── files
│   ├── other
│   └── samtools
└── prepare                       # Prepare scripts for WEScall pipeline, deprecated

