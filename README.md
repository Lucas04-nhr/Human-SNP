# Human SNP

This is a repository for the Human SNP project.

See notes at WPS cloud document: [Human SNP](https://kdocs.cn/l/ctgvmxKPKfYD).

## Mirrored Repositories

- [GitHub](https://github.com/Lucas04-nhr/Human-SNP)
- [Gitee](https://gitee.com/lucas04/Human-SNP)

---

## Dataset Location

```
data
├── 00_bwa_index                  # Index files and reference genome, known sites for the pipeline
│   ├── chm13v2                   # Reference genome and known sites for CHM13v2, deprecated
│   │   └── known_sites
│   └── GRCh38                    # Reference genome and known sites for GRCh38
│       ├── known-sites           # Known sites for GRCh38, using when performing SNP calling
│       │   ├── 1000g             #
│       │   ├── dbsnp138          # Homo_sapiens_assembly38_modified.dbsnp138.vcf(.idx)
│       │   ├── hapmap            #
│       │   └── omni              #
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
├── skin_metagenome
├── skin_microbiome_2022
└── test
```
