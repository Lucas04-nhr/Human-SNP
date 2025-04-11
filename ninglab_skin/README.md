# Process of the data in NingLab Skin

## To-do

- [x] Align the data with reference genome GRCh38
    - [x] Code
    - [x] Run
- [ ] Pre-process the data using `samtools` and `picard`
  - [x] Create BAM files from SAM files
    - [x] Code
    - [x] Run
  - [x] Add read groups to BAM files
    - [x] Code
    - [x] Run
  - [x] Sort BAM files
    - [x] Code
    - [x] Run
  - [x] Index sorted BAM files
    - [x] Code
    - [x] Run
  - [x] Mark and remove duplicates in BAM files
    - [x] Code
    - [ ] Run
- [ ] Call variants using `GATK`
  - [ ] BQSR using `GATK BaseRecalibrator`
    - [x] Code
    - [ ] Run
  - [ ] Apply BQSR using `GATK ApplyBQSR`
    - [x] Code
    - [ ] Run
  - [ ] Call variants using `GATK HaplotypeCaller`
    - [x] Code
    - [ ] Run
  - [ ] Filter variants using `GATK VariantFiltration`
    - [ ] Code
    - [ ] Run
- [ ] Annotate variants using `ANNOVAR`
    - [ ] Code
    - [ ] Run
