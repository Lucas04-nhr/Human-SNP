# Process of the data in NingLab Skin

## To-do

- [x] Align the data with reference genome GRCh38
    - [x] Code
    - [x] Run
- [ ] Pre-process the data using `samtools` and `picard`
  - [ ] Create BAM files from SAM files
    - [x] Code
    - [ ] Run
  - [ ] Add read groups to BAM files
    - [x] Code
    - [ ] Run
  - [ ] Sort BAM files
    - [x] Code
    - [ ] Run
  - [ ] Index sorted BAM files
    - [x] Code
    - [ ] Run
  - [ ] Mark and remove duplicates in BAM files
    - [ ] Code
    - [ ] Run
- [ ] Call variants using `GATK`
  - [ ] BQSR using `GATK BaseRecalibrator`
    - [ ] Code
    - [ ] Run
  - [ ] Apply BQSR using `GATK ApplyBQSR`
    - [ ] Code
    - [ ] Run
  - [ ] Call variants using `GATK HaplotypeCaller`
    - [ ] Code
    - [ ] Run
  - [ ] Filter variants using `GATK VariantFiltration`
    - [ ] Code
    - [ ] Run
