# Process of the data in NingLab Skin

## To-do

<details>
<summary>Pre-process the data using samtools and picard</summary>

  - [ ] Create BAM files from SAM files
    <details>
    <summary>Details</summary>

    - [x] Code  
    - [ ] Run  

    </details>

  - [ ] Add read groups to BAM files
    <details>
    <summary>Details</summary>

    - [x] Code  
    - [ ] Run  

    </details>

  - [ ] Sort BAM files
    <details>
    <summary>Details</summary>

    - [x] Code  
    - [ ] Run  

    </details>

  - [ ] Index sorted BAM files
    <details>
    <summary>Details</summary>

    - [x] Code  
    - [ ] Run  

    </details>

  - [ ] Mark and remove duplicates in BAM files
    <details>
    <summary>Details</summary>

    - [ ] Code  
    - [ ] Run  

    </details>

</details>

<details>
<summary>Call variants using GATK</summary>

  - [ ] BQSR using `GATK BaseRecalibrator`
    <details>
    <summary>Details</summary>

    - [ ] Code  
    - [ ] Run  

    </details>

  - [ ] Apply BQSR using `GATK ApplyBQSR`
    <details>
    <summary>Details</summary>

    - [ ] Code  
    - [ ] Run  

    </details>

  - [ ] Call variants using `GATK HaplotypeCaller`
    <details>
    <summary>Details</summary>

    - [ ] Code  
    - [ ] Run  

    </details>

  - [ ] Filter variants using `GATK VariantFiltration`
    <details>
    <summary>Details</summary>

    - [ ] Code  
    - [ ] Run  

    </details>

<summary>Performing GWAS analysis using plink</summary>

  - [ ] Convert VCF to PLINK format
    <details>
    <summary>Details</summary>

    - [ ] Code  
    - [ ] Run  

    </details>

  - [ ] Perform quality control on PLINK data
    <details>
    <summary>Details</summary>

    - [ ] Code  
    - [ ] Run  

    </details>

  - [ ] Perform GWAS analysis using PLINK
    <details>
    <summary>Details</summary>

    - [ ] Code  
    - [ ] Run  

    </details>

  - [ ] PCA analysis using PLINK
    <details>
    <summary>Details</summary>

    - [ ] Code  
    - [ ] Run  

    </details>

<summary>Using ANNOVAR to annotate variants</summary>

  - [ ] Annotate variants using `ANNOVAR`
    <details>
    <summary>Details</summary>

    - [ ] Code  
    - [ ] Run  

    </details>

  - [ ] Filter annotated variants
    <details>
    <summary>Details</summary>

    - [ ] Code  
    - [ ] Run  

    </details>

</details>
