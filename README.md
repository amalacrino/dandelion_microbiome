# Natural rubber reduces herbivory and alters the microbiome below ground 

**Laura Böttner, Antonino Malacrinò, Christian Schulze Gronover, Nicole van Deenen, Boje Müller, Shuqing Xu, Jonathan Gershenzon, Dirk Prüfer, Meret Huber**

## Abstract
Laticifers are hypothesized to mediate both plant-herbivore and plant-microbe interactions. Yet, evidence for the dual function of these secretory structures is scarce. Here, we tested whether natural rubber, a common and economically important latex polymer, alters the performance and the root microbiome of the Russian dandelion (*Taraxacum koksaghyz*) under attack of a soil-dwelling herbivore, the larva of the May cockchafer (*Melolontha melolontha*). Rubber-depleted transgenic plants lost more shoot and root biomass upon herbivory compared to rubber-bearing near-isogenic lines. *Melolontha melolontha* larvae preferred to feed on artificial diet supplemented with latex of rubber-depleted rather than of rubber-bearing plants. Likewise, adding purified cis-1,4-polyisoprene in ecologically relevant concentrations to diet deterred larval feeding and reduced larval weight gain. Furthermore, metagenomics and metabarcoding revealed that abolishing natural rubber biosynthesis alters the structure but not the diversity of the rhizosphere and root microbiota in a herbivore- and wounding-dependent manner. Roots from rubber-deficient plants, however, did not exhibit a higher pathogen load compared to the roots of rubber-bearing plants. Taken together, our data demonstrate that natural rubber biosynthesis reduces herbivory and alters the plant microbiota in a herbivore-dependent manner, which highlights the role of plant specialized metabolites and secretory structures in shaping multitrophic interactions.

# Disclaimer

This repository contains the main components used to process the raw data and to analyze it. Raw data is available at NCBI SRA under the BioProject numbers `PRJNA779274` (16S amplicon metagenomics), `PRJNA779290` (ITS amplicon metagenomics) and `PRJNA779369` (shotgun metagenomics).

Our pipeline included:
* TrimGalore [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5127899.svg)](https://doi.org/10.5281/zenodo.5127899)
* Kraken2 (Wood et al. [2019](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1891-0))
* MegaHit (Li et al. [2016](10.1016/j.ymeth.2016.02.020))
* Prokka [![DOI:10.1093/bioinformatics/btu153](https://zenodo.org/badge/DOI/10.1093/bioinformatics/btu153.svg)](https://doi.org/10.1093/bioinformatics/btu153)
* Bowtie2 (Langmead et al. [2019](https://academic.oup.com/bioinformatics/article/35/3/421/5055585))
* VSEARCH (Rognes et al. [2016](https://peerj.com/articles/2584/))
* R  (R Core Team [2022](https://www.R-project.org/))

# Code

## 1. Amplicon metagenomics

[Process data](/amplicon/1_vsearch.md)

## 2. Shotgun metagenomics

[Pre-process data](/shotgun/1_preprocess.md)

[Gene content - data processing](/shotgun/3_functional_analysis.md)

## 3. Data analysis

[Data-analysis](/amplicon/2_R_code.md)
