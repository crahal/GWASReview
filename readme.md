### The Anatomy of GWAS

[![Generic badge](https://img.shields.io/badge/Python-3.6-<red>.svg)](https://shields.io/)  [![Generic badge](https://img.shields.io/badge/License-MIT-blue.svg)](https://shields.io/)  [![Generic badge](https://img.shields.io/badge/Maintained-Yes-green.svg)](https://shields.io/)


### Introduction

This is a repository to accompany the paper 'The Anatomy of Genome Wide Association Studies' by M.C. Mills and C.Rahal. For a copy of the working paper, please e-mail *melinda dot mills at sociology dot ox dot ac dot uk*. A link to an open-access version will appear here in due course. Please see the notebook (and associated paper) for a description of what code does, but in summary, it undertakes a systematic data-driven review of all GWAS as of 01/03/2018. This repo can be cloned and ran on-the-fly as required to update the results with only minimal updates and adjustments. Dictionaries for regex based exercise will be maintained as new catalog releases introduce new terms needing normalization.


### Prerequisites

As a pre-requisite to running the iPython notebook, you will need a working Python 3 installation with all of the necessary dependancies detailed in requirements.txt (pip install the  [requirements.txt](http:\\github.com\crahal\anatomyofGWAS\requirements.txt)). We strongly recommend [Anaconda](http:\\continuum.io\anaconda). Comment out cell [14] for any issues with basemap, where condaforge --install basemap might be required for Windows. For those unfamiliar with juputer notebooks,  just ```cd``` to the cloned directory, and then run the ```jupyter notebook``` command in the terminal to launch the notebook.

### Running the Code

This code is operating system independent and should work on Windows, Linux and Mac all the same. For those unfamiliar with juputer notebooks,  just ```cd``` to the the Code subdirectory in the anatomyofGWAS clone, and then run the ```jupyter notebook``` command in the terminal to launch the notebook.

#### Repository Structure

Repo structure made using the Linux ```tree``` [utility](https://en.wikipedia.org/wiki/Tree_%28Unix%29).



├── Code  
│   ├── Notebook.ipynb  
│   └── Support  
│       ├── Additional.py  
│       ├── Ancestry.py  
│       ├── LoadData.py  
│       ├── PubMed.py  
│       └── Robustness.py  
├── Data  
│   ├── Catalogue  
│   │   ├── Raw  
│   │   │   ├── Cat_Anc.tsv  
│   │   │   ├── Cat_Full.tsv  
│   │   │   ├── Cat_Map.tsv  
│   │   │   └── Cat_Stud.tsv  
│   │   └── Synthetic  
│   │       ├── ancestry_CoR.csv  
│   │       ├── GWASCatalogue_CleanedCountry.tsv  
│   │       ├── Mapped_Traits_EFO.csv  
│   │       ├── Mapped_Traits_EFO_NoCommas.csv  
│   │       ├── new_initial_sample.csv  
│   │       └── new_replication_sample.csv  
│   ├── PUBMED  
│   │   ├── Pubmed_AbstractCount.csv  
│   │   ├── Pubmed_AbstractInfo.csv  
│   │   ├── Pubmed_AuthorInfo.csv  
│   │   ├── Pubmed_Cites.csv  
│   │   ├── Pubmed_CollectiveInfo.csv  
│   │   └── Pubmed_FunderInfo.csv  
│   ├── ShapeFiles  
│   │   ├── Country_Lookup.csv  
│   │   ├── ne_10m_admin_0_countries.cpg  
│   │   ├── ne_10m_admin_0_countries.dbf  
│   │   ├── ne_10m_admin_0_countries.prj  
│   │   ├── ne_10m_admin_0_countries.README.html  
│   │   ├── ne_10m_admin_0_countries.shp  
│   │   ├── ne_10m_admin_0_countries.shx  
│   │   └── ne_10m_admin_0_countries.VERSION.txt  
│   └── Support  
│       ├── Author_Supplmentary.csv  
│       ├── Cohorts  
│       │   ├── Dataset_Frequencies.csv  
│       │   ├── Dictionary_cohorts.csv  
│       │   ├── GWASCat1to1000 final.csv  
│       │   └── manual_cohort_for_merge.csv  
│       ├── Collectives  
│       │   ├── Pubmed_CollectiveInfo_Dictionary.csv  
│       │   └── Pubmed_CollectiveInfo_unverified.csv  
│       ├── dict_replacer_broad.xlsx  
│       └── doublehelix_mask.png  
├── Figures  
│   ├── pdf  
│   │   └── helix_wordcloud_1250_5000_black.pdf  
│   └── svg  
│       ├── Ancestral_SubPlots.svg  
│       ├── Country_Rec_robin_OrRed_N_iresolution.svg  
│       ├── Diseases_Traits_EFO_Bars.svg  
│       ├── Funder_Heatmap.svg  
│       ├── Gender_by_Subjects.svg  
│       └── GWAS_Popularity.svg  
├── piprequirements.txt  
├── readme.md  
├── Tables  
│   ├── Authors.csv  
│   ├── Broad_Ancestral_Full.csv  
│   ├── Broad_Ancestral_NoNR.csv  
│   ├── ContinentOfRecruitment.csv  
│   └── Funder_Ancestry.csv  




#### Versioning
This is a prototype version prior to full release and publication of the associated paper.

##### License
This work is free. You can redistribute it and/or modify it under the terms of the MIT license. It comes without any warranty, to the extent permitted by applicable law.

#### Acknowledgments
Research assistance for the manual data curation was provided by Pilar Wiegand, Xuejie Ding and Domantė Grendaitė.
