<p align="center">
  <img src="https://github.com/crahal/architectureofGWAS/blob/master/figures/png/helix.png" width="550"/>
</p>

# The Architecture of GWAS

[![Generic badge](https://img.shields.io/badge/Python-3.6-<red>.svg)](https://shields.io/)  [![Generic badge](https://img.shields.io/badge/License-MIT-blue.svg)](https://shields.io/)  [![Generic badge](https://img.shields.io/badge/Maintained-Yes-green.svg)](https://shields.io/)

### Introduction

This is a repository to accompany a scientometric review of all GWAS by M.C. Mills and C. Rahal. The article which corresponds to this repository is now published Open Access in [Communications Biology!](https://www.nature.com/commsbio/) A 'shortlink' to the paper can be found here:
<p align="center">
https://rdcu.be/bfLdY
</p>
Please see the notebook (and associated paper) for a description of what the code does, but in summary, it undertakes a systematic data-driven review of all GWAS as of 27/10/2018. This repo can be cloned and ran on-the-fly as required to update the results with only minimal adjustments. Dictionaries for regular expression based exercises will be maintained as new Catalog releases introduce new terms needing normalization. Please raise all issues related to the code here in this repo.

### Prerequisites

As a pre-requisite to running the iPython notebook, you will need a working Python 3 installation with all of the necessary dependencies detailed in [requirements.txt](https://github.com/crahal/architectureofGWAS/blob/master/requirements.txt). We strongly recommend virtual environments and [Anaconda](https://www.anaconda.com/distribution/). Specific requirements: you will at least need ```conda install basemap -c conda-forge```, ```conda install basemap-data-hires``` (```config --add channels conda-forge```), ```pip install gender_guesser```, ```pip install biopython```, ```pip install requests_ftp```,```pip install tqdm```, ```pip install unidecode``` and ```pip install wordcloud```. Comment out cell the choropleth map subsection for any issues with basemap compatibility.

### Running the Code

This code is operating system independent (through the ``os`` module) and should work on Windows, Linux and OS X all the same. For those unfamiliar with jupyter notebooks, just ```cd``` to the the Code subdirectory in the architectureofGWAS clone, and then run the ```jupyter notebook``` command in the terminal to launch the notebook. A recent [presentation](https://github.com/crahal/GWASReview/tree/master/presentation) accompanies the repository.

### Repository Structure

Repo structure made using the ```tree``` [utility](https://en.wikipedia.org/wiki/Tree_%28Unix%29).

├── src  
│   ├── GWASReview_Notebook.ipynb  
│   └── Support  
│   │   ├── Additional.py  
│   │   ├── Analysis.py  
│   │   ├── Ancestry.py  
│   │   ├── Figures.py  
│   │   ├── LoadData.py  
│   │   ├── PubMed.py  
│   │   └── Robustness.py  
├── data  
│   ├── Catalogue  
│   │   ├── Raw  
│   │   └── Synthetic  
│   ├── PUBMED  
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
│       ├── doublehelix_mask.png  
│       └── native_classifier_dictionary.csv  
├── Figures  
│   ├── pdf  
│   │   └── helix_wordcloud_1250_5000_black.pdf  
│   ├── png  
│   │   └── helix.png  
│   └── svg  
│       ├── Figure_1.svg  
│       ├── Figure_2.svg  
│       ├── Figure_3.svg  
│       ├── Figure_3_Blues.svg  
│       ├── Figure_4.svg  
│       ├── Figure_4_NoCommas.svg  
│       ├── Sup_Figure_1_Commas.svg  
│       └── Sup_Figure_1_NoCommas.svg  
├── piprequirements.txt  
├── readme.md  
├── GWAS_Review.html  
├── Tables  
│   ├── Authors.csv  
│   ├── Broad_Ancestral_Full.csv  
│   ├── Broad_Ancestral_NoNR.csv  
│   ├── Broad_Ancestral_Time_NoNR_PC.csv  
│   ├── ContinentOfRecruitment.csv  
│   └── CountryOfRecruitment.csv

### Versioning
This is the version which powers the conditionally accepted article. Note: some of the manually parsed observations on utilised datasets may require further cleaning and correcting in future versions of the Catalog. Please raise an issue in this repo to alert us of any necesarry changes.

### License
This work is free. You can redistribute it and/or modify it under the terms of the MIT license. It comes without any warranty, to the extent permitted by applicable law.

### Acknowledgments
Research assistance for the manual data curation was provided by Pilar Wiegand, Xuejie Ding and Domantė Grendaitė. Please see the accompanying paper for a full list of acknowledgements. A code review was undertaken by [IanKnowles](https://github.com/ianknowles)

### To Do (Long Run)

1. Custom function to create funder-specific plots for each funder
2. Complete manual parsing task for 100% of all studies.
3. Move dictionaries from support functions to text delimited files
