# brainxcan

S-BrainXcan takes GWAS as input and return the association between GWAS phenotype and a list of brain image-derived phenotypes.

BrainXcan manuscript [link](https://www.medrxiv.org/content/10.1101/2021.06.01.21258159v1).

Software documentation (Quick start) [link](https://liangyy.github.io/brainxcan-docs/docs/index.html).

BrainXcan database [link](http://doi.org/10.5281/zenodo.4895174).

Analysis scripts [link](https://github.com/liangyy/ukb_idp_genetic_arch).

**[https://github.com/liangyy/brainxcan](https://github.com/liangyy/brainxcan) is deprecated.** Please go to [https://github.com/hakyimlab/brainxcan](https://github.com/hakyimlab/brainxcan) for the latest BrainXcan software.

# Installation notes

## Software dependencies

The software is built upon both Python and R scripts along with some standalone executables.
Here we provide a conda environment containing all the Python dependencies and `snakemake`.

```
conda env create -f environment.yml
``` 

Also, install [plink 1.9](https://www.cog-genomics.org/plink/) which will be used for LD clumping in MR analysis.

By default, the pipeline call `python`, `Rscript`, and `plink` as is.
And you can provide the path to the desired executables in the configuration file. For instance,

```
# in config.[name].yaml
rscript_exe: 'path-to/Rscript' 
python_exe: 'path-to/python'
plink_exe: 'path-to/plink'
``` 

## Standalone R

R dependencies are: `ggplot2`, `dplyr`, `optparse`, `logging`, `rmarkdown`, `pdftools`, `patchwork`, `oro.nifti`, `data.table`, `pander`, `arrow`, `TwoSampleMR`, `qvalue`.

Below, we provide an example for installing R dependencies as a conda environment. 
Any standalone R installation with these dependent packages being installed should work just fine.

```
$ conda create -n r_36 -y
$ conda activate r_36
(r_36) $ conda install -c r r
(r_36) $ conda install -c conda-forge r-arrow
(r_36) $ conda install -c conda-forge r-pdftools
(r_36) $ conda install -c conda-forge r-gmp
(r_36) $ conda install -c conda-forge r-rio
(r_36) $ conda install -c conda-forge r-pander
(r_36) $ conda install -c conda-forge r-sf
(r_36) $ conda install -c conda-forge r-stars
(r_36) $ conda install -c conda-forge r-plotly
(r_36) $ conda install -c conda-forge r-ggnewscale
(r_36) $ R
# inside R
> install.packages(c('ggplot2', 'dplyr', 'logging', 'optparse', 'rmarkdown', 'patchwork', 'oro.nifti', 'data.table', 'remotes', 'raster', 'rgeos'))
> remotes::install_github("MRCIEU/TwoSampleMR")
> devtools::install_github("jdstorey/qvalue")
```






