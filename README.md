## Materials for ESA '21 ipmr Workshop

This repository hosts the slides and tutorials for the `ipmr` workshop from ESA 2021. Tutorials are stored in `ipmr_examples/` subfolder. To run

1. clone this repository into a folder on your computer and open `ipmr_esa.Rproj`. 

2. Open `ipmr_examples/ipmr_tutorial.rmd`. 

3. In RStudio, you can simply click the `Run Document` button above the script pane. To run from the command line in an R session, run `rmarkdown::run("ipmr_examples/ipmr_tutorial.rmd")`.

## Prerequisites

In order to get the most out of this tutorial, you will need to be familiar with the regression modeling and when to use different link functions (e.g. Logit, Log). You will also need to understand some basic IPM theory, such as what projection kernels are, how projection kernels are comprised of sub-kernels, and how regression models are combined to create sub-kernels. There are a variety of good publications to familiarize yourself with IPM theory and practice. [Merow et al. 2014](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12146) is a good start. Though not required for this tutorial, more advanced treatments can be found in [Ellner & Rees 2006](https://www.journals.uchicago.edu/doi/pdfplus/10.1086/499438?casa_token=hVkM-U1RHs0AAAAA:YduKcTwNVRUviS5j1soKVQ62bSTNLFN8Cx-9mQTdju4yov83XNHFJNRaXptLMfDbhQUKrWFf9HI), [Rees & Ellner 2009](https://esajournals.onlinelibrary.wiley.com/doi/10.1890/08-1474.1), and Ellner, Childs, & Rees 2016. 

You will also need to install the `ipmr`, `DiagrammeR`, and `learnr` R packages:

```
install.packages(c("ipmr", "DiagrammeR","learnr"))

```


