This is a tiny tutotial on how to run a heritability analysis on a microbiome dataset with SOLAR.
Most of it was summarised by David Emmert in a useful tutorial: 
I have contributed to that, plus I am providing a few utility functions to quickly set up a microbiome heritability project in a command-line environment.

# Download and compile SOLAR

The software was previously available at several locations. The latest one was http://solar-eclipse-genetics.org/ but I think it got discontinued in the summer of 2023. This is the information I have in my SOLAR build:

```
SOLAR Eclipse version 8.5.1 (beta), last updated on July 31, 2020
Developed at Maryland Psychiatric Research Center,
University of Maryland School of Medicine, Baltimore.
Visit our documentation and tutorial website www.solar-eclipse-genetics.org
Our download page https://www.nitrc.org/projects/se_linux
Our github page https://github.com/brian09/solar-eclipse
For questions email: pkochunov@gmail.com
Enter help for help, exit to exit, doc to browse documentation.
The software development is supported by NIH grant RO1EB015611
from The National Institute for Biomedical Imaging and Bioengineering.
Enter cite to see how to cite this software
```
Follow their instuctions on how to compile it. Once it is compiled, you should be able to start solar by typing `solar` on the command line

Now we can start setting up a microbiome heritability analysis

# IMPORTANT NOTICE on data structure

	* SOLAR internally uses Linear Mixed Models with a variance partitioning approach (https://www.frontiersin.org/articles/10.3389/fninf.2019.00016/full). It is highly sensitive to non-normality, especially bimodality. Many people use Centered-Log Ratios transformations, but that tends to generate Bimodal distributions. 
	* Covariates should be all encoded as numerical/integers, including categorical data. Point estimates won't be highly affected, but 	

# Set up a run

## Transform microbiome data

The script is rather case-specific, especially because the categorical variables should be handled carefully.
The script does a few things:

    1) Reads in the phyloseq object
    2) Filters and transforms the microbiome data

## Handle covariates

## Write phenotypes file in the location

## Copy the pedigree file in the location

## write .tcl instruction files

# Retrieve results

