This is a tiny tutotial on how to run a heritability analysis on a microbiome dataset with SOLAR.
Most of it was summarised by David Emmert in a useful tutorial: 
I have contributed to that, plus I am providing a few utility functions to quickly set up a microbiome heritability project in a command-line environment.

# Download and compile SOLAR

The software was previously available at several locations. The latest one was http://solar-eclipse-genetics.org/. Follow their instuctions on how to compile it. Once it is compiled, you should be able to start solar by typing `solar` on the command line. When I type `solar` in my command line, this is the information I see:

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
If you need to run one single trait, the interactive interface works just fine. The following is if you are testing a number of traits at the same time, which requires some degree of automation.

Now we can start setting up a microbiome heritability analysis

# IMPORTANT NOTICE on data structure

-  SOLAR internally uses Linear Mixed Models with a variance partitioning approach (https://www.frontiersin.org/articles/10.3389/fninf.2019.00016/full). It is highly sensitive to non-normality, especially bimodality. Many people use Centered-Log Ratios transformations, but that tends to generate Bimodal distributions. My choice is to go big or go home, with Inverse Rank Normal transformation.
-  Covariates should be all encoded as numerical/integers, including categorical data. Point estimates won't be highly affected, but 	

# Set up a run

## Transform microbiome data

The script is rather case-specific, especially because the categorical variables should be handled carefully.
In general, it should prepare a few elements

## Load and transform microbiome data

```
devtools::install_github("g-antonello/gautils")
library(gautils)

physeq <- readRDS("path/to/phyloseq/physeq.Rds")
physeq.transf <- physeq %>%
	# you can filter out rare taxa (prevalence below 20%, arbitrarily), because the distributions will look bad anyway and the model's estimates won't be reliable
	core(detection = 1, prevalence = 0.2) %>%
	phy_transform(physeq, "IRN")

work.directory <- "absolute/path/to/project"
dir.create(work.directory, recursive = T, showWarnings = FALSE)
```

## Handle covariates

In particular, categorical variables should be handled carefully. NOTE: NAs are best coded as empty fields. I think you can use NA too, but be it's at your own risk
```
physeq.transf <- physeq.transf %>%
    mutate_sample_data(

	age.rounded = round(age, 0),

        income_groups = as.numeric(income), # assuming that for example you have an ordered factor like c("< 20k", "20-50 k", "50-100k", "> 100k")

	smoking = case_when(smoking %in% c("Never", "Former") ~ 0,
		            smoking == "Current", 1,
        TRUE ~ ""
	)
    ) %>%
    
    select_sample_data(id, age.rounded, income_groups, smoking)

```
This makes you keep only the covariates you want to include in the dataset, nothing more

## Write phenotypes file in the location

```
phy_OtuMetaTable(physeq.transf) %>%
	write.csv(file.path(work.directory, "phenotypes.csv"), quote = F)
```

## Copy the pedigree file in the location

```
file.copy("where/is/the/pedigree.ped", file.path(working.directory, "pedigree.csv"))
```
The pedigree format is important! It should be structured as follows:

```
id	fa	mo	sex	
id1	id2	id3	1
id4	id5	id6	2
...			1
...			2
```

`id` is the individual's unique identifier, fa is the father, mo is the mother, sex must be coded as `males = 1`. `females = 2`

Additionally, one extra column, named `hhid` or `HHID` can be added to include the household partitioning to the moded. They can be household numbers, but I prefer to put a "h" before each house code, to be sure SOLAR encodes them as non-numeric (but it's superfluous I think)
NB: this is tab-delimited, but the documentation states that comma-delimited is the suggested format.

## write .tcl instruction files

This step is time-consuming, so make sure you automate it:

```
traits <- taxa_names(physeq.transf)
covariates <- sample_variables(physeq.transf) %>% .[!grepl("id", ., fixed = T)]

for (trait in traits){
  # create directory  
  writing.dir <- file.path(working.dir, trait)
  dir.create(writing.dir, showWarnings = F)
  
  # write tcl file
  cat(paste("proc", trait, "{} {"), 
      paste("\tload", "pedigree", file.path(working.dir, "pedigree.csv")),
      paste("\tload", "phenotypes", file.path(working.dir,"phenotypes.csv")),
      "\tmodel new",
      paste("\ttrait", trait),
      paste("\tcovariates", paste(covariates, collapse = " "), collapse = " "),
      file = file.path(writing.dir, "heritability.tcl")),
	
 # if you want the household model, uncomment the next line
      #"\thouse", 
      "\tpolygenic",
      "}",
      sep = "\n",
      file = file.path(writing.dir, "heritability.tcl")
      )
  
  cat(
    "#!/bin/bash",
    paste("cd", file.path(working.dir, trait)),
    paste("solar", trait), 
    sep = "\n",
      file = file.path(writing.dir, "solar.sh")
    )

  cat("bash", file.path(writing.dir, "solar.sh"), "\n", sep = " ", append = T, file = file.path(working.dir, "run_solar_all_array.txt")
      )
}

```
An example of a .tcl command is the following. Keep in mind that the order of the commands is important for SOLAR

```
user$ cat trait.x/heritability.tcl

proc x__ASV1 {} {
        load pedigree /home/gantonello/CHRISMB/PAPER_DRAFTS/saliva_geography_genetics/results/heritability/SOLAR/run4//chrisDataPedigree.ped
        load phenotypes /home/gantonello/CHRISMB/PAPER_DRAFTS/saliva_geography_genetics/results/heritability/SOLAR/run4//phenotypes.csv
        model new
        trait x__ASV1
        covariates x0_sex x0_ager smoking.binary x0oh01
        house
        polygenic
}
```

# Run all jobs in parallel from the command line (example with slurm)

```
cd <working.dir>
cat run_solar_all_array.txt | sarrayscript -c1 --mem-per-cpu=4G -pslow
```

each run takes 1-5 h depending on the goodness of the convergence, allocate jobs in the queues accordingly, if you have a queuing system with maximum time 
Also, by experience, for > 1k individuals it's best to allocate 4G per job, to load everything in memory. If you use less, your job may be killed. 1 thread per job is mandatory, because SOLAR doesn't support multithreading

# Retrieve results

Once results are obtained, SOLAR will write some objects in the `working.directory/traitX/traitX/polygenic.out` file.
To retrieve all at once, I designed a function that you can use after installing my R package, `gautils`

```
results.df <- get_SOLAR_results(working.directory, prefix = "")
# this next function is also useful to quickly look at the overall results. it generates a kable object
SOLAR_results_kbl(results.df, highligh_significant = 0.05)
```
the prefix parameter is useful if you put some prefix to taxa names, to tell them apart from the covariates. In this example I didn't do it, so I changed the default "x__" to "".
