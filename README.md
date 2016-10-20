Welcome to MEGA-V Repository
-------------------------------

1. General Notes and How to cite MEGA-V
=======================================

MEGA-V (Mutation Enrichment Gene set Analysis of Variants) is a development of the method used in Cereda et al. (2016) Nature Comm. 7; doi:10.1038/ncomms12072.

MEGA-V is a R application with a Shiny web interface that allows its execution from a web-based environment (section 4). The stand-alone version is also supported (section 5). Examples on how to use MEGA-V are provided (section 6).


2. GNU General Public License
=============================

MEGA-V is a free software that can be redistributed and/or modified under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should receive a copy of the GNU General Public License along with the program. If not, please refer to <http://www.gnu.org/licenses/>.


3. MEGA-V Description
=====================

MEGA-V identifies gene sets with a significantly higher number of variants in a cohort of interest (cohort A) as compared to (1) a control cohort (cohort B) or (2) a random distribution generated using Monte Carlo.

Gene sets are predefined by the user and they can be group of genes involved in the same biological processes, Gene Ontology groups, or genes associated to the same disease. 

The type of genetic variants to be compared are also predefined by the user and can be damaging alterations, rare variants, or other types of genetic modifications.

(1) If a control cohort B is used, the user can choose to compare the two distributions of mutation counts between the two cohorts using either the Wilcoxon-rank sum test or the Kolmogorov–Smirnov test according to the data distribution. If multiple gene sets are tested, the resulting p-values are corrected for multiple testing. When the sample size between the two cohorts differs substantially, MEGA-V can run a random sampling with replacement where the larger cohort is randomly down-sampled to the size of the smaller cohort.

(2) If a random distribution of variants in cohort A is used, Monte Carlo permutations are applied to estimate the expected mutation counts within each gene set across all individuals of cohort A. An empirical p-value is then measured as the fraction of observed mutation counts that is greater or equal than the estimated mutation counts. 





4. How to run MEGA-V as a R shiny app in your browser
=====================================================

Install the **Shiny** package (and the required dependencies) in R, and use the function `runGithub()`. See the example below,
```
install.packages("shiny")
install.packages("shinyjs")
install.packages("shinythemes")
install.packages("DT")
install.packages("parallel") # for Monte Carlo

library(shiny)
library(shinyjs)
library(shinythemes)
library(DT)

shiny::runGitHub('ciccalab/MEGA')
```

5. How to run MEGA-V as a stand-alone app
=========================================

Clone the repository on your local machine and open R from the repository folder. Load the following files in the R Global Environment:

``` 
source(“functions/MEGA.R”)
source(“functions/MEGA_MC_libs.R”)
```

6. Exemplar usage of MEGA-V
===========================

Two lists of 186 biological processes (ncomm.cereda.186.KEGG.gmt) and 346 disease-associated gene sets (ncomm.cereda.346.gwas.gmt) are provided in the folder ‘gene_sets’.

Two example datasets of gene variants for cohort A (ncomm.cereda.syCRCs.tsv) and cohort B (ncomm.cereda.1000.genomes.tsv) are 
provided in the folder “example_dataset”. These derive from Cereda et al. Nature Comm 2016 and can be used as exemplar files to run MEGA-V. 

Example 1: Use MEGA-V to identify altered gene sets using a control cohort 
```
1. Load the MEGA-V functions in the Global Environment
> source("./functions/MEGA.R")

2. Load the predefined list of gene sets. As example, here we use KEGG gene sets
> gene.sets.kegg = read.gmt.file("gene_sets/ncomm.cereda.186.KEGG.gmt")

3. Load the two sets of patients with synchronous CRCs and individuals of 1000 Genomes Project. A and B: Data frame object containing the mutations counts. Columns = samples; rows = mutations. The first column must always contain the name of the mutated gene.

> A <- read.delim("./example_dataset/ncomm.cereda.syCRCs.tsv.gz",stringsAsFactors = F)
> B <- read.delim("./example_dataset/ncomm.cereda.1000.genomes.tsv.gz",stringsAsFactors = F)

4. Run MEGA-V and identify which gene sets are significantly mutated in the
group of samples A as compared to B
> r = MEGA(A,B,gene.sets.kegg, bootstrapping=TRUE, nsim=1000)
```

Example 2: Use MEGA-V to identify altered gene sets without a control cohort
```
1. Load the MEGA-V functions in the Global Environment
> source("./functions/MEGA.R")

2. Load the predefined list of gene sets. As example, here we use KEGG gene sets
> gene.sets.kegg = read.gmt.file("gene_sets/ncomm.cereda.186.KEGG.gmt")

3. Load the set of patients with synchronous CRCs 
> A <- read.delim("./example_dataset/ncomm.cereda.syCRCs.tsv.gz",stringsAsFactors = F)

4. Run MEGA-V and identify which gene sets are significantly mutated in the group of samples A 
> r = MEGA(A,gene.sets.kegg, montecarlo=TRUE, nsim=1000)
```

