Welcome to MEGA-V Repository
-------------------------------

1. Notes and How to cite MEGA-V
=================================

MEGA-V is a user friendly version of the method developed and used in Cereda et al. Patients with genetically heterogeneous synchronous colorectal cancer carry rare damaging germline mutations in immune-related genes Nature Comm. 2016; doi:10.1038/ncomms12072.

MEGA-V is implemented as a R shiny application that allows its execution also from a web-based environment. Original stand-alone version is still supported and can be found in this repository (see section 5).


2. GNU General Public License
==============================

MEGA-V is a free software that can be redistributed and/or modified under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should receive a copy of the GNU General Public License along with the program. If not, please refer to <http://www.gnu.org/licenses/>.


3. Repository Description
==========================

MEGA-V (Mutation Enrichment Gene set Analysis of Variants) identifies gene sets that show a significantly higher number of variants in a cohort of interest (cohort A) as compared to control (cohort B).

The gene sets are predefined by the user and they can be genes involved in the same biological processes from curated databases, Gene Ontology groups, or genes associated to the same disease. Two pre-processed gene sets of biological processes and diseases are provided with the software.

Similarly, also the type of gene variants to be compared are predefined by the user and 
can be damaging mutations as well as other types of genetic alterations.


4. How to run MEGA-V as a R shiny app in your browser
===============================================

Install the **Shiny** package (and the required dependencies) in R, and use the function `runGithub()`. See the example below,
```
install.packages("shiny")
install.packages("shinyjs")
install.packages("shinythemes")
install.packages("DT")
install.packages("parallel")

library(shiny)
library(shinyjs)
library(shinythemes)
library(DT)

shiny::runGitHub('ciccalab/MEGA')
```


5. How to run MEGA-V as a stand-alone app
===================================================

input:
A and B: Data frame object containing the mutations counts. 
Columns = samples; rows = mutations. The first column must always contain the name of the mutated gene.

Example:
```
1. Load the MEGA-V functions in the Global Environment
> source("./MEGA.R")

2. Load the predifined list of gene sets. As example, here we use KEGG gene sets
> load("./example_dataset/KEGG.186.gene.sets.Rdata")

3. Load the two sets of individuals A and B.
> A <- read.delim("./example_dataset/A.tsv.gz",stringsAsFactors = F)
> B <- read.delim("./example_dataset/B.tsv.gz",stringsAsFactors = F)

4. Run MEGA-V and identify which gene sets are significantly mutated in the
group of samples A as compared to B
> r = MEGA(A,B,gene.sets.kegg)

5. A summary with the input parameters used will be showed before MEGA-V start.

+----------------------------------------+
 Input parameters:
 FDR threshold: 0.1
 Number of Gene Sets: 186
 Bootstrapping: YES
 Number of iterations: 1000
+----------------------------------------+

Step 1: Enrichment Gene Set Enrichment Analysis
|===============================================| 100%

Step 2: Bootstrapping for 4 significant gene sets
|===============================================| 100%

Results:
Significant Gene sets before FDR: 26
Significant Gene sets after FDR: 4

6. Show the 4 significant pathways
> head(r,4)
```

