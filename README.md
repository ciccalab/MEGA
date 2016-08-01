Welcome to MEGA-RVs Repository
-------------------------------

0. How to cite us
==============================

MEGA-RVs is an evolution of the original method MEGA pubblished at the following
address http://www.nature.com/ncomms/2016/160705/ncomms12072/full/ncomms12072.html 
by Cereda et al. (Nature Comm. 2016)


1. GNU General Public License
==============================

This repo contains free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This software is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.


2. Repository Description
==========================

MEGA-RVs (Mutation Enrichment Gene set Analysis of Rare Variants) was developed to 
identify predefined gene sets (e.g. genes involved in the same
pathway, or predisposing to specific diseases) that show a
significantly higher number of mutations in a group of samples A
as compared to another group of samples B.


3. How to run it in your browser as R shiny app
===============================================

The easiest way to run MEGA-RVs is to install **Shiny** package (and the required dependencies) in R, and use the function `runGithub()`. See the example below,
```
install.packages("shiny")
install.packages("shinyjs")
install.packages("shinythemes")
install.packages("DT")

library(shiny)
library(shinyjs)
library(plotly)
library(shinythemes)
library(DT)

shiny::runGitHub('ciccalab/MEGA')
```


4. How to run it in from R shell as stand-alone app
===================================================

input:
A and B: Data frame object contaning the mutations counts. Coloums are samples,
while rows are mutations. The first coloumn must always contain the name of the
gene in which the mutation fall.

Example:
```
1. Load the MEGA-RVs functions in the Global Enviroment
> source("./MEGA.R")

2. Load the predifined list of gene sets. As example, here we use KEGG gene sets
> load("./example_dataset/MEGA.example.imput.Rdata")

3. Load the two sets of individuals A and B.
> load("./example_dataset/KEGG.186.gene.sets.Rdata")

4. Run MEGA-RVs and identify which gene sets are significanly mutated in the
group of samples A as compared to B
> r = MEGA(A,B,gene.sets.kegg)

5. A summary with the imput parameters used will be showed befor MEGA-RVs start.

+----------------------------------------+
 Input parameters:
 FDR threshold: 0.1
 Number of Gene Sets: 186
 Bootstrapping: YES
 Number of iterations: 1000
+----------------------------------------+

Step 1: Enrichement Gene Set Enrichement Anlysis
|===============================================| 100%

Step 2: Bootstrapping for 4 significant gene sets
|===============================================| 100%

Results:
Significant Gene sets before FDR: 26
Significant Gene sets after FDR: 4

6. Show the 4 significant pathways
> head(r,4)
```