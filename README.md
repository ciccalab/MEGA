Welcome to MEGA-RVs Repository
-------------------------------

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

runGitHub("MEGA","ciccalab")
```


3. How to run it in from R shell as standalon app
=================================================

input:
A and B: Data frame object contaning the mutations counts. Coloums are samples,
while rows are mutations. The first coloumn must always contain the name of the
gene in which the mutation fall.

Example:
```
1. Load the MEGA-RVs functions in the Global Enviroment <br />
source("./MEGA.R") <br />

2. Load the predifined list of gene sets. As example, here we use KEGG gene sets <br />
load("./example_dataset/MEGA.example.imput.Rdata") <br />

3. Load the two sets of individuals A and B.
load("./example_dataset/KEGG.186.gene.sets.Rdata") <br />

4. Run MEGA-RVs and identify which gene sets are significanly mutated in the group of samples A as compared to B <br />
r = MEGA(A,B,gene.sets.kegg) <br />

5. A summary with the imput parameters used will be showed befor MEGA-RVs start.

+----------------------------------------+<br />
 Input parameters:<br />
 FDR threshold: 0.1 <br />
 Number of Gene Sets: 186 <br />
 Bootstrapping: YES<br />
 Number of iterations: 1000 <br />
+----------------------------------------+<br />

Step 1: Enrichement Gene Set Enrichement Anlysis<br />
|===============================================| 100%

Step 2: Bootstrapping for 4 significant gene sets<br />
|===============================================| 100%

Results:<br />
Significant Gene sets before FDR: 26<br />
Significant Gene sets after FDR: 4<br />

6. Show the 4 significant pathways<br />
head(r,4) <br />
```