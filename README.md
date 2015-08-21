Welcome to MEGA Repository
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

MEGA (Mutation Enrichment Gene set Analysis ) was developed to 
identify predefined gene sets (e.g. genes involved in the same
pathway, or predisposing to specific diseases) that show a
significantly higher number of mutations in a group of samples A
as compared to another group of samples B.

3. How to use it
==================
From R shell

inunput:
A and B: Boolean matrices of mutations. Coloums are samples, while rows are
mutations. The first coloumn must always contain the name of the gene in which
the mutation fall.

Example:

1. Load the MEGA functions in the Global Enviroment <br />
source("./MEGA.R") <br />

2. Load the predifined list of gene sets. As example, here we use KEGG gene sets <br />
load("./example_dataset/MEGA.example.imput.Rdata") <br />

3. Load the two sets of individuals A and B.
load("./example_dataset/KEGG.186.gene.sets.Rdata") <br />

4. Run MEGA and identify which gene sets are significanly mutated in the group of samples A as compared to B <br />
r = MEGA(A,B,gene.sets.kegg) <br />

5. A summary with the imput parameters used will be showed befor MEGA start.

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

5. Show the 4 significant pathways <br />
head(r,4) <br />

+-----------------------------------------------------------------------------------------------------------+
gene.set      									p.value      fdr 						success_percentage
KEGG\_CYTOKINE\_CYTOKINE\_RECEPTOR\_INTERACTION		4.903350e-06 0.0009120231               95.2
KEGG\_BIOSYNTHESIS\_OF\_UNSATURATED\_FATTY\_ACIDS	1.802473e-04 0.0111753301               72.8
KEGG\_TOLL\_LIKE_RECEPTOR\_SIGNALING\_PATHWAY		1.323134e-04 0.0111753301               79.2
KEGG\_CYTOSOLIC\_DNA\_SENSING\_PATHWAY				1.945476e-03 0.0904646149               59.6
+-----------------------------------------------------------------------------------------------------------+