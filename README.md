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
pathway,or predisposing to specific diseases) that show a
significantly higher number of mutations in a group of samples
as compared to another group of samples.

3. How to use it
==================
From R shell:
source(".MEGA.R")
load("./example_dataset/MEGA.example.imput.Rdata")
load("./example_dataset/KEGG.186.gene.sets.Rdata")
r = MEGA(A,B,gene.sets.kegg)
head(r)
