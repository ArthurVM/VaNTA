AUTHORS NOTE : REPO IN WORKS. UNFINISHED.

VaNTA
-----
v0.1.1 May 2020 by Arthur V. Morris

Automated discovery and analysis of Tandem Repeats.

This pipeline was developed to automate the process of VNTR discovery. A description of the automated workflow is as follows:

- **Step 1.** Identify Tandem Repeats (TRs) within coding regions of a reference genome.
- **Step 2.** Generate a reference library of these TRs using Tandem Repeats Finder (Benson, 1999).
- **Step 3.** Detect the presence of these TRs within a dataset of query annotated genomes.
- **Step 4.** Align the TRs detected within the query annotated genomes with the corresponding reference TR in a pairwise  manner using EMBOSS-water (part of the EMBOSS: The European Molecular Biology Open Software Suite (2000) Rice, Longden and Bleasby).
- **Step 5.** Identify and score Variable Number Tandem Repeats (VNTRs) using the alignment data.

Dependancies
---------------------
Python Dependencies:
- Python3
- Biopython
- argsparse

Other Dependencies:
- This pipeline uses Tandem Repeats Finder and EMBOSS-water.
- Tandem Repeats Finder (Benson, 1999) can be downloaded from https://github.com/Benson-Genomics-Lab/TRF
- EMBOSS-water (part of the EMBOSS: The European Molecular Biology Open Software Suite (2000) Rice, Longden and Bleasby) can be downloaded from http://emboss.sourceforge.net/
    
Notes on input:
---------------
- All annotation files must be in .gff format.
- All annotation files MUST have a consistent naming convention for their features. There is no sequence homology inference in this pipeline, it instead relies on homologous features being named identically across the whole dataset. It may therefore be prudent to transfer annotations to the query dataset using an annotation liftover tool. Alternatively, you may wish to use a homolog clustering tool (e.g. MCScan) to identify homologous genes, and then assign them some kind of consistent naming convention.

**TR scoring:**
The V-score for each TR locus is calculated to bias smaller polymorphic VNTR's with low repeat motif length (2 to 6 bases) and high flank conservation.

**Output:**
A .tsv file of all Tandem Repeats (VNTRs or otherwise) ordered by V-score.
