# iBio-GATS (Inverted Bio-GATS): 
## A workflow for building high-quality structural models of insect ORs, based on the Bio-GATS (Jabeen et al.) approach.
### Getting started 
__ Install BLAST __
The NCBI provides a suite of command-line tools to run BLAST called BLAST+. This allows users to perform BLAST searches on their own server without size, volume, and database restrictions. BLAST+ can be used with a command line so it can be integrated directly into your workflow.
-  Download and install BLAST command line tools here
BLAST Help
-  BLAST documentation here
-  Running BLAST from Python here
-  Download and install the BLAST+ package from here
-  BLAST manual can be found here
-  Standalone BLAST Setup for Windows PC here
-  Standalone BLAST Setup for Unix here
-  Blast FTP Site here
-  BLAST Command Line Applications User Manual here

__ Install Modeller 10.4 __

MODELLER 10.4 is a tool for modelling protein three-dimensional structures based on homology or omparative analysis. An alignment file containing the alignment of a sequence to be modelled with known related structure/structures is given as input for the Modeller for it to automatically calculate a model for the given target sequence. Modeller is available for download for Unix/Linux, Windows, and Mac operating systems.                                       

- Current version Modeller 10.4 released on Nov 2022 can be downloaded from here 
- Modeller Help 
- Modeller documentation is found here
- Modeller manual is found here
- Modeller tutorial is found here
- Modeller registration is found here 

__ Install dependencies for iBio-GATS __
__ Checkout project folders __
__ pip install -r requirements.txt __


__Prerequisites: __
Input file: A target sequence from _Drosophila melanogaster_ odorant receptor sequence (sequence.fasta) to be kept in iBio-GATS folder by the user

__ Run bash script ibio-gats.sh __

__ Outputs: __
1. Result summaries showing hydrophobicity plots and helical plots for each template-target alignment.
1. Total of 5 Model output .pdb files for the target sequence modelled using selected template stored in the modeller-files folder

