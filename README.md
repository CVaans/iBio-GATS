# iBio-GATS (Inverted Bio-GATS): 
## A semi-automated workflow for building high-quality structural models of insect ORs, based on the Bio-GATS (Jabeen et al.) approach.
### Getting started 
**Install BLAST**
The NCBI provides a suite of command-line tools to run BLAST called BLAST+. This allows users to perform BLAST searches on their own server without size, volume, and database restrictions. BLAST+ can be used with a command line so it can be integrated directly into your workflow.

-  Download and install BLAST command line tools [here](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)

**BLAST Help**
-  BLAST documentation [here](https://biopython.readthedocs.io/en/latest/chapter_blast.html)
-  Running BLAST from Python [here](https://gtpb.github.io/PPB18/assets/15_Running-BLAST_sys.argv)
-  Download and install the BLAST+ package from [here](http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
-  BLAST manual can be found [here](http://www.ncbi.nlm.nih.gov/books/NBK1762/)
-  Standalone BLAST Setup for Windows PC [here](https://www.ncbi.nlm.nih.gov/books/NBK52637/)
-  Standalone BLAST Setup for Unix [here](https://www.ncbi.nlm.nih.gov/books/n/helpblast/chapter1/)
-  Blast FTP Site [here](https://www.ncbi.nlm.nih.gov/books/NBK62345/)
-  BLAST Command Line Applications User Manual [here](https://www.ncbi.nlm.nih.gov/books/NBK279688/)

**Install Modeller 10.4**

MODELLER 10.4 is a tool for modelling protein three-dimensional structures based on homology or comparative analysis. An alignment file containing the alignment of a sequence to be modelled with known related structure/structures is given as input for the Modeller for it to automatically calculate a model for the given target sequence. Modeller is available for download for Unix/Linux, Windows, and Mac operating systems.                                       

- Current version Modeller 10.4 released on Nov 2022 can be downloaded from [here](https://salilab.org/modeller/download_installation.html)

**Modeller Help** 
- Modeller documentation is found [here](https://salilab.org/modeller/documentation.html)
- Modeller manual is found [here](https://salilab.org/modeller/manual/)
- Modeller tutorial is found [here](https://salilab.org/modeller/tutorial/)
- Modeller registration is found [here](https://salilab.org/modeller/registration.html) 

**Install dependencies for iBio-GATS**
**Checkout project folders**
**pip install -r requirements.txt**


**Prerequisites:**
**Input file:**
A target sequence from _Drosophila melanogaster_ odorant receptor sequence (sequence.fasta) to be kept in iBio-GATS folder by the user

**Run bash script ibio-gats.sh**

**Outputs:**
1. Result summaries showing hydrophobicity plots and helical plots for each template-target alignment.
1. Total of 5 Model output .pdb files for the target sequence modelled using selected template stored in the modeller-files folder

