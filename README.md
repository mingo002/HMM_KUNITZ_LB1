# PROFILE HMM BUILDING FOR BPTI/Kunitz Domain annotation
## Preliminar requirements
Download the following dependecies useful to manage and process data
- cd-hit v.4.8.1: a fast program for clustering and comparing large sets of protein or nucleotide sequences
  conda install -c bioconda cd-hit
- hmmer v.3.3.2: search a protein sequence against a protein database (phmmer)
  conda install -c bioconda hmmer
-blast v.2.12.0+: Basic Local Alignment Search Tool, a suite of programs developed by NCBI to compare biological sequences (blastp)
  conda install -c bioconda blast
-biopython v.1.85: python tool for bioinformatics 
  conda install -c conda-forge biopython

## Step 1 - Data gathering
the first trainig set used for this model was retrieved from RCSB PDB advanced serach using the following criteria:
  -Data Collection Resolution: ≤ 3.5 Å
  -Pfam Annotation: PF00014
  -Polymer Entity Sequence Length: between 45 and 80 residues

the dataset was stored in a custom report including the following information about the protein structures:
  -Structure Data: PDB ID, Data Collection Resolution
  -Polymer Entity Data: Sequence, Auth Asym ID, Annotation Identifier,Entity ID
  
the custom report was converted in fasta format using the csv2fasta.sh script

Using Cd-hit - a tool for ***clustering*** and filtering sequences based on sequence identity - with the standard threshold (sequence identity = 90%) allowed the identification of the set of proteins (N = 25) on which to perform the multiple structural alignment while avoiding redundancy. The clustering was performed running the following command:

cd-hit -i pdb_kunitz.fasta -o pdb_kunitz_cluster.txt
this command will generate two files: 


  
  

  


