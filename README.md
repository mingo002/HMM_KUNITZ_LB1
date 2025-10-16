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
- kunitz_cluster_info -->  is a cluster report containing for each sequence:

Its length (in amino acids).

Its identifier (taken from FASTA headers).

The % identity to the cluster representative (marked with '*')
- pdb_kunitz_clustered.clstr --> contains the clustered sequences. 

The retrieved sequences were cleaned using the script seq_filter.py according to these constraints:
-sequence length: between 40 and 100 residues (taken as safe range since Kuntiz domain = 50/60aa)
-unknown characters warning


##Multiple Structural Alignment (MSA)
After collecting a cleaned set of clustered sequences, MSA was calculated using PDBeFold tool.
the input file needed for this MSA should contain just the PDB ids, to do so, a new file containing just the ids of the clusterd squences was generated with this command:

grep "^>" filtered_kunitz.fasta | sed 's/^>//' > filtered_ids.txt

Then to perform the alignment the option *multiple* and *List of PDB codes* were picked in the PDBeFold website, filtered_ids was uploaded.

OUTPUTS
 superimposed_kunitz --> a FASTA file containing all the superimposed sequences
 MSA_results --> a text file with the results of the alignment (RMSD, Q-score)


 
