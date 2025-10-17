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

the dataset was stored in a custom report rcsb_pdb_custom_report_kunitz.csv including the following information about the protein structures:
  -Structure Data: PDB ID, Data Collection Resolution
  -Polymer Entity Data: Sequence, Auth Asym ID, Annotation Identifier,Entity ID
  

rcsb_pdb_custom_report_kunitz.csv was converted in fasta format using the csv2fasta.sh script 
the output fasta file is kunits_sequencs.fasta containing PDBID_chain and the sequences (160)

Using Cd-hit - a tool for ***clustering*** and filtering sequences based on sequence identity - with the standard threshold (sequence identity = 90%) allowed the identification of the set of proteins (N = 25) on which to perform the multiple structural alignment while avoiding redundancy. The clustering was performed running the following command:

cd-hit -i kunitz_sequences.fasta -o pdb_kunitz_cluster.txt
this command will generate two files: 
- kunitz_cluster_info.clstr -->  is a cluster report containing for each sequence:

Its length (in amino acids).

Its identifier (taken from FASTA headers).

The % identity to the cluster representative (marked with '*')
- kunitz_cluster.txt --> contains the clustered sequences. 

The retrieved sequences were cleaned using the script seq_filter.py according to these constraints:
-sequence length: between 40 and 100 residues (taken as safe range since Kuntiz domain = 50/60aa)
-unknown characters warning

with this script two files have been generated:
- filtered_kunitz.fasta -->  contains the filterd sequences more reliable to kunitz doamin annotation
- filterd_out --> contains the sequences which did not respected the constrints cited above (2ODY-->  too long 127aa, 5JBT--> too short 38aa)

Another entry was manually removed (4BQD:A) because its A chain lacks one or both β-strands typical of Kunitz fold; likely a partial or degraded domain rather than a full, functional one.

##Multiple Structural Alignment (MSA)
After collecting a cleaned set of clustered sequences, MSA was calculated using PDBeFold tool.
the input file needed for this MSA should contain just the PDB ids, to do so, a new file containing just the ids of the clusterd squences was generated with this command:

grep '^>' filtered_kunitz.fasta | sed 's/^>//' | sed 's/_/:/'  > filtered_ids.txt

the output file will contain PDBID:chain

Then to perform the alignment the option *multiple* and *List of PDB codes* were picked in the PDBeFold website, filtered_ids was uploaded.

OUTPUTS
 MSA_22.fasta --> a FASTA file containing all the 22 superimposed sequences


The MSA output file was cleaned: alla the residues were converted into uppercase characters, the 'PDB' prefix was removed, and the output file clean_MSA_22.fasta was re-printed without adding a trailing newline

awk '{if (substr($1,1,1)==">") {print "\n"toupper($1)} else {printf "%s" toupper($1)}}' MSA_22.fasta | sed s/PDB://| tail -n +2 > clean_MSA_22.fasta
 
## HMM Building
  profile HMM construction from multiple sequence alignments

  hmmbuild pdb_kunitz.hmm clean_MSA_22.fasta

## Validation set creation
  from uniprot website human and non human kunitz proteins were downloaded in fasta fromat using the following queries:
  -Human: (taxonomy_id:9606) AND (reviewed:true) AND (xref:pfam-PF00014) --> 18 entries

  -Non-human: NOT (taxonomy_id:9606) AND (reviewed:true) AND (xref:pfam-PF00014) --> 380 entries

  the files were then merged in all_kuntiz.fasta
  cat human_kunitz.fasta non-human.fasta > all_kunitz.fasta

  Create a blast database with the kunitz proteins from SwissProt
  makeblastdb -in all_kunitz.fasta -input_type fasta -dbtype prot -out all_kunitz.fasta

  Perform BLAST search on the 22 representative sequences against the wholw kunitz dataset
  blastp -query clean_MSA_22.fasta -db all_kunitz.fasta -out bpti_kunitz_clean.blast -outfmt 7