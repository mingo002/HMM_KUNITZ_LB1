# PROFILE HMM BUILDING FOR BPTI/Kunitz Domain annotation
## Preliminar requirements
Download the following dependecies useful to manage and process data
- cd-hit v.4.8.1: a fast program for clustering and comparing large sets of protein or nucleotide sequences
  conda install -c bioconda cd-hit
- hmmer v.3.3.2: search a protein sequence against a protein database (phmmer)
  conda install -c bioconda hmmer
-blast v.2.12.0+: Basic Local Alignment Search Tool, a suite of programs developed by NCBI to compare biological sequences 
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
-sequence length: between 40 and 100 residues (taken as safe range since Kunitz domain = 50/60aa)
-unknown characters warning

with this script two files have been generated:
- filtered_kunitz.fasta -->  contains the filterd sequences more reliable to kunitz doamin annotation
- filterd_out --> contains the sequences which did not respected the constrints cited above (2ODY:E -->  too long 127aa, 5JBT:Y--> too short 38aa)

Another entry was manually removed (4BQD:A) because its A chain lacks one or both β-strands typical of Kunitz fold; likely a partial or degraded domain rather than a full, functional one.

## Multiple Structural Alignment (MSA)
After collecting a cleaned set of clustered sequences, MSA was calculated using PDBeFold tool.
the input file needed for this MSA should contain just the PDB ids, to do so, a new file containing just the ids of the clusterd squences was generated with this command:

grep '^>' filtered_kunitz.fasta | sed 's/^>//' | sed 's/_/:/'  > filtered_ids.txt

the output file will contain PDBID:chain

Then to perform the alignment the option *multiple* and *List of PDB codes* were picked in the PDBeFold website, filtered_ids was uploaded.

OUTPUTS
 MSA_22.fasta --> a FASTA file containing all the 22 superimposed sequences


The MSA output file was cleaned: all the residues were converted into uppercase characters, the 'PDB' prefix was removed, and the output file clean_MSA_22.fasta was re-printed without adding a trailing newline

awk '{if (substr($1,1,1)==">") {print "\n"toupper($1)} else {printf "%s" toupper($1)}}' MSA_22.fasta | sed s/PDB://| tail -n +2 > clean_MSA_22.fasta
 
## HMM Building
  profile HMM construction from multiple structural alignments

  hmmbuild pdb_kunitz.hmm clean_MSA_22.fasta

## Validation set creation
  from uniprot website human and non human kunitz proteins were downloaded in fasta fromat using the following queries:
  (reviewed:true) AND (xref:pfam-PF00014) --> 380 entries

  the output was exported in all_kuntiz.fasta

  Create a blast database with the kunitz proteins from SwissProt
  makeblastdb -in all_kunitz.fasta -input_type fasta -dbtype prot -out all_kunitz.fasta

  Perform BLAST search on the 22 representative sequences against the whole kunitz dataset
  blastp -query clean_MSA_22.fasta -db all_kunitz.fasta -out bpti_kunitz_clean.blast -outfmt 7

  
## Method optimization and assessment.

  At this point balanced positive and negative test sets for Kunitz protein validation were prepared by removing redundant sequences, shuffling, splitting, and extracting the relevant FASTA entries.

  HMM search on positive and negative test sets, relevant E-values extraction, sequence labeling, undetected negatives with high E-values insertion and combination of results into final classification files for further analysis.
  ### 1. Remove Highly Similar Sequences
      Purpose: Identify and remove sequences with high identity (≥95%) and alignment length ≥50 from the dataset to avoid redundancy.
      How:
          Parse bpti_kunitz_clean.blast.
          Use awk to filter for identity and alignment length.
          Extract the second column (subject ID), sort uniquely, and extract the Uniprot ID.
          Save these IDs to ids_to_remove.txt.
      
      grep -v "^#" bpti_kunitz_clean.blast | awk '{if ($3>=95 && $4>=50) print $2}' | sort -u | cut -d "|" -f 2 > ids_to_remove.txt

  ### 2. Extract All Kunitz IDs
      Purpose: Get all Uniprot IDs from the Kunitz dataset.
          How:
              Search for header lines in all_kunitz.fasta.
              Extract the Uniprot ID from each header.
              Save to all_kunitz_ids.txt.
      grep ">" all_kunitz.fasta | cut -d "|" -f 2 > all_kunitz_ids.txt
  
  ### 3. Create Positive Test Sets
      Purpose: Generate two positive test sets (set1 and set2) of Kunitz proteins, excluding those with high similarity.
      How:
      Remove IDs in ids_to_remove.txt from all_kunitz_ids.txt using comm -23.
      Shuffle the remaining IDs.
      Split into two sets: first 183 IDs (pos_set_1.txt), last 182 IDs (pos_set_2.txt).
      Extract corresponding sequences from uniprot_sprot.fasta using get_seq.py, saving to pos_1.fasta and pos_2.fasta.

      Download the whole Swissprot database from Uniprot website and export it in uniprot_sprot.fasta

      comm -23 <(sort all_kunitz_ids.txt) <(sort ids_to_remove.txt) > ids_to_keep.txt
      sort -R ids_to_keep.txt > kunitz_shuffled.txt

      total=$(wc -l < kunitz_shuffled.txt)
      half=$(( (total + 1) / 2 ))   # This rounds up since kunitz ids are odd
      head -n $half kunitz_shuffled.txt > pos_set_1.txt
      tail -n $((total - half)) kunitz_shuffled.txt > pos_set_2.txt kunitz_shuffled.txt > pos_set_2.txt
      python3 get_seq.py pos_set_1.txt uniprot_sprot.fasta > pos_1.fasta
      python3 get_seq.py pos_set_2.txt uniprot_sprot.fasta > pos_2.fasta

  ### 4. Create Negative Test Sets
      Purpose: Generate two negative test sets (set1 and set2) of non-Kunitz Swiss-Prot proteins.
      How:
      Extract all Uniprot IDs from uniprot_sprot.fasta.
      Remove all Kunitz IDs to get negative candidates.
      Shuffle the negative IDs.
      Split into two sets: first 286,286 IDs (neg_set_1.txt), last 286,286 IDs (neg_set_2.txt).
      Extract corresponding sequences from SwissProt using get_seq.py, saving to neg_1.fasta and neg_2.fasta.
  
      grep ">" uniprot_sprot.fasta | cut -d "|" -f 2 > all_sp_ids.txt
      comm -23 <(sort all_sp_ids.txt) <(sort all_kunitz_ids.txt) > negative_ids.txt
      sort -R negative_ids.txt > negative_shuffled.txt
      head -n 286286 negative_shuffled.txt > neg_set_1.txt
      tail -n 286286 negative_shuffled.txt > neg_set_2.txt
      python3 get_seq.py neg_set_1.txt uniprot_sprot.fasta > neg_1.fasta
      python3 get_seq.py neg_set_2.txt uniprot_sprot.fasta > neg_2.fasta
### 5. Run HMM Search on Test Sets
      Purpose: Search for Kunitz domains in positive and negative FASTA files using a structural alignment HMM.
      How:
      Runs hmmsearch with the pdb_kunitz.hmm model on each of the four FASTA files (pos_1.fasta, pos_2.fasta, neg_1.fasta, neg_2.fasta).
      Outputs tabular results to .out files (e.g., pos_1_strali.out).
      hmmsearch -Z 1000 --max --tblout pos_1_strali.out pdb_kunitz.hmm pos_1.fasta
      hmmsearch -Z 1000 --max --tblout pos_2_strali.out pdb_kunitz.hmm pos_2.fasta
      hmmsearch -Z 1000 --max --tblout neg_1_strali.out pdb_kunitz.hmm neg_1.fasta
      hmmsearch -Z 1000 --max --tblout neg_2_strali.out pdb_kunitz.hmm neg_2.fasta
### 6. Extract E-values and Build Classification Files
      Purpose: Create files summarizing HMM results for each sequence, labeling positives and negatives.
      How:
      For each .out file, uses awk to:
      Extract the Uniprot ID from the sequence name.
      Assign a label: 1 for positives, 0 for negatives.
      Extract the full sequence E-value ($5) and best domain E-value ($8).
      Write to .class files (e.g., pos_1_strali.class).

      grep -v "^#" pos_1_strali.out | awk '{split($1,a,"|"); print a[2]"\t"1"\t"$5"\t"$8}' > pos_1_strali.class
      grep -v "^#" pos_2_strali.out | awk '{split($1,a,"|"); print a[2]"\t"1"\t"$5"\t"$8}' > pos_2_strali.class
      grep -v "^#" neg_1_strali.out | awk '{split($1,a,"|"); print a[2]"\t"0"\t"$5"\t"$8}' > neg_1_strali.class
      grep -v "^#" neg_2_strali.out | awk '{split($1,a,"|"); print a[2]"\t"0"\t"$5"\t"$8}' > neg_2_strali.class
### 7. Add True Negatives Not Detected by HMM
      Purpose: Ensure all negative sequences are represented, even if not detected by HMM.
      How:
      Finds negative IDs in neg_set_1.txt and neg_set_2.txt that are missing from the corresponding .class files.
      Adds them to the .class files with a fake high E-value (10.0), indicating no detection.

      comm -23 <(sort neg_set_1.txt) <(cut -f 1 neg_1_strali.class | sort) | awk '{print $1"\t"0"\t"10.0"\t"10.0}' >> neg_1_strali.class
      comm -23 <(sort neg_set_1.txt) <(cut -f 1 neg_2_strali.class | sort) | awk '{print $1"\t"0"\t"10.0"\t"10.0}' >> neg_2_strali.class
### 8. Combine Positive and Negative Results
      Purpose: Prepare final datasets for downstream analysis.
      How:
      Concatenates positive and negative .class files to create set_1_strali.class and set_2_strali.class.
      Combines both sets into all_strali.class.

      cat pos_1_strali.class neg_1_strali.class > set_1_strali.class
      cat pos_2_strali.class neg_2_strali.class > set_2_strali.class
      cat set_1_strali.class set_2_strali.class > all_strali.class


## Performance Evaluation 
For this final step the script performance.py script was used to:
- Build a confusion matrix.
- Calculate performance metrics such as:
    q2: overall accuracy
    MCC: Matthews Correlation Coefficient
    TPR: true positive rate (sensitivity/recall)
    PPV: positive predictive value (precision)
- Evaluate the classification performance twice:

  Once using the full-sequence E-value (column 3 in the .class file).
  Once using the best domain E-value (column 4 in the .class file).

This dual evaluation allows to compare whether using the full-sequence E-value or the best domain hit is more accurate for classifying.
Using both approaches increases coverage and robustness but requires careful threshold selection and interpretation to avoid false positives/negatives.

Note: A third argument can also be specified when running the performance.py script:

1 → only evaluate full sequence
2 → only evaluate best domain
0 (or undefined) → evaluate both

performace.py has been run in performance_eval.sh which is filtering the output using grep 'True' and grep 'False' to separate full sequence and single domain results. This works because performance.py prints a line with fullseq= True and another with fullseq= False for each threshold.



## Visualization