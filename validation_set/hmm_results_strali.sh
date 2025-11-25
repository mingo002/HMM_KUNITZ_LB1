#!usr/bin/bash
#STRUCTURAL ALIGNMENT HMM RESULTS
# Run hmmsearch on all positive and negative FASTA files and create tabular output. Generates .out files that contains the E-values computed on the HMM of the sequence alignment.
hmmsearch -Z 1000 --max --tblout pos_1_strali.out pdb_kunitz.hmm pos_1.fasta
hmmsearch -Z 1000 --max --tblout pos_2_strali.out pdb_kunitz.hmm pos_2.fasta
hmmsearch -Z 1000 --max --tblout neg_1_strali.out pdb_kunitz.hmm neg_1.fasta
hmmsearch -Z 1000 --max --tblout neg_2_strali.out pdb_kunitz.hmm neg_2.fasta

# Extract E-values and build .class files with: ID, label (1/0), full sequence E-value ($5), and best domain E-value ($8)
grep -v "^#" pos_1_strali.out | awk '{split($1,a,"|"); print a[2]"\t"1"\t"$5"\t"$8}' > pos_1_strali.class
grep -v "^#" pos_2_strali.out | awk '{split($1,a,"|"); print a[2]"\t"1"\t"$5"\t"$8}' > pos_2_strali.class
grep -v "^#" neg_1_strali.out | awk '{split($1,a,"|"); print a[2]"\t"0"\t"$5"\t"$8}' > neg_1_strali.class
grep -v "^#" neg_2_strali.out | awk '{split($1,a,"|"); print a[2]"\t"0"\t"$5"\t"$8}' > neg_2_strali.class
# Add true negatives not detected by HMM (missing from the .out), with a fake high E-value (10.0)
comm -23 <(sort neg_set_1.txt) <(cut -f 1 neg_1_strali.class | sort) | awk '{print $1"\t"0"\t"10.0"\t"10.0}' >> neg_1_strali.class
comm -23 <(sort neg_set_1.txt) <(cut -f 1 neg_2_strali.class | sort) | awk '{print $1"\t"0"\t"10.0"\t"10.0}' >> neg_2_strali.class

# Combine positive and negative sets into final training/testing files
cat pos_1_strali.class neg_1_strali.class > set_1_strali.class
cat pos_2_strali.class neg_2_strali.class > set_2_strali.class
cat set_1_strali.class set_2_strali.class > all_strali.class