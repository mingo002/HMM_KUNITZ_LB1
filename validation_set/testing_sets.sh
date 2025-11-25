# Extract Uniprot IDs of sequences with high identity (≥95%) and alignment length ≥50 to remove them from the training/testing pool
grep -v "^#" bpti_kunitz_clean.blast | awk '{if ($3>=95 && $4>=50) print $2}' | sort -u | cut -d "|" -f 2 > ids_to_remove.txt

# Extract all IDs from the full Kunitz dataset
grep ">" all_kunitz.fasta | cut -d "|" -f 2 > all_kunitz_ids.txt

echo 'Creating positive test files (set1 and set2)...'

# Remove the overlapping IDs and keep only the remaining Kunitz proteins
comm -23 <(sort all_kunitz_ids.txt) <(sort ids_to_remove.txt) > ids_to_keep.txt
sort -R ids_to_keep.txt > kunitz_shuffled.txt
head -n 184 kunitz_shuffled.txt > pos_set_1.txt
tail -n 184 kunitz_shuffled.txt > pos_set_2.txt
python3 get_seq.py pos_set_1.txt uniprot_sprot.fasta > pos_1.fasta
python3 get_seq.py pos_set_2.txt uniprot_sprot.fasta > pos_2.fasta

echo 'Creating negative test files (set1 and set2)...'

# Extract all Swiss-Prot IDs
grep ">" uniprot_sprot.fasta | cut -d "|" -f 2 > all_sp_ids.txt

# Remove the Kunitz proteins to get the negative candidates
comm -23 <(sort all_sp_ids.txt) <(sort all_kunitz_ids.txt) > negative_ids.txt
sort -R negative_ids.txt > negative_shuffled.txt
head -n 286286 negative_shuffled.txt > neg_set_1.txt
tail -n 286286 negative_shuffled.txt > neg_set_2.txt
python3 get_seq.py neg_set_1.txt uniprot_sprot.fasta > neg_1.fasta
python3 get_seq.py neg_set_2.txt uniprot_sprot.fasta > neg_2.fasta

echo 'FASTA files created.'
