tr -d '"' < rcsb_pdb_custom_report_kunitz.csv | \
awk -F ',' '{if (length($2)>0) {name=$2}; if (length($5)>0) {chain=$5}; if (length($4)>0) {seq=$4}; if ($6 ~ /PF00014/ && length(seq)>0) print ">"name"_"chain"\n"seq}' > kunitz_sequences.fasta
