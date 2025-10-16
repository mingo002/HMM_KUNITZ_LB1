#!/usr/bin/env bash
    
    INPUT="rcsb_pdb_custom_report_kunitz.csv"
    OUT="pdb_kunitz_customreported.fasta"
    
    tail -n +2 "$INPUT" | tr -d '"' \
    | awk -F ',' '
        $7 ~ /PF00014/ {
            id = $NF   # Always get the last column (Entry Id)
            seq = $4   # Sequence
    
            if (id == "" || seq == "") next  # skip if missing
    
            print ">" id
            for (i=1; i<=length(seq); i+=80)
                print substr(seq,i,80)
        }' > "$OUT"
    

echo "rote FASTA to $OUT successfully!"






