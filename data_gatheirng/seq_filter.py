from Bio import SeqIO

# Parameters
input = "kunitz_cluster.txt"
output = "filtered_kunitz.fasta"
bin_seq = "filtered_out.txt"

min_len = 40    # minimum length of sequence
max_len = 100   # maximum length of sequence
allowed_chars = set("ACDEFGHIKLMNPQRSTVWY-")  # standard amino acids + gap


kept = []
removed = []

for record in SeqIO.parse(input, "fasta"):
    seq_str = str(record.seq).upper()
    seq_len = len(seq_str)
    
    # conditions
    too_short = seq_len < min_len
    too_long = seq_len > max_len
    invalid_chars = False
    for c in seq_str:
        if c not in allowed_chars:
            invalid_chars = True
            break

    if too_short or too_long or invalid_chars:
        reason = []
        if too_short: reason.append(f"too short ({seq_len})")
        if too_long: reason.append(f"too long ({seq_len})")
        if invalid_chars: reason.append("invalid chars")
        removed.append((record.id, ", ".join(reason)))
    else:
        kept.append(record)

# write filtered fasta
SeqIO.write(kept, output, "fasta")

# write log
with open(bin_seq, "w") as bin:
    for rid, reason in removed:
        bin.write(f"{rid}\t{reason}\n")

print(f"Kept {len(kept)} sequences -> {output}")
print(f"Removed {len(removed)} sequences -> {bin_seq}")
