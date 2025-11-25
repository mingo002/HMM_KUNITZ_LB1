# Usage: bash run_mcc_fullseq.sh INPUT_FILE.txt

INPUT=$1
OUTFILE="best_dom_graph.tsv"

echo -e "threshold\tmcc" >> $OUTFILE

for i in $(seq 1 12)
do
    th="1e-$i"
    MCC=$(python performance.py "$INPUT" $th 2 | grep "MCC=" | head -n 1 | sed -E 's/.*MCC= ([^ ]+).*/\1/')
    echo -e "$th\t$MCC" >> $OUTFILE
done

echo "Results saved to $OUTFILE"
