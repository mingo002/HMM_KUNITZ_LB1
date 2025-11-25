#PERFORMANCES EVALUATION
#For each class file (set_1_strali.class, set_2_strali.class), and for each coverage type ("True" = full sequence, "False" = one domain):
#-Runs performance evaluation for thresholds 1e-1 to 1e-12. 
#-Filters results by coverage type. (full sequence or one domain)
#-Sorts by the performance metric (column 6).
#-Selects the best threshold. (keeps only the top line)
#-Prints the result with a label.

for i in $(seq 1 12); do
    python3 performance.py set_1_strali.class 1e-"$i"
done | grep 'threshold' | grep 'True' | sort -nrk 6 | head -n 1 | awk '{print "set_1_strali.class (full_seq):", $2}'

for i in $(seq 1 12); do
    python3 performance.py set_1_strali.class 1e-"$i"
done | grep 'threshold' | grep 'False' | sort -nrk 6 | head -n 1 | awk '{print "set_1_strali.class (one_domain):", $2}'

for i in $(seq 1 12); do
    python3 performance.py set_2_strali.class 1e-"$i"
done | grep 'threshold' | grep 'True' | sort -nrk 6 | head -n 1 | awk '{print "set_2_strali.class (full_seq):", $2}'

for i in $(seq 1 12); do
    python3 performance.py set_2_strali.class 1e-"$i"
done | grep 'threshold' | grep 'False' | sort -nrk 6 | head -n 1 | awk '{print "set_2_strali.class (one_domain):", $2}'

