
'''THIS SCRIPT EXTRACT PDB ID AND SEQUENCE FROM FASTA FILES '''
#!/usr/bin/python

import sys

def get_ids(idlist):
    f = open(idlist)
    return f.read().strip().split('\n')

def get_sequence(pidlist,seqfile):
    f = open(seqfile)
    for line in f: 
        if line.startswith('>'): 
           pid = line.split('|')[1].strip()
           print_seq = pid in pidlist
        if print_seq:
            print(line.rstrip())


if __name__ == '__main__':
    idlist = sys.argv[1]
    seqfile = sys.argv[2]
    pidlist = get_ids(idlist)
    get_sequence(pidlist,seqfile)