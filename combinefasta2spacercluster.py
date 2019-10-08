import sys
fasta = {}
with open(sys.argv[1]) as file_one:
    for line in file_one:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            active_sequence_name = line[5:]
            if active_sequence_name not in fasta:
                fasta[active_sequence_name] = []
            continue
        sequence = line
        fasta[active_sequence_name].append(sequence)
#print fasta
with open(sys.argv[2]) as spacer2cluster:
    for line in spacer2cluster:
	if line.split('\t')[0] in fasta:
	    print line.strip()+'\t'+str(fasta[line.split('\t')[0]][0])
#print fasta
