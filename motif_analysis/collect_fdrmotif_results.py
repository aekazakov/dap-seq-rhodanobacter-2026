#!/usr/bin/env python
import os

def main():
    regulons_file = '../regulators_targets.txt'
    out_file = 'fdrMotif_results.txt'

    with open(out_file, 'w') as outfile:
        outfile.write('\t'.join(['Regulon', 'Motif consensus', 'FDR upper limit', 'Computed FDR', 'Sites', 'E-value']) + '\n')
        with open(regulons_file, 'r') as listfile:
            listfile.readline()
            for line in listfile:
                row = line.rstrip('\n\r').split('\t')
                regulon = row[0]
                in_file = os.path.join(regulon, regulon + '-motif.summary.txt')
                if os.path.exists(in_file):
                    with open(in_file, 'r') as infile:
                        row = infile.readline().split('\t')
                        outfile.write(regulon + '\t' + row[1] + '\t' + '\t'.join(row[3:]))
                else:
                    outfile.write(regulon + '\tN/A\n')
if __name__ == "__main__":
    main()
