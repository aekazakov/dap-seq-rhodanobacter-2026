#!/usr/bin/env python
import os
from pathlib import Path
from Bio import SeqIO

def main():
    genomes_file = '../genomes.tsv'
    in_file = '../dereplicate_genome_set/Wdb.csv'
    out_dir = './proteins'
    Path(out_dir).mkdir(exist_ok=True)
    
    genomes = {}
    with open(genomes_file, 'r') as infile:
        infile.readline()
        for line in infile:
            genome_id, _, gbk_path = line.rstrip('\n\r').split('\t')
            genomes[genome_id] = os.path.join('..', gbk_path)
    
    with open(in_file, 'r') as infile:
        infile.readline()
        for line in infile:
            fasta_name, _, _ = line.rstrip('\n\r').split(',')
            genome_id = fasta_name.split('.upstreams')[0]
            gbk_path = genomes[genome_id]
            out_file = os.path.join(out_dir, genome_id + '.faa')
            if os.path.exists(out_file):
                print(out_file, 'already exists')
                continue
            feature_index = 0
            with open(out_file, 'w') as outfile:
                with open(gbk_path, 'r') as infile:
                    for seq_record in SeqIO.parse(infile, "genbank"):
                        for feature in seq_record.features:
                            if feature.type == 'CDS':
                                if 'translation' in feature.qualifiers:
                                    if 'locus_tag' in feature.qualifiers:
                                        outfile.write('>' + genome_id + '|' + str(feature.qualifiers['locus_tag'][0]) + '\n')
                                    elif 'protein_id' in feature.qualifiers:
                                        outfile.write('>' + genome_id + '|' + str(feature.qualifiers['protein_id'][0]) + '\n')
                                    else:
                                        feature_index += 1
                                        outfile.write('>' + genome_id + '|' + 'gene_' + str(feature_index) + '\n')
                                    outfile.write(''.join(feature.qualifiers['translation']) + '\n')
            
if __name__=='__main__':
    main()

