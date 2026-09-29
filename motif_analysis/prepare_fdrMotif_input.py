#!/usr/bin/env python
import os
import sys
from Bio import SeqIO
from Bio import motifs
from Bio.Seq import Seq
from pathlib import Path
from collections import defaultdict

def autovivify(levels=1, final=dict):
    return (defaultdict(final) if levels < 2 else
            defaultdict(lambda: autovivify(levels - 1, final)))


def get_interregions(genbank_path, up_limit, down_limit, min_size, genome_id):
    intergenic_regions = {}
    if genbank_path.endswith('.gz'):
        gbk_fh = gzip.open(genbank_path, 'rt')
    else:
        gbk_fh = open(genbank_path, 'r')
    max_size = down_limit - up_limit
    for seq_record in SeqIO.parse(gbk_fh, "genbank"):
        cds_list = []
        # Loop over the genome file, get the CDS features 
        for feature in seq_record.features:
            try:
                if feature.type == "CDS":
                    mystart = int(feature.location.start)
                    myend = int(feature.location.end)
                    cds_list.append((mystart, myend, feature.location.strand, feature.qualifiers['locus_tag'][0]))
            except AttributeError as error:
                print(error)
                print('AttributeError', str(feature))
                continue
        for i, (start, end, strand, locus_tag) in enumerate(cds_list):
            if i == 0 and strand == 1:
                # First CDS of the record, positive strand
                strand_string = '+'
                this_start = start + down_limit
                if this_start > max_size:
                    last_end = this_start - max_size
                else:
                    last_end = 0
                if this_start - last_end > min_size:
                    intergene_seq = str(seq_record.seq[last_end:this_start])
                    intergenic_regions[locus_tag] = f'>{genome_id}|{locus_tag}|{str(seq_record.name)}|{str(last_end)}|{str(this_start)}' + '\n' + intergene_seq

            elif i == len(cds_list) - 1 and strand == -1:
                # Last CDS of the record, negative strand
                strand_string = '-'
                this_start = end - down_limit
                last_end = len(seq_record.seq)
                if last_end - this_start > max_size:
                    last_end = this_start + max_size
                if last_end - this_start > min_size:
                    intergene_seq = str(seq_record.seq[this_start:last_end].reverse_complement())
                    intergenic_regions[locus_tag] = f'>{genome_id}|{locus_tag}|{str(seq_record.name)}|{str(this_start + 1)}|{str(last_end)}' + '\n' + intergene_seq
            elif strand == 1:
                # All other CDSs on the positive strand
                strand_string = '+'
                this_start = start + down_limit
                last_end = cds_list[i-1][1] + 1
                if this_start - last_end > max_size:
                    last_end = this_start - max_size
                if this_start - last_end > min_size:
                    intergene_seq = str(seq_record.seq[last_end:this_start])
                    intergenic_regions[locus_tag] = f'>{genome_id}|{locus_tag}|{str(seq_record.name)}|{str(last_end)}|{str(this_start)}' + '\n' + intergene_seq
            elif strand == -1:
                # All other CDSs on the negative strand
                strand_string = '-'
                this_start = end - down_limit
                last_end = cds_list[i+1][0]
                if last_end - this_start > max_size:
                    last_end = this_start + max_size
                if last_end - this_start > min_size:
                    intergene_seq = str(seq_record.seq[this_start:last_end].reverse_complement())
                    intergenic_regions[locus_tag] = f'>{genome_id}|{locus_tag}|{str(seq_record.name)}|{str(this_start + 1)}|{str(last_end)}' + '\n' + intergene_seq
    gbk_fh.close()
    return intergenic_regions


def read_targets(in_file):
    ret = defaultdict(dict)
    with open(in_file, 'r') as infile:
        header = infile.readline().rstrip('\n\r').split('\t')[1:]
        for line in infile:
            row = line.rstrip('\n\r').split('\t')
            if row[2] == '':
                ret[row[0]]['genome'] = 'FW510-R12'
                ret[row[0]]['regulator'] = row[3]
            else:
                ret[row[0]]['genome'] = 'FW104-10B01'
                ret[row[0]]['regulator'] = row[2]
            ret[row[0]]['target'] = row[4]
    return ret


def read_orthologues(ort_files):
    ret = autovivify(2,dict)
    genomes = set()
    for genome, in_file in ort_files.items():
        genomes.add(genome)
        with open(in_file, 'r') as infile:
            header = infile.readline().rstrip('\n\r').split('\t')[1:]
            for line in infile:
                _, species, gene1, gene2 = line.rstrip('\n\r').split('\t')
                ret[genome][gene1][species] = gene2
                ret[genome][gene1][genome] = gene1
                genomes.add(species)
    return ret, sorted(list(genomes))


def read_upstreams(in_dir, genomes):
    ret = defaultdict(dict)
    for genome in genomes:
        fna_path = os.path.join(in_dir, genome + '.upstreams.fna')
        for seq_record in SeqIO.parse(fna_path, "fasta"):
            ret[genome][seq_record.id.split(' ')[0]] = str(seq_record.seq)
    return ret


def main():
    genomes_file = '../genomes.tsv'
    targets_file = '../regulators_targets.txt'
    orthologues_files = {'FW104-10B01':'../orthologs/FW104-10B01_Orthologues.tsv', 'FW510-R12':'../orthologs/FW510-R12_Orthologues.tsv'}
    profiles_dir = '../profiles'
    tf_conservation_file = 'regulator_orthologs.tsv'
    target_conservation_file = 'target_orthologs.tsv'
    min_upstream_size = 80
    up_limit = -350
    down_limit = 0
    
    genome_files = {}
    with open(genomes_file, 'r') as infile:
        infile.readline()
        for line in infile:
            genome_id, _, gbk_path = line.rstrip('\n\r').split('\t')
            genome_files[genome_id] = os.path.join('..', gbk_path)
    
    regulons = read_targets(targets_file)
    orthologues, genomes = read_orthologues(orthologues_files)
    
    with open(tf_conservation_file, 'w') as outfile:
        outfile.write('Genome\t' + '\t'.join(sorted(regulons.keys())) + '\n')
        for genome in genomes:
            outfile.write(genome)
            for regulon in sorted(list(regulons.keys())):
                outfile.write('\t')
                target_genome = regulons[regulon]['genome']
                target_regulator = regulons[regulon]['regulator']
                if genome in orthologues[target_genome][target_regulator]:
                    outfile.write(orthologues[target_genome][target_regulator][genome])
            outfile.write('\n')
    
    target_upstreams = defaultdict(list)
    
    with open(target_conservation_file, 'w') as outfile:
        outfile.write('Genome\t' + '\t'.join(sorted(regulons.keys())) + '\n')
        for genome in genomes:
            print('Processing', genome)
            outfile.write(genome)
            for regulon in sorted(list(regulons.keys())):
                outfile.write('\t')
                target_genome = regulons[regulon]['genome']
                target_regulator = regulons[regulon]['regulator']
                target_genes = regulons[regulon]['target'].split(';')
                if regulon == 'AcrR':
                    down_limit = 10
                else:
                    down_limit = 0
                if genome in orthologues[target_genome][target_regulator]:
                    genes = []
                    upstreams = get_interregions(genome_files[genome], up_limit, down_limit, min_upstream_size, genome)
                    for target_gene in target_genes:
                        if genome in orthologues[target_genome][target_gene]:
                            gene = orthologues[target_genome][target_gene][genome]
                            if gene in upstreams:
                                gene_upstream = upstreams[gene]
                                genes.append(gene)
                                target_upstreams[regulon].append(gene_upstream)
                    if genes:
                        outfile.write(';'.join(genes))
            outfile.write('\n')
    for regulon in regulons:
        output_dir = regulon
        Path(output_dir).mkdir(exist_ok=True)
        with open(os.path.join(output_dir, regulon + '.target_upstreams.fna'), 'w') as outfile:
            outfile.write('\n\n'.join(target_upstreams[regulon]))
            outfile.write('\n')
        profile_file = os.path.join(profiles_dir, regulon + '.profile')
        
        sites = []
        with open(profile_file, 'r') as infile:
            for line in infile:
                sites.append(Seq(line.rstrip('\n\r').upper()))
        m = motifs.create(sites)
        pwm = m.counts.normalize(pseudocounts={"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25})
        out_file = os.path.join(output_dir, regulon + '.mx')
        with open(out_file, 'w') as outfile:
            outfile.write('4  ' + str(len(m)) + '\n')
            for line in str(pwm).split('\n')[1:]:
                outfile.write(line[2:] + '\n')


if __name__ == "__main__":
    main()
