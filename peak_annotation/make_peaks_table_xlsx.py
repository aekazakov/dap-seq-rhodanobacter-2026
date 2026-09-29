#!/usr/bin/env python
import os
import csv
import gzip
from Bio import SeqIO
from openpyxl import Workbook
from collections import defaultdict

def autovivify(levels=1, final=dict):
    return (defaultdict(final) if levels < 2 else
            defaultdict(lambda: autovivify(levels - 1, final)))

def read_genes(genome, gbff_file):
    '''
        Reads a genome file in GenBank format, makes a dictionary of genes,
        finds UTR coordinates
    '''
    result = defaultdict(dict)
    with open(gbff_file, 'r') as infile:
        for seq_record in SeqIO.parse(infile, "genbank"):
            genes = []
            products = {}
            for feature in seq_record.features:
                if feature.type == 'gene':
                    gene_id = str(feature.qualifiers['locus_tag'][0])
                    genes.append([gene_id, feature.location.strand, feature.location.start, feature.location.end, seq_record.id])
                else:
                    if 'product' in feature.qualifiers and 'locus_tag' in feature.qualifiers:
                        products[str(feature.qualifiers['locus_tag'][0])] = str(feature.qualifiers['product'][0])
                    
            for gene_ind, gene in enumerate(genes):
                gene_id = gene[0]
                result[gene_id]['strand'] = gene[1]
                result[gene_id]['start'] = gene[2]
                result[gene_id]['end'] = gene[3]
                result[gene_id]['contig'] = gene[4]
                if gene_id in products:
                    result[gene_id]['product'] = products[gene_id]
                else:
                    result[gene_id]['product'] = ''
                if gene[1] == 1:
                    utr_end = gene[2] - 1
                    if gene_ind == 0:
                        utr_start = 1
                    else:
                        utr_start = genes[gene_ind - 1][3] + 1
                elif gene[1] == -1:
                    utr_start = gene[3] + 1
                    if gene_ind == len(genes) - 1:
                        utr_end = len(seq_record.seq)
                    else:
                        utr_end = genes[gene_ind + 1][2] - 1
                else:
                    print(gene_id, 'strand not defined:', gene[1])
                result[gene_id]['utr_start'] = utr_start
                result[gene_id]['utr_end'] = utr_end
    return result


def read_orthologs(ortholog_file, genomes):
    orthologs = defaultdict(list)
    with open(ortholog_file, 'r') as infile:
        infile.readline()
        for line in infile:
            row = line.rstrip('\n\r').split('\t')
            if row[1] in genomes:
                for gene_id in row[3].split(', '):
                    if gene_id != '':
                        orthologs[row[2]].append(gene_id)
                        orthologs[gene_id].append(row[2])
    return orthologs
    

def read_tf_list(tf_file, orthologs):
    result = {}
    with open(tf_file, 'r') as infile:
        for line in infile:
            tf = line.rstrip('\n\r')
            if tf in orthologs:
                result[tf] = orthologs[tf]
            else:
                result[tf] = []
    return result
    

def annotate_peak(contig, start, end, summit, genes):
    '''
    Finds genes around a peak. 
    Returns list of genes, if the summit of the peak is in the 
    gene's UTR or less than 50 bp from the UTR
    '''
    result = []
    for gene_id, gene_data in genes.items():
        if gene_data['contig'] != contig:
            continue
        if summit > gene_data['utr_start'] - 50:
            if summit < gene_data['utr_end'] + 50:
                result.append([gene_id, gene_data['product']])
    return result

    
def read_peaks_file(peaks_file, genes):
    '''
        Reads peaks file and populates the peaks dictionary
    '''
    peaks = []
    header_row = []
    with open(peaks_file, 'r') as infile:
        for line in infile:
            line = line.rstrip('\n\r')
            if line == '' or line.startswith('#'):
                continue
            elif line.startswith('chr\t'):
                header_row = line.rstrip('\n\r').split('\t')
                #header_row.append('Genes')
                #peaks.append(header_row)
            else:
                row = line.split('\t')
                locus_tags = []
                for (locus_tag, product) in annotate_peak(row[0], int(row[1]), int(row[2]), int(row[4]), genes):
                    locus_tags.append(locus_tag)
                row[1] = int(row[1])
                row[2] = int(row[2])
                row[3] = int(row[3])
                row[4] = int(row[4])
                row[5] = int(row[5])
                row[6] = float(row[6])
                row[7] = float(row[7])
                row[8] = float(row[8])
                row.append(';'.join(locus_tags))
                peaks.append(row)
    return peaks
    

def main():
    # Genome names
    genome1 = 'FW104-10B01'
    genome2 = 'FW510-R12'

    # The name of a file with orthologous genes calculated for genome1 and genome2
    ortholog_file = '../orthologs/FW104-10B01_Orthologues.tsv'
    # Location of genome files
    genomes = {'FW104-10B01':'../genomes/FW104-10B01.ncbi.1.genome.gbff', 'FW510-R12':'../genomes/FW510-R12.ncbi.1.genome.gbff'}

    # Read genes 
    genes = {}
    for genome in genomes:
        genes[genome] = read_genes(genome, genomes[genome])
    # Read the table of orthologs
    orthologs = read_orthologs(ortholog_file, genomes)
            
    xlsx_outfile = 'DAP-seq_peaks_dataset.xlsx'
    wb = Workbook()
    ws = wb.active
    ws.title = 'DAP-seq_peaks'
    ws.append(['Protein_source', 'TF locus tag', 'Genomic DNA','chr','start','end','length','abs_summit','pileup','-log10(pvalue)','fold_enrichment','-log10(qvalue)','Peak_name','Genes'])
    for genome_tf in genomes:
        for genome_dna in genomes:
            # Read lists of TFs
            tfs = read_tf_list(genome_tf + '_TFs.txt', orthologs)

            # Files with combined statistics, tab-separated
            summary_file = genome_dna + '_dap_stats.tsv'

            good_samples = []
            with open(summary_file, 'r') as infile:
                header_row = infile.readline().rstrip('\n\r').split('\t')
                for line in infile:
                    row = line.rstrip('\n\r').split('\t')
                    if row[-2] in tfs:
                        if row[7] != '0':
                            peaks_file = os.path.join(genome_dna + '_peaks', row[0] + '_peaks.xls')
                            for peak_row in read_peaks_file(peaks_file, genes[genome_dna]):
                                ws.append([genome_tf, row[-2], genome_dna] + peak_row)
    ws.column_dimensions['A'].width = 15
    ws.column_dimensions['N'].width = 60
    wb.save(xlsx_outfile)
    

if __name__=='__main__':
    main()
