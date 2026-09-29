#!/usr/bin/env python
import gzip
import os
import sys
import Bio
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord
from pathlib import Path


def get_interregions(genbank_path, up_limit, down_limit, min_size):
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
                    intergene_seq = seq_record.seq[last_end:this_start]
                    intergenic_regions[(str(seq_record.name), last_end, this_start)] = (locus_tag, intergene_seq)

            elif i == len(cds_list) - 1 and strand == -1:
                # Last CDS of the record, negative strand
                strand_string = '-'
                this_start = end - down_limit
                last_end = len(seq_record.seq)
                if last_end - this_start > max_size:
                    last_end = this_start + max_size
                if last_end - this_start > min_size:
                    intergene_seq = seq_record.seq[this_start:last_end].reverse_complement()
                    intergenic_regions[(str(seq_record.name), this_start + 1, last_end)] = (locus_tag, intergene_seq)
            elif strand == 1:
                # All other CDSs on the positive strand
                strand_string = '+'
                this_start = start + down_limit
                last_end = cds_list[i-1][1] + 1
                if this_start - last_end > max_size:
                    last_end = this_start - max_size
                if this_start - last_end > min_size:
                    intergene_seq = seq_record.seq[last_end:this_start]
                    intergenic_regions[(str(seq_record.name), last_end, this_start)] = (locus_tag, intergene_seq)
            elif strand == -1:
                # All other CDSs on the negative strand
                strand_string = '-'
                this_start = end - down_limit
                last_end = cds_list[i+1][0]
                if last_end - this_start > max_size:
                    last_end = this_start + max_size
                if last_end - this_start > min_size:
                    intergene_seq = seq_record.seq[this_start:last_end].reverse_complement()
                    intergenic_regions[(str(seq_record.name), this_start + 1, last_end)] = (locus_tag, intergene_seq)
    gbk_fh.close()
    intergenic_records = []
    for key, val in intergenic_regions.items():
        intergenic_records.append(SeqRecord(val[1], id="%s" % (val[0]), description="|%s|%d|%d" % (key[0], key[1], key[2])))
    return intergenic_records


def main():
    in_file = '../genomes.tsv'
    out_dir = './upstreams_nr'
    Path(out_dir).mkdir(exist_ok=True)
    with open(in_file, 'r') as infile:
        infile.readline()
        for line in infile:
            genome_id, _, gbk_path = line.rstrip('\n\r').split('\t')
            intergenic_regions = get_interregions(os.path.join('..',gbk_path), -1000000, 0, 50)
            out_file = os.path.join(out_dir, genome_id + '.upstreams.fna')
            SeqIO.write(intergenic_regions, out_file, "fasta")

    
if __name__ == "__main__":
    main()
