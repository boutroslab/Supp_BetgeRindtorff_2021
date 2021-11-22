#!/usr/bin/python

import os
from pyliftover import LiftOver

manifest = snakemake.input[0]
#Chain file downloaded from: http://hgdownload.cse.ucsc.edu/gbdb/hg19/liftOver/
chain_file = "/Users/valentini/Documents/amplicon_analysis/hg19ToHg38.over.chain.gz"
# op = open('conversion.txt', 'w')
# op.write('sample_id'+'\tsample'+'\tchrom_h19'+'\tpos_h19'+'\tref_h19'+'\talt_h19'\
# +'\tconversion_h38'+'\n')
manifest_name_output = snakemake.output[0]
manifest_output = snakemake.output[1]
man_names = open(manifest_name_output, 'w')
man = open(manifest_output, 'w')

lo = LiftOver(chain_file)
amplicon_list = []

with open(manifest) as f:
    for line in f:
        if "chr" in line.split("\t")[5]:
            TARGET = line.split("\t")[0]
            CHROM = line.split("\t")[5]
            START = line.split("\t")[6]
            STOP = line.split("\t")[7]
            STRAND = line.split("\t")[8]
            print(CHROM, START, STOP)
            print("Converted to:")
            conversion_start = lo.convert_coordinate(CHROM, int(START), STRAND)
            con_start = conversion_start[0]
            CHROM_H38=con_start[0]
            START_H38=con_start[1]
            conversion_stop = lo.convert_coordinate(CHROM, int(STOP), STRAND)
            con_stop = conversion_stop[0]
            STOP_H38=con_stop[1]
            print(TARGET,CHROM_H38,":",START_H38,"-",STOP_H38)
            man_names.write(TARGET+"\t"+CHROM_H38+":"+str(START_H38)+"-"+str(STOP_H38)+"\n")
            amplicon = CHROM_H38+":"+str(START_H38)+"-"+str(STOP_H38)
            if amplicon not in amplicon_list:
                amplicon_list.append(amplicon)
                man.write(CHROM_H38+":"+str(START_H38)+"-"+str(STOP_H38)+"\n")

man_names.close()
man.close()
