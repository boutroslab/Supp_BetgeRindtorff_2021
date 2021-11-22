#!/usr/bin/python

import os

vcf_directory = snakemake.input[0]
op = open(snakemake.output[0], 'w')
op.write('SAMPLE'+'\tCHROM'+'\tPOS'+'\tID'+'\tREF'+'\tALT'\
+'\tQUAL'+'\tFILTER')

GATK_fields = ['BaseQRankSum', 'ECNT', 'FS', 'HCNT', 'MAX_ED', 'MIN_ED', 'TLOD', 'PON',
'GT', 'AD', 'AF']

for field in GATK_fields:
    op.write('\t' + field)

VEP_fields = ['Feature', 'Feature_type', 'Consequence', 'Protein_position', \
'Amino_acids', 'Codons', 'Existing_variation', 'IMPACT', 'SYMBOL', 'BIOTYPE', \
'SIFT', 'PolyPhen', 'EXON', 'gnomAD_AF', 'CLIN_SIG', 'PUBMED']

# unwanted_consequences = ['downstream_gene_variant', 'upstream_gene_variant',\
# 'synonymous_variant', 'non_coding_transcript_variant', '3_prime_UTR_variant', \
# 'intron_variant', '3_prime_UTR_variant', 'non_coding_transcript_exon_variant']

for field in VEP_fields:
    op.write('\t' + field)

op.write('\n')

for root, dirs, files in os.walk(vcf_directory):
    for name in files:
        if "mutect2_pon_qc_VEP.vcf" in name and "html" not in name:
            name = os.path.join(root,name)
            print(name)
            with open(name) as f:
                for line in f:
                    if line[0] != "#":
                        gatk_list = ['-'] * len(GATK_fields)
                        #Define all fields:
                        SAMPLE = name.split('_')[0].split('/')[1]
                        field = line.split("\t")
                        CHROM = field[0]
                        POS = field[1]
                        ID = field[2]
                        REF = field[3]
                        ALT = field[4]
                        QUAL = field[5]
                        FILTER = field[6]
                        INFO = field[7]
                        FORMAT = field[8]
                        DATA = field[9]
                        #Add info to gatk list (last 3 fields in FORMAT)
                        for index, data in enumerate(DATA.split(":")):
                            if FORMAT.split(":")[index] in GATK_fields:
                                if data != '.':
                                    gatk_list[GATK_fields.index(FORMAT.split(":")[index])] = data.strip()
                        #Split info in 2: gatk and vep
                        if "CSQ=ENST" in INFO:
                            gatk_info = INFO.split("CSQ=")[0]
                            #Add gatk info  to list
                            for gatk in gatk_info.split(";"):
                                if gatk.split('=')[0] in GATK_fields:
                                    if gatk.split('=')[1].strip() != '.':
                                        gatk_list[GATK_fields.index(gatk.split('=')[0])] = gatk.split('=')[1].strip()
                            vep_info = INFO.split("CSQ=")[1]
                            #Add vep info to list
                            for vep in vep_info.split(","):
                                vep_list = ['-'] * len(VEP_fields)
                                for index, elem in enumerate(vep.split("|")):
                                    if elem != '':
                                        if '&' in elem:
                                            elem = elem.replace('&', ',')
                                        vep_list[index] = elem
                                #Print all the info in output file
                                op.write(SAMPLE + '\t' + CHROM + '\t' + POS + '\t'+ ID + '\t' + REF + '\t' + ALT + '\t' + QUAL + '\t' + FILTER)
                                for elem in gatk_list:
                                    op.write('\t' + str(elem))
                                for elem in vep_list:
                                    op.write('\t' + str(elem))

                                op.write('\n')
                        #If the info part has a different structure then you will have to do something else
                        else:
                            print(gatk_info)
                            op.write(SAMPLE + '\t' + CHROM + '\t' + POS + '\t'+ ID + '\t' + REF + '\t' + ALT + '\t' + QUAL + '\t' + FILTER)
                            for elem in gatk_list:
                                op.write('\t' + str(elem))
                            op.write('\n')
