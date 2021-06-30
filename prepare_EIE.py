#!/usr/bin/env python

import os
import sys
import optparse

optParser = optparse.OptionParser()
(opts, args) = optParser.parse_args()
genome = args[0]
gtf = args[1]
EIE_path = args[2]
branchpoint_flag = args[3]

if not os.path.isdir(EIE_path):
    os.makedirs(EIE_path)

### step1: get EIE_primer ###
def merge(features):
    features.sort( key = lambda x:(x[1],x[2]) )
    start = 0
    lenf = len(features)
    merge_feature = []
    while start < lenf:
        chro = features[start][0]
        e_s = features[start][1]
        e_e = features[start][2]
        strand = features[start][3]
        if start == lenf-1: # last one
            merge_feature.append( [chro,e_s,e_e,strand] )
            break
        for i,val in enumerate(features[start+1:lenf]):
            start += 1
            if val[1] > e_e:
                merge_feature.append( [chro,e_s,e_e,strand] )
                break
            else:
                e_e = max(e_e,val[2]) # the key of merge function
            if start == lenf-1: # last one, this is the condition of exon-exon, merge
                start = lenf
                merge_feature.append( [chro,e_s,e_e,strand] )
                break
    return merge_feature

def unique_intron(start,end,origin_gene_chro):
    count = 0
    overlap_gene = ''
    strand = ''
    for gene in origin_gene_chro:
        if start >= gene[1] and end <= gene[2]:
            count += 1
            overlap_gene = gene[4]
            strand = gene[3]
            if count >= 2:
                break
    if count == 1:
        return overlap_gene,strand
    else:
        return 'None',strand

def get_intron(gtf,chrosomes=list(map(str,range(1,23)))+['X','Y']):
    print('Get intron from gtf')
    gene_coor = {}
    exon_coor = {}
    intron_coor = {}
    gene_order = []
    ### get gene and exon feature
    f = open(gtf)
    for line in f:
        tl = line.split('\t')
        if line[0] != '#' and tl[0] in chrosomes:  # skip comment line and insure standard chrosome
            chro = tl[0]
            gene_id  = line.split('gene_id')[1].split('"')[1]
            if tl[2] == 'gene':
                gene_tmp = [chro,int(tl[3])-1,int(tl[4]),tl[6],gene_id] # 1-based convert to 0-based
                if chro in gene_coor:
                    gene_coor[chro].append(gene_tmp)
                else:
                    gene_coor[chro] = [ gene_tmp ]
            elif tl[2] == 'exon':
                exon_tmp = [chro,int(tl[3])-1,int(tl[4]),tl[6]]
                if chro in exon_coor:
                    exon_coor[chro].append(exon_tmp)
                else:
                    exon_coor[chro] = [ exon_tmp ]
    f.close()

    ## get intron feature
    for chro in chrosomes:
        if chro not in gene_coor.keys():
            continue
        merge_gene_chro = merge(gene_coor[chro])
        merge_exon_chro = merge(exon_coor[chro])
        for gene in merge_gene_chro:
            exon_in_gene = [] # all exon in a single gene
            for exon in merge_exon_chro:
                if (exon[1] >= gene[1] and exon[1] < gene[2]) == True or (exon[2] > gene[1] and exon[2] < gene[2]) == True:
                    if exon[1] <= gene[1]: # recover exon
                        exon[1] = gene[1]
                    if exon[2] >= gene[2]:
                        exon[2] = gene[2]
                    exon_in_gene.append(exon)
            len_exon = len(exon_in_gene)
            if len_exon < 2:
                continue
            # get intron feature: upstream_exon (e1) - intron (i) - downstream_exon (e2)
            for i in range(len_exon-1): # get intron, foreach each gene
                e1 = exon_in_gene[i]
                e2 = exon_in_gene[i+1]
                ix1 = e1[1] # up_exon_start
                ix2 = e1[2] # up_exon_end
                ix3 = e2[1] # down_exon_start
                ix4 = e2[2] # down_eoxon_end
                if ix3-ix2>=20 and ix2-ix1>=3 and ix4-ix3>=3: # for MaxEntScan,exon:3base,intron:20bp
                    gene_id,strand = unique_intron( ix2, ix3, gene_coor[chro] ) # origin gene_coor_chro with gene_id
                    if gene_id != 'None': # unique intron with gene_id != 'None',otherwise 'None'
                        intron_tmp = [ e1[0],ix1,ix2,ix3,ix4,strand,gene_id ]
                        if gene_id not in gene_order: # gene added order
                            gene_order.append(gene_id)
                        if gene_id in intron_coor:
                            intron_coor[gene_id].append(intron_tmp)
                        else:
                            intron_coor[gene_id] = [ intron_tmp ]
    return intron_coor,gene_order

def get_gene_type(gtf_file,chrosomes=list(map(str,range(1,23)))+['X','Y']):
    gtf = open(gtf_file)
    gene_type_dict = {} # for gene_biotype, if protein_coding gene:2, lincRNA gene:1, otherwise gene:0
    for line in gtf:
        tl = line.split('\t')
        if line[0] != '#' and tl[0] in chrosomes and tl[2] == 'gene': # skip comment line and insure gene within standard chrosome
            gene_id  = line.split('gene_id')[1].split('"')[1]
            gene_biotype = line.split('gene_biotype')[1].split('"')[1]
            gene_start = int(tl[3])
            gene_end = int(tl[4])
            if gene_biotype == 'protein_coding':
                gene_type_dict[gene_id] = [ 2, gene_start, gene_end ]
            elif gene_biotype == 'lincRNA':
                gene_type_dict[gene_id] = [ 1, gene_start, gene_end ]
            else:
                gene_type_dict[gene_id] = [ 0, gene_start, gene_end ]
    gtf.close()
    return gene_type_dict

def get_EIE_primer(gtf):
    intron_coor,gene_order = get_intron(gtf)
    gene_type_dict = get_gene_type(gtf)
    EIE_primer = []
    ie_ratio = 0.5 # threshold for short intron
    ie_ratio_middle = 1.0 # threshold for middle intron
    for gene_id in gene_order:
        gene_intron_data = intron_coor[gene_id]
        gene_type = gene_type_dict[gene_id][0]
        gene_start = gene_type_dict[gene_id][1]
        gene_end = gene_type_dict[gene_id][2]
        gene_length = gene_end - gene_start
        num = len(gene_intron_data)
        for j,i_f in enumerate(gene_intron_data):
            chro,ix1,ix2,ix3,ix4,strand = i_f[0],i_f[1],i_f[2],i_f[3],i_f[4],i_f[5]
            i_len = ix3-ix2
            e_len = float( (ix2-ix1) + (ix4-ix3) )
            intron_id_meta = [chro,ix1,ix2,ix3,ix4,strand,gene_id]
            intron_id = '|'.join( map(str,intron_id_meta) )
            if i_len / e_len > ie_ratio_middle:
                intron_type = 'A' # long intron
            elif i_len / e_len > ie_ratio:
                intron_type = 'B' # middle intron
            else:
                intron_type = 'C' # short intron
            if strand == '+':
                # relative intron start position in gene, ranging from 0 to 1
                rel_i_start_pos = float(ix2-gene_start) / gene_length
                # relative intron end position in gene, ranging from 0 to 1
                rel_i_end_pos =  float(ix3-gene_start) / gene_length
            else:
                rel_i_start_pos = float( abs(ix3-gene_end) ) / gene_length
                rel_i_end_pos =  float( abs(ix2-gene_end) ) / gene_length
            if strand == '+':
                i_rank = j + 1
            else:
                i_rank = num - j
            eie = [ intron_id,ix1,ix2,ix3,ix4,strand,intron_type,i_rank,num,rel_i_start_pos,rel_i_end_pos,gene_type ]
            # j is jth intron in the gene;'num' is the number of intron in the gene
            EIE_primer.append(eie)
    return EIE_primer

EIE_primer = get_EIE_primer(gtf)
### step1 end: get EIE_primer ###

### step2: get EIE and EIE_fa ###
def get_EIE_and_EIE_fa(EIE_primer):
    print('Get EIE bed file and EIE fasta file')
    with open(EIE_path+'EIE.tab','w') as EIE:
        for eie in EIE_primer:
            EIE.write( '\t'.join( map(str,eie) ) + '\n' )

    with open(EIE_path+'EIE.bed','w') as EIE_bed:
        for eie in EIE_primer:
            chro = eie[0].split('|')[0]
            l = [chro,eie[1],eie[4],eie[0],'0',eie[5]]
            EIE_bed.write( '\t'.join( map(str,l) ) + '\n' )

    os.system('bedtools getfasta -s -name -fo '+EIE_path+'EIE.fa'+' -fi '+genome+' -bed '+EIE_path+'EIE.bed')

get_EIE_and_EIE_fa(EIE_primer)

# MaxEntScan
os.system('get_max_ent_ss5_ss3.py '+EIE_path) # appended column index of EIE.tab: -3 and -2

if branchpoint_flag == 'off':
    pass # not get feature of branchpoint
else:
    # branchpoint_flag == 'on'
    # branchpoint fa -> L = 70
    print('Get branchpoint')
    os.system('get_branchpoint.py '+EIE_path+' '+genome) # appended column index of EIE.tab: -1

# build model intron feature
os.system('build_model_intron_feature.py '+EIE_path+' '+branchpoint_flag) # get static model feature
### step2 end: get EIE and fa###

print('EIE model is built successful!')
