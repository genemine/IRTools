#!/usr/bin/env python

import optparse
import os
import sys

optParser = optparse.OptionParser()
(opts, args) = optParser.parse_args()

EIE_path = args[0]
genome = args[1]

L = 70 # input sequences are required to be 70 bases by labranchor

branchpoint_fa = EIE_path+'branchpoint.fa'
branchpoint_bed = EIE_path+'branchpoint.bed'
branchpoint_predict = EIE_path+'branchpoint_predict.bed'
branchpoint_feature_bed = EIE_path+'branchpoint_feature.bed'
branchpoint_feature_fa = EIE_path+'branchpoint_feature.fa'

def read_file(filename):
    read_data = []
    with open(filename,'r') as fn:
        for l in fn:
            read_data.append(l.strip())
    return read_data

def get_branchpoint_fa():
    EIE_data = read_file(EIE_path+'EIE.tab')
    with open(branchpoint_bed,'w') as b_b:
        for l in EIE_data:
            l = l.split('\t')
            chro = l[0].split('|')[0]
            strand = l[5]
            if strand == '+':
                name = l[0]+':'+l[3]+':'+strand # same as labranchor
                l = [chro,str(int(l[3])-L),l[3],name,'0',strand]
            elif strand == '-':
                name = l[0]+':'+l[2]+':'+strand
                l = [chro,l[2],str(int(l[2])+L),name,'0',strand]
            b_b.write( '\t'.join(l) + '\n')
    os.system('bedtools getfasta -s -name -fo '+branchpoint_fa+' -fi '+genome+' -bed '+branchpoint_bed)
    return EIE_data

def write_branchpoint(EIE_data):
    # re-write branchpoint
    branchpoint_data = read_file(branchpoint_predict)
    branchpoint_feature = []
    with open(EIE_path+'EIE.tab','w') as EIE:
        for i,tl in enumerate(EIE_data):
            tl = tl.split('\t')
            chro = tl[0].split('|')[0]
            name = tl[0]
            strand = tl[5]
            bp_pos = branchpoint_data[i].split('\t')[1] # position of branchpoint
            tl.append(bp_pos)
            EIE.write( '\t'.join(tl) + '\n' )

            # get upstream 5 <--0(bp)--> downstream 5 bases of brachpoint, total 11 nt
            bp_start = int(bp_pos) - 5
            bp_end = int(bp_pos)+1 + 5
            feature = [chro, str(bp_start), str(bp_end), name, '0', strand]
            branchpoint_feature.append(feature)

    with open(branchpoint_feature_bed,'w') as b_f_b:
        for tl in branchpoint_feature:
            b_f_b.write( '\t'.join(tl) + '\n' )
    os.system('bedtools getfasta -s -name -fo '+branchpoint_feature_fa+' -fi '+genome+' -bed '+branchpoint_feature_bed)

    os.system('rm '+branchpoint_bed)
    os.system('rm '+branchpoint_fa)
    os.system('rm '+branchpoint_predict)
    os.system('rm '+branchpoint_feature_bed)

def get_branchpoint():
    EIE_data = get_branchpoint_fa()
    pre_exe_path = os.getcwd()

    pre_EIE_exe = sys.argv[0].split('/')[-1] # current execute file
    labranchor_path = sys.argv[0].replace(pre_EIE_exe,'')+'labranchor/' # go to labranchor execute path
    os.chdir(labranchor_path) # for labranchor execute
    labranchor_path = os.getcwd()+'/' # full path
    os.system('python labranchor.py '+labranchor_path+'2layer.h5 '+'top-bed '+branchpoint_fa+' '+branchpoint_predict)

    os.chdir(pre_exe_path) # recover work_dir
    write_branchpoint(EIE_data)

get_branchpoint()

