#!/usr/bin/env python

import optparse
import os
import sys

optParser = optparse.OptionParser()
(opts, args) = optParser.parse_args()

EIE_path = args[0]

def read_file(filename):
    read_data = []
    with open(filename,'r') as fn:
        for l in fn:
            read_data.append(l.strip())
    return read_data

def get_seq(sign,e1_len,i_len,el,il,seq):
    # importment param: el and il
    # splice sites(ss): exon1-intron intron-exon2
    # el is the length of base in exon part of ss
    # il is the length of base in intron part of ss
    # seq have been convert into rightly sequence with strand
    f_seq = ''
    if sign == 'ss5':
        f_seq = seq[e1_len-el:e1_len-el+el+il]
    elif sign == 'ss3':
        f_seq = seq[e1_len+i_len-il:e1_len+i_len-il+il+el]
    if 'N' in f_seq:
        f_seq = f_seq.replace('N','A') # 'N' is replaced to 'A'
    return f_seq

def get_ss5_ss3():
    EIE_fa = read_file(EIE_path+'EIE.fa')
    EIE_data = read_file(EIE_path+'EIE.tab')
    with open(EIE_path+'ss5.fa','w') as ss5_f:
        with open(EIE_path+'ss3.fa','w') as ss3_f:
            for i,tl in enumerate(EIE_data):
                tl = tl.split('\t')
                eie_fa = EIE_fa[2*i+1] # have been convert into rightly sequence with strand
                ix1,ix2,ix3,ix4 = int(tl[1]),int(tl[2]),int(tl[3]),int(tl[4])
                strand = tl[5]
                i_len = ix3 - ix2
                if strand == '+':
                    e1_len = ix2 - ix1
                    e2_len = ix4 - ix3
                elif strand == '-':
                    e1_len = ix4 - ix3
                    e2_len = ix2 - ix1

                ex,in5,in3 = 3,6,20 # exon 3 bases,ss5 intron 6 bases,ss3 intron 20 bases
                name = '>'+tl[0]
                ss5_sub_seq = get_seq('ss5',e1_len,i_len,ex,in5,eie_fa)
                ss5_f.write(name+'\n')
                ss5_f.write(ss5_sub_seq+'\n')

                ss3_sub_seq = get_seq('ss3',e1_len,i_len,ex,in3,eie_fa)
                ss3_f.write(name+'\n')
                ss3_f.write(ss3_sub_seq+'\n')

    pre_EIE_exe = sys.argv[0].split('/')[-1] # current execute file
    pl_path = sys.argv[0].replace(pre_EIE_exe,'')+'MaxEntScan/' # go to MaxEntScan execute path
    os.chdir(pl_path) # for MaxEntScan execute

    os.system('perl score5.pl '+EIE_path+'ss5.fa'+' > '+EIE_path+'ss5.out')
    os.system('perl score3.pl '+EIE_path+'ss3.fa'+' > '+EIE_path+'ss3.out')
    ss5 = []
    ss3 = []

    with open(EIE_path+'ss5.out','r') as ss5_o:
        for tl in ss5_o:
            tl = tl.split('\t')
            ss5.append(tl[1])
    with open(EIE_path+'ss3.out','r') as ss3_o:
        for tl in ss3_o:
            tl = tl.split('\t')
            ss3.append(tl[1])

    with open(EIE_path+'EIE.tab','w') as EIE: # overwrite EIE for ss5,ss3
        for i,tl in enumerate(EIE_data):
            tl = tl.split('\t')
            tl.append(ss5[i])
            tl.append(ss3[i])
            EIE.write( '\t'.join( map(str,tl) ) + '\n')

    os.system('rm '+EIE_path+'ss5.fa')
    os.system('rm '+EIE_path+'ss3.fa')
    os.system('rm '+EIE_path+'ss5.out')
    os.system('rm '+EIE_path+'ss3.out')

get_ss5_ss3()
