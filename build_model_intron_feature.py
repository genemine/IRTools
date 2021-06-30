#!/usr/bin/env python

import os
import math
import optparse
from extract_intron_feature_tools import read_file,my_round,compute_GC_content,entroy,compute_dinucleotide_content,splice_site_type,get_right_seq_bp

optParser = optparse.OptionParser()
(opts, args) = optParser.parse_args()
EIE_path = args[0]
branchpoint_flag = args[1]

decimal = 2 # keep 2 decimal
###                            ###
### prepare for intron feature ###
###                            ###
EIE_fa = read_file(EIE_path+'EIE.fa')
if branchpoint_flag == 'off':
    pass # not get feature of branchpoint
else:
    # branchpoint_flag == 'on'
    bp_feature_fa = read_file(EIE_path+'branchpoint_feature.fa')

def extract_model_feature():
    #######                 #######
    #######     feature     #######
    ######                  #######
    feature_data = []
    EIE_tab = open(EIE_path+'EIE.tab')
    for i,eie in enumerate(EIE_tab): # get feature of intron
        eie = eie.strip().split('\t')
        eie_fa = EIE_fa[2*i+1] # sequence of each EIE, have been convert into rightly sequence with strand
        ix1,ix2,ix3,ix4 = int(eie[1]),int(eie[2]),int(eie[3]),int(eie[4])
        strand,i_type_len,i_rank,num = eie[5],eie[6],int(eie[7]),int(eie[8])
        tl = [] # temp storing feature of each eie
        tl.append( eie[0] ) # id

        i_len = ix3 - ix2 # 1 feature
        if strand == '+':
            e1_len = ix2 - ix1 # 1 feature
            e2_len = ix4 - ix3 # 1 feature
        elif strand == '-':
            e1_len = ix4 - ix3
            e2_len = ix2 - ix1
        type_i = 0 # 1 feature
        if eie[6] == 'C': type_i = 2 # short intron
        elif eie[6] == 'B': type_i = 1 # middle intron
        else: type_i = 0 # long intron
        rel_i_rank =  my_round( float(i_rank) / num, decimal )  # 1 feature
        first_i = 1 if i_rank==1 else 0 # 1 feature
        middle_i = 1 if i_rank!=1 and i_rank!=num else 0 # 1 feature
        last_i = 1 if i_rank==num else 0 # 1 feature
        tl.append( i_len )
        tl.append( e1_len )
        tl.append( e2_len )
        tl.append( type_i )
        tl.append( rel_i_rank )
        tl.append( first_i )
        tl.append( middle_i )
        tl.append( last_i )
        # 8 feature above

        rel_e1i_len = float(i_len) / e1_len # 1 feature
        rel_ie2_len = float(i_len) / e2_len # 1 feature
        rel_e1e2_len = float(e1_len) / e2_len # 1 featur
        rel_log_e1i_len = my_round( math.log(rel_e1i_len,decimal), decimal ) # 1 feature
        rel_log_ie2_len = my_round( math.log(rel_ie2_len,decimal), decimal ) # 1 feature
        rel_log_e1e2_len = my_round( math.log(rel_e1e2_len,decimal), decimal ) # 1 feature
        rel_e1i_len = my_round( rel_e1i_len, decimal )
        rel_ie2_len = my_round( rel_ie2_len, decimal )
        rel_e1e2_len = my_round( rel_e1e2_len, decimal )
        tl.append( rel_e1i_len )
        tl.append( rel_ie2_len )
        tl.append( rel_e1e2_len )
        tl.append( rel_log_e1i_len )
        tl.append( rel_log_ie2_len )
        tl.append( rel_log_e1e2_len )
        # 6 feature above

        i_GC = compute_GC_content('ss5',e1_len,i_len,0,i_len,eie_fa) # 1 feature
        e1_GC = compute_GC_content('ss5',e1_len,i_len,e1_len,0,eie_fa) # 1 feature
        e2_GC = compute_GC_content('ss3',e1_len,i_len,e2_len,0,eie_fa) # 1 feature

        el = 10 if 10 < e1_len else e1_len # exon 10 bases
        il = 20 # intron 20 bases
        e1i_GC = compute_GC_content('ss5',e1_len,i_len,el,il,eie_fa) # 1 feature
        el = 10 if 10 < e2_len else e2_len # exon 10 bases
        il = 20 # intron 20 bases
        ie2_GC = compute_GC_content('ss3',e1_len,i_len,el,il,eie_fa) # 1 feature
        se = 300 if 300 < i_len else i_len # 300 bases in intron
        i_f_GC = compute_GC_content('ss5',e1_len,i_len,0,se,eie_fa) # 1 feature
        i_l_GC = compute_GC_content('ss3',e1_len,i_len,0,se,eie_fa) # 1 feature
        tl.append( i_GC )
        tl.append( e1_GC )
        tl.append( e2_GC )
        tl.append( e1i_GC )
        tl.append( ie2_GC )
        tl.append( i_f_GC )
        tl.append( i_l_GC )
        # 7 feature above

        rel_i_g_start_pos = my_round( float(eie[9]), decimal ) # 1 feature
        rel_i_g_end_pos = my_round( float(eie[10]), decimal )  # 1 feature
        gene_type = int(eie[11]) # 1 feature,the type of gene where intron located in, protein_coding:2, lincRNA:1, otherwise:0
        ss5_type = splice_site_type('ss5',e1_len,i_len,0,2,eie_fa) # 1 feature, 'GT' essential, otherwise cryptic
        ss3_type = splice_site_type('ss3',e1_len,i_len,0,2,eie_fa) # 1 feature, 'AG' essential, otherwise cryptic
        ex,in5,in3 = 3,6,20 # exon 3 bases, 9-3 = 6 bases of intron in ss5, 23-3 = 20 bases of intron in ss3
        ss5_ent = entroy('ss5',e1_len,i_len,ex,in5,eie_fa) # 1 feature
        ss3_ent = entroy('ss3',e1_len,i_len,ex,in3,eie_fa) # 1 feature
        ss5_max_ent = my_round( float(eie[-3]), decimal ) # 1 feature
        ss3_max_ent = my_round( float(eie[-2]), decimal ) # 1 feature
        if branchpoint_flag == 'off':
            pass # not get feature of branchpoint
        else:
            # branchpoint_flag == 'on'
            bp_pos = int(eie[-1]) # position of branchpoint
            if strand == '+': rel_bp_i_pos = my_round( (bp_pos - ix2) / float(i_len), decimal ) # 1 feature
            elif strand == '-': rel_bp_i_pos = my_round( abs(bp_pos - ix3) / float(i_len), decimal )
            bp_fa = bp_feature_fa[2*i+1] # sequence of branchpoint of each EIE
            bps = get_right_seq_bp(bp_fa) # 11 feature, branchpoint (bps[5]) and nearby 10 bases
            tl.append( rel_bp_i_pos )
            for f in bps:
                tl.append( f )
        tl.append( rel_i_g_start_pos )
        tl.append( rel_i_g_end_pos )
        tl.append( gene_type )
        tl.append( ss5_type )
        tl.append( ss3_type )
        tl.append( ss5_ent )
        tl.append( ss3_ent )
        tl.append( ss5_max_ent )
        tl.append( ss3_max_ent )
        # 21 feature above

        i_din = compute_dinucleotide_content('ss5',e1_len,i_len,0,i_len,eie_fa) # 16 feature
        for f in i_din:
            tl.append( f )
        e1_din = compute_dinucleotide_content('ss5',e1_len,i_len,e1_len,0,eie_fa) # 16 feature
        for f in e1_din:
            tl.append( f )
        e2_din = compute_dinucleotide_content('ss3',e1_len,i_len,e2_len,0,eie_fa) # 16 feature
        for f in e2_din:
            tl.append( f )
        se = 300 if 300 < i_len else i_len
        i_f_din = compute_dinucleotide_content('ss5',e1_len,i_len,0,se,eie_fa) # 16 feature
        for f in i_f_din:
            tl.append( f )
        i_l_din = compute_dinucleotide_content('ss3',e1_len,i_len,0,se,eie_fa) # 16 featue
        for f in i_l_din:
            tl.append( f )
        el = 10 if 10 < e1_len else e1_len # exon 10 bases
        il = 20 # intron 20 bases
        e1i_din = compute_dinucleotide_content('ss5',e1_len,i_len,el,il,eie_fa) # 16 feature
        for f in e1i_din:
            tl.append( f )
        el = 10 if 10 < e2_len else e2_len # exon 10 bases
        il = 20 # intron 20 bases
        ie2_din = compute_dinucleotide_content('ss3',e1_len,i_len,el,il,eie_fa) # 16 feature
        for f in ie2_din:
            tl.append( f )
        # 16 * 7 = 112 feature above

        feature_data.append(tl)
    EIE_tab.close()
    return feature_data

intron_feature = extract_model_feature()

def write_model_feature(feature_data):
    with open(EIE_path+'model_intron_feature.txt','w') as f_t:
        for tl in feature_data:
            f_t.write( '\t'.join( map(str,tl) ) + '\n' )

write_model_feature(intron_feature)

