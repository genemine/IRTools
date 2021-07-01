#!/usr/bin/env python

import os
import optparse
from extract_intron_feature_tools import get_map_reads,my_round,get_read_count_and_coverage

optParser = optparse.OptionParser()
(opts, args) = optParser.parse_args()
EIE_path = args[0]
BAM_file = args[1]
output = args[2]
threads = int(args[3])
decimal = 2 # keep 2 decimal

if int(os.popen('ls '+output+'| grep trunt_bamtobed_ | wc -l').read()) > 0:
    os.system('rm '+output+'trunt_bamtobed_*.bed') # must be removed, this is a+ or >>
###                                        ###
### prepare for intron feature of bam file ###
###                                        ###
### map reads ###
map_reads_file = output+'map_reads_file.txt'
os.system('samtools flagstat -@ '+str(threads)+' '+BAM_file+' > '+map_reads_file)
map_reads, layout = get_map_reads(map_reads_file) # all map reads, if layout is paired, fpkm / 2, othwise not change.
# for fpkm calculation
bias_intron = 6+20 # bias: calculated for intron length, 5ss:9-3=6 bases,3ss:23-3=20 bases
bias_exon = 3 # 3 bases for exon
# pre-filter bam file must overlap with EIE.bed
print('-----> Filtering bam file')
bamtobed_file = output+'bamtobed.bed'
os.system('samtools view -b -q 255 -@ '+str(threads)+' -L '+EIE_path+'EIE.bed'+' '+BAM_file+\
    ' | bedtools bamtobed -split -i stdin '+\
    " | awk '$1~/^[0-9]$/ || $1~/^[0-9][0-9]$/ || $1 ~/^[X Y]$/'"+\
    ' > '+bamtobed_file)

def extract_feature_bam(line_num):
    #######                 #######
    #######     feature     #######
    ######                  #######
    part, l_s, l_e = line_num.split('-')[:] # the part number of intron, line start, line end
    intersect_file = output+'intersect_feature_'+part+'.txt'

    # trunt the intersect file, avoid load all bamtobed file
    part_chro = os.popen("awk 'NR=="+l_s+", NR=="+l_e+"' "+EIE_path+'EIE.bed'+" | awk '{print$1}' | uniq").read().strip().split()
    part_bamtobed = output+'trunt_bamtobed_'+part+'.bed'
    # only append the read data of chrosome may be intersect
    for i_chro in part_chro:
        if i_chro not in ['X', 'Y']:
            os.system("awk '$1=="+i_chro+"{print}'"+' '+bamtobed_file+' >> '+part_bamtobed)
        else:
            os.system("awk '$1==\""+i_chro+"\"{print}'"+' '+bamtobed_file+' >> '+part_bamtobed)
    os.system("awk 'NR=="+l_s+", NR=="+l_e+"' "+EIE_path+'EIE.bed'+\
    ' | bedtools intersect -a stdin -b '+part_bamtobed+' -loj > '+intersect_file)
    introns_count,up_exons_count,down_exons_count,introns_coverage,up_exons_coverage,down_exons_coverage = get_read_count_and_coverage(intersect_file)

    feature_data = []
    EIE_bed = os.popen("awk 'NR=="+l_s+", NR=="+l_e+"' "+EIE_path+'EIE.bed')
    for val in EIE_bed: # get feature of intron
        l = val.strip().split('\t')
        eie_id = l[3]
        eie = l[3].split('|')
        ix1,ix2,ix3,ix4,strand = int(eie[1]),int(eie[2]),int(eie[3]),int(eie[4]),eie[5]
        i_len = ix3 - ix2
        if strand == '+':
            e1_len = ix2 - ix1
            e2_len = ix4 - ix3
        else:
            e1_len = ix4 - ix3
            e2_len = ix2 - ix1
        tl = [] # storing feature of bam file

        i_r_c = introns_count[eie_id] # 1 feature, intron read counts
        i_fpkm = my_round( 1000000000*i_r_c/float(i_len + bias_intron)/map_reads/layout, decimal ) # 1 feature
        i_cover = introns_coverage[eie_id] # 1 feature, intron coverage (bin size = 30 bases)
        up_ex_r_c = up_exons_count[eie_id] # 1 feature, read counts of upstream exon
        up_ex_fpkm = my_round( 1000000000*up_ex_r_c/float(e1_len + bias_exon)/map_reads/layout, decimal ) # 1 feature
        up_ex_cover = up_exons_coverage[eie_id] # 1 feature, upstream exon coverage (bin size=30 bases)
        down_ex_r_c = down_exons_count[eie_id] # 1 feature, read counts of downstream exon
        down_ex_fpkm = my_round( 1000000000*down_ex_r_c/float(e2_len + bias_exon)/map_reads/layout, decimal ) # 1 feature
        down_ex_cover = down_exons_coverage[eie_id] # 1 feature, downstream exon coverage (bin size = 30 bases)
        tl.append( i_r_c )
        tl.append( up_ex_r_c )
        tl.append( down_ex_r_c )
        tl.append( i_fpkm )
        tl.append( up_ex_fpkm )
        tl.append( down_ex_fpkm )
        tl.append( i_cover )
        tl.append( up_ex_cover )
        tl.append( down_ex_cover )
        # 9 feature above

        feature_data.append(tl)
    return feature_data

from multiprocessing import Pool
if __name__=="__main__":
    print('-----> Calculating expression feature')
    EIE_num = int(os.popen('cat '+EIE_path+'model_intron_feature.txt | wc -l').read().strip())
    single_num = EIE_num // threads
    line_start = []
    line_end = []
    for line in range(threads):
        tl = line + 1
        if tl == 1:
            line_start.append(1)
        else:
            line_start.append(((tl-1)*single_num+1))
        if tl == threads:
            line_end.append(EIE_num)
        else:
            line_end.append(tl*single_num)
    line_start_end = [str(i+1)+'-'+str(line_start[i])+'-'+str(line_end[i]) for i in range(threads)]
    pool = Pool(processes=threads)
    bam_feature = pool.map(extract_feature_bam, line_start_end)
    pool.close()
    pool.join()
    os.system('rm '+output+'trunt_bamtobed_*.bed') # must be removed, this is a+ or >>
    os.system('rm '+output+'intersect_feature_*.txt')

    feature_data_bam = []
    for ft in bam_feature:
        feature_data_bam += ft
    bp_flag = len(os.popen("awk 'NR==1' "+EIE_path+'model_intron_feature.txt').read().strip().split('\t'))
    first_line = ['intron_id','i_len','e1_len','e2_len','type_i','rel_i_rank','first_i','middle_i','last_i',
        'rel_e1i_len','rel_ie2_len','rel_e1e2_len','rel_log_e1i_len','rel_log_ie2_len','rel_log_e1e2_len',
        'i_GC','e1_GC','e2_GC','e1i_GC','ie2_GC','i_f_GC','i_l_GC','rel_bp_i_pos','up_5_branchpoint',
        'up_4_branchpoint','up_3_branchpoint','up_2_branchpoint','up_1_branchpoint', 'branchpoint',
        'down_1_branchpoint','down_2_branchpoint','down_3_branchpoint','down_4_branchpoint','down_5_branchpoint',
        'rel_i_g_start_pos','rel_i_g_end_pos','gene_type','ss5_type','ss3_type','ss5_ent','ss3_ent',
        'ss5_max_ent','ss3_max_ent','i_din_AA','i_din_AC','i_din_AG','i_din_AT','i_din_CA','i_din_CC',
        'i_din_CG','i_din_CT','i_din_GA','i_din_GC','i_din_GG','i_din_GT','i_din_TA','i_din_TC','i_din_TG',
        'i_din_TT','e1_din_AA','e1_din_AC','e1_din_AG','e1_din_AT','e1_din_CA','e1_din_CC','e1_din_CG',
        'e1_din_CT','e1_din_GA','e1_din_GC','e1_din_GG','e1_din_GT','e1_din_TA','e1_din_TC','e1_din_TG',
        'e1_din_TT','e2_din_AA','e2_din_AC','e2_din_AG','e2_din_AT','e2_din_CA','e2_din_CC','e2_din_CG',
        'e2_din_CT','e2_din_GA','e2_din_GC','e2_din_GG','e2_din_GT','e2_din_TA','e2_din_TC','e2_din_TG',
        'e2_din_TT','i_f_din_AA','i_f_din_AC','i_f_din_AG','i_f_din_AT','i_f_din_CA','i_f_din_CC','i_f_din_CG',
        'i_f_din_CT','i_f_din_GA','i_f_din_GC','i_f_din_GG','i_f_din_GT','i_f_din_TA','i_f_din_TC','i_f_din_TG',
        'i_f_din_TT','i_l_din_AA','i_l_din_AC','i_l_din_AG','i_l_din_AT','i_l_din_CA','i_l_din_CC','i_l_din_CG',
        'i_l_din_CT','i_l_din_GA','i_l_din_GC','i_l_din_GG','i_l_din_GT','i_l_din_TA','i_l_din_TC','i_l_din_TG',
        'i_l_din_TT','e1i_din_AA','e1i_din_AC','e1i_din_AG','e1i_din_AT','e1i_din_CA','e1i_din_CC','e1i_din_CG',
        'e1i_din_CT','e1i_din_GA','e1i_din_GC','e1i_din_GG','e1i_din_GT','e1i_din_TA','e1i_din_TC','e1i_din_TG',
        'e1i_din_TT','ie2_din_AA','ie2_din_AC','ie2_din_AG','ie2_din_AT','ie2_din_CA','ie2_din_CC','ie2_din_CG',
        'ie2_din_CT','ie2_din_GA','ie2_din_GC','ie2_din_GG','ie2_din_GT','ie2_din_TA','ie2_din_TC','ie2_din_TG',
        'ie2_din_TT','i_r_c','e1_r_c','e2_r_c','i_fpkm','e1_fpkm','e2_fpkm','i_cover','e1_cover','e2_cover']
    if bp_flag == 155: # intron_id + 154 feature, with branpoint
        pass
    elif bp_flag == 143: # intron_id + 142 feature, no branpoint
        first_line.remove('rel_bp_i_pos');first_line.remove('up_5_branchpoint');first_line.remove('up_4_branchpoint');first_line.remove('up_3_branchpoint')
        first_line.remove('up_2_branchpoint');first_line.remove('up_1_branchpoint');first_line.remove('branchpoint');first_line.remove('down_1_branchpoint')
        first_line.remove('down_2_branchpoint');first_line.remove('down_3_branchpoint');first_line.remove('down_4_branchpoint');first_line.remove('down_5_branchpoint')
    else: # should not be occurred
        print('-----> Index occurred an error, wrong number of features!')
    m_i_f = open(EIE_path+'model_intron_feature.txt')
    f_i = open(output+'feature_intron.txt','w')
    f_i.write('\t'.join(first_line)+'\n')
    print('-----> Merging all feature')
    for i,val in enumerate(m_i_f):
        l_l = val.strip().split('\t')
        l_r = feature_data_bam[i]
        l = l_l+l_r
        f_i.write( '\t'.join( map(str,l) ) + '\n' )
    m_i_f.close()
    f_i.close()

    print('-----> Intron features are extracted successful')
