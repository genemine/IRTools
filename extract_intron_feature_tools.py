import math
decimal = 2 # keep same as extract_intron_feature.py

def read_file(filename):
    read_data = []
    with open(filename,'r') as fn:
        for l in fn:
            read_data.append(l.strip())
    return read_data

def my_round(_float, _len):
    if isinstance(_float, float):
        if str(_float)[::-1].find('.') <= _len:
            return(_float)
        if str(_float)[-1] == '5':
            return(round(float(str(_float)[:-1]+'6'), _len))
        else:
            return(round(_float, _len))
    else:
        return(round(_float, _len))

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

def compute_GC_content(sign,e1_len,i_len,el,il,seq):
    # el, il is the number of got base
    f_seq = get_seq(sign,e1_len,i_len,el,il,seq)
    GC = f_seq.count('G')+f_seq.count('C')
    f_len = el+il
    GC_content = my_round( (float(GC) / f_len), decimal )
    return GC_content

def entroy(sign,e1_len,i_len,el,il,seq):
    # el and il is the number of got base
    f_seq = get_seq(sign,e1_len,i_len,el,il,seq)
    f_len = el + il
    G_cont = f_seq.count('G') / float(f_len)
    C_cont = f_seq.count('C') / float(f_len)
    A_cont = f_seq.count('A') / float(f_len)
    T_cont = f_seq.count('T') / float(f_len)
    ent = 0
    if G_cont != 0:
        ent += (-G_cont)*math.log(G_cont,2)
    if C_cont != 0:
        ent += (-C_cont)*math.log(C_cont,2)
    if A_cont != 0:
        ent += (-A_cont)*math.log(A_cont,2)
    if T_cont != 0:
        ent += (-T_cont)*math.log(T_cont,2)
    return my_round( ent, decimal )

def compute_dinucleotide_content(sign,e1_len,i_len,el,il,seq):
    # el and il is the number of got base
    f_seq = get_seq(sign,e1_len,i_len,el,il,seq)
    f_len = el + il
    AA = f_seq.count('AA')
    AC = f_seq.count('AC')
    AG = f_seq.count('AG')
    AT = f_seq.count('AT')
    CA = f_seq.count('CA')
    CC = f_seq.count('CC')
    CG = f_seq.count('CG')
    CT = f_seq.count('CT')
    GA = f_seq.count('GA')
    GC = f_seq.count('GC')
    GG = f_seq.count('GG')
    GT = f_seq.count('GT')
    TA = f_seq.count('TA')
    TC = f_seq.count('TC')
    TG = f_seq.count('TG')
    TT = f_seq.count('TT')

    all_din = AA+AC+AG+AT+CA+CC+CG+CT+GA+GC+GG+GT+TA+TC+TG+TT
    AA_c = my_round( (float(AA) / all_din) , decimal )
    AC_c = my_round( (float(AC) / all_din) , decimal )
    AG_c = my_round( (float(AG) / all_din) , decimal )
    AT_c = my_round( (float(AT) / all_din) , decimal )
    CA_c = my_round( (float(CA) / all_din) , decimal )
    CC_c = my_round( (float(CC) / all_din) , decimal )
    CG_c = my_round( (float(CG) / all_din) , decimal )
    CT_c = my_round( (float(CT) / all_din) , decimal )
    GA_c = my_round( (float(GA) / all_din) , decimal )
    GC_c = my_round( (float(GC) / all_din) , decimal )
    GG_c = my_round( (float(GG) / all_din) , decimal )
    GT_c = my_round( (float(GT) / all_din) , decimal )
    TA_c = my_round( (float(TA) / all_din) , decimal )
    TC_c = my_round( (float(TC) / all_din) , decimal )
    TG_c = my_round( (float(TG) / all_din) , decimal )
    TT_c = my_round( (float(TT) / all_din) , decimal )
    arr_din = [AA_c,AC_c,AG_c,AT_c,CA_c,CC_c,CG_c,CT_c,GA_c,GC_c,GG_c,GT_c,TA_c,TC_c,TG_c,TT_c]
    return arr_din

def splice_site_type(sign,e1_len,i_len,el,il,seq): # ss5:GT or ss3:AG essential, otherwise cryptic
    s_s = get_seq(sign,e1_len,i_len,el,il,seq) # splice site
    if sign == 'ss5':
        if s_s == 'GT': return 0 # essential
        else: return 1 # cryptic
    if sign == 'ss3':
        if s_s == 'AG': return 0 # essential
        else: return 1 # cryptic

def get_right_seq_bp(seq):
    # convert A 0, C 1, G 2, T 3
    seq = seq.replace('N','A') # 'N' is replaced to 'A'
    bps = [-1] * len(seq) # init
    for i,base in enumerate(seq):
        if base == 'A': bps[i] = 0
        elif base == 'C': bps[i] = 1
        elif base == 'G': bps[i] = 2
        elif base == 'T': bps[i] = 3
    return bps

def get_map_reads(map_file):
    map_results = read_file(map_file)
    flag = True # True is paired reads,False is single reads
    map_reads = 0
    for result in map_results:
        if result.find('properly') == -1:
            continue
        result = result.split(' ')
        if int(result[0]) > 0: # paired reads
            map_reads = result[0]
            break
        else: # single reads
            flag = False
    if flag == True: # paired reads, fpkm / 2
        return int(map_reads), 2
    else: # single reads, fpkm not change
        for result in map_results:
            if result.find('mapped') == -1:
                continue
            result = result.split(' ')
            map_reads = result[0]
            break
        return int(map_reads), 1

def calculate_coverage(features, start, end, bin_size=30): # the bin size default is 30 bases
    # ceil operation not identify between python2 and python3,so need use ' -(-x//y) '
    bin_num = int( -( -(end - start) // bin_size ) ) # equal to ceil opertion, the number of bins
    bins = [0] *  bin_num
    for f in features:
        s = int(f[0]) # start of overlap feature
        e = int(f[1]) # end of overlap feature
        s_bin = int( (s-start) // float(bin_size) ) # start bin position of coverage or overlap feature
        e_bin = int( (e-start) // float(bin_size) ) # end bin position of coverage or overlap feature
        if e_bin > bin_num -1: # this condition will occur when e_bin locate in the end of position
            e_bin = bin_num -1 # for example bin_num = 3; index 0,1,2; 90/30 = 3; acutal can be classified as 2nd bin
        for pos in range(s_bin,e_bin+1,1): # e_bin+1 because e_bin also is coverage
            bins[pos] = 1
    bin_coverage_count = 0 # the number of bins was coverage
    for bin_c in bins:
        if bin_c == 1:
            bin_coverage_count += 1
    coverage = my_round( bin_coverage_count / float(bin_num), decimal)
    return coverage

def get_read_count_and_coverage(intersect_feature_file):
    introns_count= {}
    up_exons_count = {}
    down_exons_count = {}
    introns_coverage = {}
    up_exons_coverage = {}
    down_exons_coverage = {}
    introns_feature = [] # store current intron_id
    up_exons_feature = [] # store current intron_id
    down_exons_feature = [] # store current intron_id
    intersect_data = read_file(intersect_feature_file)
    len_inter = len(intersect_data)
    for i, val in enumerate(intersect_data):
        l = val.split('\t')
        intron_id = l[3]
        strand = l[5]
        if introns_count.get(intron_id,-1) == -1: # avoid skipped when not intersect with this intron_id
            introns_count[intron_id] = 0
            introns_coverage[intron_id] = 0
            up_exons_count[intron_id] = 0
            up_exons_coverage[intron_id] = 0
            down_exons_count[intron_id] = 0
            down_exons_coverage[intron_id] = 0
        if l[6] == '.':
            continue
        eie = l[3].split('|')
        # intron
        ix1,ix2,ix3,ix4 = int(eie[1]),int(eie[2]),int(eie[3]),int(eie[4])
        start_read,end_read = int(l[7]),int(l[8])
        if (ix2 <= start_read < ix3) or (ix2 < end_read <= ix3) \
            or (start_read < ix2 and end_read > ix3):
            if start_read < ix2:
                start_read = ix2
            if end_read > ix3:
                end_read = ix3
            read_coor = [start_read, end_read]
            introns_feature.append(read_coor)
        # up '+', 'down' '-'
        start_read,end_read = int(l[7]),int(l[8])
        if (ix1 <= start_read < ix2) or (ix1 < end_read <= ix2) \
            or (start_read < ix1 and end_read > ix2):
            if start_read < ix1:
                start_read = ix1
            if end_read > ix2:
                end_read = ix2
            read_coor = [start_read, end_read]
            if strand == '+':
                up_exons_feature.append(read_coor)
            else:
                down_exons_feature.append(read_coor)
        # down '+', up '-'
        start_read,end_read = int(l[7]),int(l[8])
        if (ix3 <= start_read < ix4) or (ix3 < end_read <= ix4) \
            or (start_read < ix3 and end_read > ix4):
            if start_read < ix3:
                start_read = ix3
            if end_read > ix4:
                end_read = ix4
            read_coor = [start_read, end_read]
            if strand == '+':
                down_exons_feature.append(read_coor)
            else:
                up_exons_feature.append(read_coor)
        
        # calculated read counts and coverage
        if ((i+1) < len_inter and intron_id != intersect_data[i+1].split('\t')[3]) or ((i+1) == len_inter):
            i_r_c = len(introns_feature)
            introns_count[intron_id] = i_r_c
            introns_coverage[intron_id] = calculate_coverage(introns_feature, ix2, ix3) if i_r_c != 0 else 0
            ind1 = ix1 if strand=='+' else ix3 # swap index
            ind2 = ix2 if strand=='+' else ix4
            ind3 = ix3 if strand=='+' else ix1
            ind4 = ix4 if strand=='+' else ix2
            u_r_c = len(up_exons_feature)
            up_exons_count[intron_id] = u_r_c
            up_exons_coverage[intron_id] = calculate_coverage(up_exons_feature, ind1, ind2) if u_r_c != 0 else 0
            d_r_c = len(down_exons_feature)
            down_exons_count[intron_id] = d_r_c
            down_exons_coverage[intron_id] = calculate_coverage(down_exons_feature, ind3, ind4) if d_r_c != 0 else 0
            introns_feature = []
            up_exons_feature = []
            down_exons_feature = []

    return introns_count,up_exons_count,down_exons_count,introns_coverage,up_exons_coverage,down_exons_coverage
