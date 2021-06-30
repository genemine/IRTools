#!/usr/bin/python
# This is the command line entry point to extract intron feature.
# coded at: Central South University, Changsha 410083, P.R. China
# coded by: Zhenpeng Wu.
# contact: zhenpeng@csu.edu.cn
# Mar. 8, 2020

import os
import sys
import optparse

prompt_information = 'IRTools.py [-m build] [-g <genome.fa>] [-a <annotation.gtf>] [-i <index_EIE_path>] [-p <predict_branchpoint_flag>]'\
    +'\n'+'IRTools.py [-m extract] [-i <index_EIE_path>] [-b <extracted.bam>] [-o <output_path>] [-n <number_of_threads>]'
optParser = optparse.OptionParser(
    usage = prompt_information,
    )

if os.path.exists('/proc/cpuinfo'): # linux
    all_threads = int(os.popen('cat /proc/cpuinfo| grep "processor"| wc -l').read().strip())
else: # mac, not support windows
    all_threads = int(os.popen('sysctl -n machdep.cpu.thread_count').read().strip())
threads = all_threads // 2 if all_threads > 2 else 1

# mode:
# first build EIE model index, -m build;
# then extract intron feature, -m extract
optParser.add_option('-m', '--mode', action='store', type='string', dest='mode_flag', help='build: build EIE model index, extract: extract intron feature.')
optParser.add_option('-g', '--genome_file', action='store', type='string', dest='genome_file', help='genome fasta file of Ensembl.') # genome
optParser.add_option('-a', '--annotation_file', action='store', type='string', dest='annotation_file', help='an annotation file in Ensembl GTF format.') # annotation in GTF format file
optParser.add_option('-b', '--bam_file', action='store', type='string', dest='bam_file', help='extracted bam file.') # bam
optParser.add_option('-o', '--output_path', action='store', type='string', dest='output_path', help='output path.') # output
optParser.add_option('-i', '--index_EIE_path', action='store', type='string', dest='EIE_path', help='the path of EIE idnex.') # EIE index
optParser.add_option('-p', '--predicted_branchpoint_flag', action='store', type='string', dest='branchpoint_flag', help='on: predicted branchpoint, off: not predicted branchpoint.Default: off. If want to predicted branchpoint, Python need install keras, tensorflow, numpy.')
optParser.add_option('-n', '--number_of_threads', action='store', type='int', default=threads, dest='threads', help='If user not specificed, default 1/2 of all threads of machine')
(opts, args) = optParser.parse_args()

def get_abspath(rel_path,tmp_name='get_abspath_tmp.txt'):
    rel_path_file = rel_path+tmp_name
    os.system('echo "This file for produce abpath" > '+rel_path_file)
    abspath = os.path.abspath(rel_path_file).replace(tmp_name,'')
    os.system('rm '+rel_path_file)
    return abspath

def mkdir_function(path):
    if path != None and not os.path.isdir(path):
        os.makedirs(path)

### handle path of output
mkdir_function(opts.output_path)
if opts.output_path != None and opts.output_path.split('/')[-1] != '': # output path
    opts.output_path =  opts.output_path+'/'
if opts.output_path != None and os.path.isabs(opts.output_path) == False: # get abspath
    opts.output_path = get_abspath(opts.output_path)

# handle path of EIE
mkdir_function(opts.EIE_path)
if opts.EIE_path != None and opts.EIE_path.split('/')[-1] != '':
    opts.EIE_path =  opts.EIE_path+'/'
if opts.EIE_path != None and os.path.isabs(opts.EIE_path) == False: # get abspath
    opts.EIE_path = get_abspath(opts.EIE_path)

if opts.branchpoint_flag == None or opts.branchpoint_flag == 'off':
    opts.branchpoint_flag = 'off' # not get feature of branchpoint default
elif opts.branchpoint_flag != None and opts.branchpoint_flag != 'on':
    print('Please note that want to get branchpoint, you need set "-p on" next time')
    sys.exit(1)
#######      #######
####### main #######
#######      #######
if opts.mode_flag == 'build':
    prepare_EIE = 'prepare_EIE.py '+opts.genome_file+' '+opts.annotation_file+' '+opts.EIE_path+' '+opts.branchpoint_flag
    os.system(prepare_EIE)
elif opts.mode_flag == 'extract':
    extract_intron_feature = 'extract_intron_feature.py '+opts.EIE_path+' '+opts.bam_file+' '+opts.output_path+' '+str(opts.threads)
    os.system(extract_intron_feature)
else:
    print(prompt_information)
    print('Please check the model!')
