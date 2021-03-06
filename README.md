# 1 IRTools
## 1.1 Description

IRTools (A feature extraction tool for intron retention analysis) is a tool to extract intron features from the genome file, annotation file and RNA-seq data file. The use of IRTools is divided into two steps: (1) the first step is to build a reference, that is, to extract intron sequence features based on the species genome and annotation. (2) the second step is to extract intron expression features using the RNA-Seq BAM file. IRTools requires three inputs: the RNA-Seq file, reference genome and annotation file. The output file stores a list of introns together with their sequence and expression features. [Feature_list.xlsx](https://github.com/genemine/IRTools/blob/master/Feature_list.xlsx) describes all the features in detail.

## 1.2 Download

- Software

IRTools is implemented in Python and can run on Linux and Mac operating systems. They are freely available for non-commercial use.

| Version        | Changes           |
| :-------------: | :----------------: |
|  [IRTools.zip](https://github.com/genemine/IRTools/archive/refs/heads/master.zip)     | Multi-thread computing implemented |



- Ensembl genome and annotation file

You can download Ensembl genome and annotation file from [the website of Ensembl](http://asia.ensembl.org/index.html "Ensembl")

# 2 Dependencies

- Bedtools (v2.25.0)
- Samtools (v1.6)
- Python module: os, sys, math, optparse v1.5.3, [ If you want to extract branchpoint sequence feature of intron, you need install tensorflow (v2.3.0), keras (v2.4.3), h5py (v2.10.0), numpy (v1.18.5) for Python3, or you need install tensorflow (v1.10.0), keras (v2.0.4), h5py (v2.10.0), numpy (v1.14.5) for Python2 ]
- Perl module : strict

# 3 Install

After downloading the source file, unzip it, and add the IRTools package path to your environmental variable PATH by modifying your .zshrc or .bashrc or .bash_profile file in your home directory.
Note 1: In the first line of the Python scripts in IRTools, the path for Python is set as follows by default:

`/usr/bin/env python`

Just in case, if your Python is not in the default path, please __change the first line (e.g. #!/usr/bin/env python) in the scripts__ so that it correctly points to Python in your machine.

Note 2: And, check that all scripts in the package are excecutable.

# 4 Usage
## 4.1 Run IRTools
Using IRTools is very simple. To facilitate testing, we have included the mouse genome (test_genome/mouse.fa), the mouse annotation file (test_genome/mouse.gtf ) and an RNA-seq data (test_data/mouse.bam) in the IRTools folder.

__Step 1: building reference to extract sequence feature__

` IRTools.py -m build -g test_genome/mouse.fa -a test_genome/mouse.gtf -i test_index_mouse `

or

` IRTools.py -m build -g test_genome/mouse.fa -a test_genome/mouse.gtf -i test_index_mouse -p off `

If want to extract branchpoint feature:

`  IRTools.py -m build -g test_genome/mouse.fa -a test_genome/mouse.gtf -i test_index_mouse -p on `

__Step 2: calculate expression feature and merge all feature__

` IRTools.py -m extract -i test_index_mouse -b test_data/mouse.bam -o test_output `

A BAM file that is generated by aligning reads to a reference genome using tools such as STAR.


Then, a file called ` feature_intron.txt ` file is generated in the ` test_output ` folder, where a list of introns together with their sequence and expression features are stored. 

## 4.2 Help document
To get more information, use:

` IRTools.py -h `

