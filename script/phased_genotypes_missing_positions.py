#!/usr/local/bin/python
# Build by Chalermpong Intarat

import os, sys, requests, time, progressbar, Bio
import itertools

from Bio import SeqIO
from Bio.Alphabet import generic_dna, generic_protein
from texttable import Texttable
from termcolor import colored
from itertools import combinations

# Module: Case: Phased genotypes with missing positions
button = ""

p1a1 = ''
p2a1 = ''
p3a1 = ''
p4a1 = ''
p5a1 = ''
p6a1 = ''

p1a2 = ''
p2a2 = ''
p3a2 = ''
p4a2 = ''
p5a2 = ''
p6a2 = ''

p1a3 = ''
p2a3 = ''
p3a3 = ''
p4a3 = ''
p5a3 = ''
p6a3 = ''

p1a4 = ''
p2a4 = ''
p3a4 = ''
p4a4 = ''
p5a4 = ''
p6a4 = ''

p1a5 = ''
p2a5 = ''
p3a5 = ''
p4a5 = ''
p5a5 = ''
p6a5 = ''

p1a6 = ''
p2a6 = ''
p3a6 = ''
p4a6 = ''
p5a6 = ''
p6a6 = ''

p1a7 = ''
p2a7 = ''
p3a7 = ''
p4a7 = ''
p5a7 = ''
p6a7 = ''

p1a8 = ''
p2a8 = ''
p3a8 = ''
p4a8 = ''
p5a8 = ''
p6a8 = ''

p1a9 = ''
p2a9 = ''
p3a9 = ''
p4a9 = ''
p5a9 = ''
p6a9 = ''

p1a10 = ''
p2a10 = ''
p3a10 = ''
p4a10 = ''
p5a10 = ''
p6a10 = ''

p1a11 = ''
p2a11 = ''
p3a11 = ''
p4a11 = ''
p5a11 = ''
p6a11 = ''

p1a12 = ''
p2a12 = ''
p3a12 = ''
p4a12 = ''
p5a12 = ''
p6a12 = ''

pos_1 = '1'
pos_2 = '2'
pos_3 = '3'
pos_4 = '4'
pos_5 = '5'
pos_6 = '6'

all_1 = '*1'
all_2 = '*2'
all_3 = '*3'
all_4 = '*4'
all_5 = '*5'
all_6 = '*6'
all_7 = '*7'
all_8 = '*8'
all_9 = '*9'
all_10 = '*10'
all_11 = '*11'
all_12 = '*12'

reg_all_1 = ''
reg_all_2 = ''
reg_all_3 = ''
reg_all_4 = ''
reg_all_5 = ''
reg_all_6 = ''
reg_all_7 = ''
reg_all_8 = ''
reg_all_9 = ''
reg_all_10 = ''
reg_all_11 = ''
reg_all_12 = ''

reg_1 = ''
reg_2 = ''
reg_3 = ''
reg_4 = ''
reg_5 = ''
reg_6 = ''
reg_7 = ''
reg_8 = ''
reg_9 = ''
reg_10 = ''
reg_11 = ''
reg_12 = ''

mis_pos_1 = ''
mis_pos_2 = ''
mis_pos_3 = ''
mis_pos_4 = ''
mis_pos_5 = ''
mis_pos_6 = ''

pos_alleles_com_1_reg_1 = ''
pos_alleles_none_com_1_reg_1 = ''

pos_alleles_2_com_2_reg_1 = ''
pos_alleles_3_com_2_reg_1 = ''
pos_alleles_4_com_2_reg_1 = ''
pos_alleles_5_com_2_reg_1 = ''
pos_alleles_6_com_2_reg_1 = ''
pos_alleles_7_com_2_reg_1 = ''
pos_alleles_8_com_2_reg_1 = ''
pos_alleles_9_com_2_reg_1 = ''
pos_alleles_10_com_2_reg_1 = ''
pos_alleles_11_com_2_reg_1 = ''
pos_alleles_12_com_2_reg_1 = ''
pos_alleles_none_com_2_reg_1 = ''

pos_alleles_2_com_3_reg_1 = ''
pos_alleles_3_com_3_reg_1 = ''
pos_alleles_4_com_3_reg_1 = ''
pos_alleles_5_com_3_reg_1 = ''
pos_alleles_6_com_3_reg_1 = ''
pos_alleles_7_com_3_reg_1 = ''
pos_alleles_8_com_3_reg_1 = ''
pos_alleles_9_com_3_reg_1 = ''
pos_alleles_10_com_3_reg_1 = ''
pos_alleles_11_com_3_reg_1 = ''
pos_alleles_12_com_3_reg_1 = ''
pos_alleles_none_com_3_reg_1 = ''

pos_alleles_com_1_reg_2 = ''
pos_alleles_none_com_1_reg_2 = ''

pos_alleles_2_com_2_reg_2 = ''
pos_alleles_3_com_2_reg_2 = ''
pos_alleles_4_com_2_reg_2 = ''
pos_alleles_5_com_2_reg_2 = ''
pos_alleles_6_com_2_reg_2 = ''
pos_alleles_7_com_2_reg_2 = ''
pos_alleles_8_com_2_reg_2 = ''
pos_alleles_9_com_2_reg_2 = ''
pos_alleles_10_com_2_reg_2 = ''
pos_alleles_11_com_2_reg_2 = ''
pos_alleles_12_com_2_reg_2 = ''
pos_alleles_none_com_2_reg_2 = ''

pos_alleles_2_com_3_reg_2 = ''
pos_alleles_3_com_3_reg_2 = ''
pos_alleles_4_com_3_reg_2 = ''
pos_alleles_5_com_3_reg_2 = ''
pos_alleles_6_com_3_reg_2 = ''
pos_alleles_7_com_3_reg_2 = ''
pos_alleles_8_com_3_reg_2 = ''
pos_alleles_9_com_3_reg_2 = ''
pos_alleles_10_com_3_reg_2 = ''
pos_alleles_11_com_3_reg_2 = ''
pos_alleles_12_com_3_reg_2 = ''
pos_alleles_none_com_3_reg_2 = ''

# Create region ##################################################
def create_region():
	# Create region

    global reg_all_1, reg_1
    reg_all_1 = all_1+' '+p1a1+p2a1+p3a1+p4a1+p5a1+p6a1
    reg_1 = p1a1+p2a1+p3a1+p4a1+p5a1+p6a1
    print(reg_all_1)

    if p1a2 == '':
        reg_p1a2 = p1a1
    else:
        reg_p1a2 = colored(p1a2, 'red')
    if p2a2 == '':
       	reg_p2a2 = p2a1
    else:
        reg_p2a2 = colored(p2a2, 'red')
    if p3a2 == '':
        reg_p3a2 = p3a1
    else:
        reg_p3a2 = colored(p3a2, 'red')
    if p4a2 == '':
        reg_p4a2 = p4a1
    else:
        reg_p4a2 = colored(p4a2, 'red')
    if p5a2 == '':
        reg_p5a2 = p5a1
    else:
        reg_p5a2 = colored(p5a2, 'red')
    if p6a2 == '':
        reg_p6a2 = p6a1
    else:
        reg_p6a2 = colored(p6a2, 'red')
    global reg_all_2, reg_2
    reg_all_2 = all_2+' '+reg_p1a2+reg_p2a2+reg_p3a2+reg_p4a2+reg_p5a2+reg_p6a2
    reg_2 = reg_p1a2+reg_p2a2+reg_p3a2+reg_p4a2+reg_p5a2+reg_p6a2
    print(reg_all_2)

    if p1a3 == '':
        reg_p1a3 = p1a1
    else:
       	reg_p1a3 = colored(p1a3, 'red')
    if p2a3 == '':
        reg_p2a3 = p2a1
    else:
        reg_p2a3 = colored(p2a3, 'red')
    if p3a3 == '':
        reg_p3a3 = p3a1
    else:
        reg_p3a3 = colored(p3a3, 'red')
    if p4a3 == '':
        reg_p4a3 = p4a1
    else:
        reg_p4a3 = colored(p4a3, 'red')
    if p5a3 == '':
        reg_p5a3 = p5a1
    else:
        reg_p5a3 = colored(p5a3, 'red')
    if p6a3 == '':
        reg_p6a3 = p6a1
    else:
        reg_p6a3 = colored(p6a3, 'red')
    global reg_all_3, reg_3
    reg_all_3 = all_3+' '+reg_p1a3+reg_p2a3+reg_p3a3+reg_p4a3+reg_p5a3+reg_p6a3
    reg_3 = reg_p1a3+reg_p2a3+reg_p3a3+reg_p4a3+reg_p5a3+reg_p6a3
    print(reg_all_3)

    if p1a4 == '':
        reg_p1a4 = p1a1
    else:
        reg_p1a4 = colored(p1a4, 'red')
    if p2a4 == '':
        reg_p2a4 = p2a1
    else:
        reg_p2a4 = colored(p2a4, 'red')
    if p3a4 == '':
        reg_p3a4 = p3a1
    else:
        reg_p3a4 = colored(p3a4, 'red')
    if p4a4 == '':
        reg_p4a4 = p4a1
    else:
        reg_p4a4 = colored(p4a4, 'red')
    if p5a4 == '':
        reg_p5a4 = p5a1
    else:
        reg_p5a4 = colored(p5a4, 'red')
    if p6a4 == '':
        reg_p6a4 = p6a1
    else:
        reg_p6a4 = colored(p6a4, 'red')
    global reg_all_4, reg_4
    reg_all_4 = all_4+' '+reg_p1a4+reg_p2a4+reg_p3a4+reg_p4a4+reg_p5a4+reg_p6a4
    reg_4 = reg_p1a4+reg_p2a4+reg_p3a4+reg_p4a4+reg_p5a4+reg_p6a4
    print(reg_all_4)

    if p1a5 == '':
        reg_p1a5 = p1a1
    else:
        reg_p1a5 = colored(p1a5, 'red')
    if p2a5 == '':
        reg_p2a5 = p2a1
    else:
        reg_p2a5 = colored(p2a5, 'red')
    if p3a5 == '':
        reg_p3a5 = p3a1
    else:
        reg_p3a5 = colored(p3a5, 'red')
    if p4a5 == '':
        reg_p4a5 = p4a1
    else:
        reg_p4a5 = colored(p4a5, 'red')
    if p5a5 == '':
        reg_p5a5 = p5a1
    else:
        reg_p5a5 = colored(p5a5, 'red')
    if p6a5 == '':
        reg_p6a5 = p6a1
    else:
        reg_p6a5 = colored(p6a5, 'red')
    global reg_all_5, reg_5
    reg_all_5 = all_5+' '+reg_p1a5+reg_p2a5+reg_p3a5+reg_p4a5+reg_p5a5+reg_p6a5
    reg_5 = reg_p1a5+reg_p2a5+reg_p3a5+reg_p4a5+reg_p5a5+reg_p6a5
    print(reg_all_5)

    if p1a6 == '':
        reg_p1a6 = p1a1
    else:
        reg_p1a6 = colored(p1a6, 'red')
    if p2a6 == '':
        reg_p2a6 = p2a1
    else:
        reg_p2a6 = colored(p2a6, 'red')
    if p3a6 == '':
        reg_p3a6 = p3a1
    else:
        reg_p3a6 = colored(p3a6, 'red')
    if p4a6 == '':
        reg_p4a6 = p4a1
    else:
        reg_p4a6 = colored(p4a6, 'red')
    if p5a6 == '':
        reg_p5a6 = p5a1
    else:
        reg_p5a6 = colored(p5a6, 'red')
    if p6a6 == '':
        reg_p6a6 = p6a1
    else:
        reg_p6a6 = colored(p6a6, 'red')
    global reg_all_6, reg_6
    reg_all_6 = all_6+' '+reg_p1a6+reg_p2a6+reg_p3a6+reg_p4a6+reg_p5a6+reg_p6a6
    reg_6 = reg_p1a6+reg_p2a6+reg_p3a6+reg_p4a6+reg_p5a6+reg_p6a6
    print(reg_all_6)

    if p1a7 == '':
        reg_p1a7 = p1a1
    else:
        reg_p1a7 = colored(p1a7, 'red')
    if p2a7 == '':
        reg_p2a7 = p2a1
    else:
        reg_p2a7 = colored(p2a7, 'red')
    if p3a7 == '':
        reg_p3a7 = p3a1
    else:
        reg_p3a7 = colored(p3a7, 'red')
    if p4a7 == '':
        reg_p4a7 = p4a1
    else:
        reg_p4a7 = colored(p4a7, 'red')
    if p5a7 == '':
        reg_p5a7 = p5a1
    else:
        reg_p5a7 = colored(p5a7, 'red')
    if p6a7 == '':
        reg_p6a7 = p6a1
    else:
        reg_p6a7 = colored(p6a7, 'red')
    global reg_all_7, reg_7
    reg_all_7 = all_7+' '+reg_p1a7+reg_p2a7+reg_p3a7+reg_p4a7+reg_p5a7+reg_p6a7
    reg_7 = reg_p1a7+reg_p2a7+reg_p3a7+reg_p4a7+reg_p5a7+reg_p6a7
    print(reg_all_7)

    if p1a8 == '':
        reg_p1a8 = p1a1
    else:
        reg_p1a8 = colored(p1a8, 'red')
    if p2a8 == '':
        reg_p2a8 = p2a1
    else:
        reg_p2a8 = colored(p2a8, 'red')
    if p3a8 == '':
        reg_p3a8 = p3a1
    else:
        reg_p3a8 = colored(p3a8, 'red')
    if p4a8 == '':
        reg_p4a8 = p4a1
    else:
        reg_p4a8 = colored(p4a8, 'red')
    if p5a8 == '':
        reg_p5a8 = p5a1
    else:
        reg_p5a8 = colored(p5a8, 'red')
    if p6a8 == '':
        reg_p6a8 = p6a1
    else:
        reg_p6a8 = colored(p6a8, 'red')
    global reg_all_8, reg_8
    reg_all_8 = all_8+' '+reg_p1a8+reg_p2a8+reg_p3a8+reg_p4a8+reg_p5a8+reg_p6a8
    reg_8 = reg_p1a8+reg_p2a8+reg_p3a8+reg_p4a8+reg_p5a8+reg_p6a8
    print(reg_all_8)

    if p1a9 == '':
        reg_p1a9 = p1a1
    else:
        reg_p1a9 = colored(p1a9, 'red')
    if p2a9 == '':
        reg_p2a9 = p2a1
    else:
        reg_p2a9 = colored(p2a9, 'red')
    if p3a9 == '':
        reg_p3a9 = p3a1
    else:
        reg_p3a9 = colored(p3a9, 'red')
    if p4a9 == '':
        reg_p4a9 = p4a1
    else:
        reg_p4a9 = colored(p4a9, 'red')
    if p5a9 == '':
        reg_p5a9 = p5a1
    else:
        reg_p5a9 = colored(p5a9, 'red')
    if p6a9 == '':
        reg_p6a9 = p6a1
    else:
        reg_p6a9 = colored(p6a9, 'red')
    global reg_all_9, reg_9
    reg_all_9 = all_9+' '+reg_p1a9+reg_p2a9+reg_p3a9+reg_p4a9+reg_p5a9+reg_p6a9
    reg_9 = reg_p1a9+reg_p2a9+reg_p3a9+reg_p4a9+reg_p5a9+reg_p6a9
    print(reg_all_9)

    if p1a10 == '':
        reg_p1a10 = p1a1
    else:
        reg_p1a10 = colored(p1a10, 'red')
    if p2a10 == '':
        reg_p2a10 = p2a1
    else:
        reg_p2a10 = colored(p2a10, 'red')
    if p3a10 == '':
        reg_p3a10 = p3a1
    else:
        reg_p3a10 = colored(p3a10, 'red')
    if p4a10 == '':
        reg_p4a10 = p4a1
    else:
        reg_p4a10 = colored(p4a10, 'red')
    if p5a10 == '':
        reg_p5a10 = p5a1
    else:
        reg_p5a10 = colored(p5a10, 'red')
    if p6a10 == '':
        reg_p6a10 = p6a1
    else:
        reg_p6a10 = colored(p6a10, 'red')
    global reg_all_10, reg_10
    reg_all_10 = all_10+' '+reg_p1a10+reg_p2a10+reg_p3a10+reg_p4a10+reg_p5a10+reg_p6a10
    reg_10 = reg_p1a10+reg_p2a10+reg_p3a10+reg_p4a10+reg_p5a10+reg_p6a10
    print(reg_all_10)

    if p1a11 == '':
        reg_p1a11 = p1a1
    else:
        reg_p1a11 = colored(p1a11, 'red')
    if p2a11 == '':
        reg_p2a11 = p2a1
    else:
        reg_p2a11 = colored(p2a11, 'red')
    if p3a11 == '':
        reg_p3a11 = p3a1
    else:
        reg_p3a11 = colored(p3a11, 'red')
    if p4a11 == '':
        reg_p4a11 = p4a1
    else:
        reg_p4a11 = colored(p4a11, 'red')
    if p5a11 == '':
        reg_p5a11 = p5a1
    else:
        reg_p5a11 = colored(p5a11, 'red')
    if p6a11 == '':
        reg_p6a11 = p6a1
    else:
        reg_p6a11 = colored(p6a11, 'red')
    global reg_all_11, reg_11
    reg_all_11 = all_11+' '+reg_p1a11+reg_p2a11+reg_p3a11+reg_p4a11+reg_p5a11+reg_p6a11
    reg_11 = reg_p1a11+reg_p2a11+reg_p3a11+reg_p4a11+reg_p5a11+reg_p6a11
    print(reg_all_11)

    if p1a12 == '':
        reg_p1a12 = p1a1
    else:
        reg_p1a12 = colored(p1a12, 'red')
    if p2a12 == '':
        reg_p2a12 = p2a1
    else:
        reg_p2a12 = colored(p2a12, 'red')
    if p3a12 == '':
        reg_p3a12 = p3a1
    else:
        reg_p3a12 = colored(p3a12, 'red')
    if p4a12 == '':
        reg_p4a12 = p4a1
    else:
        reg_p4a12 = colored(p4a12, 'red')
    if p5a12 == '':
        reg_p5a12 = p5a1
    else:
        reg_p5a12 = colored(p5a12, 'red')
    if p6a12 == '':
        reg_p6a12 = p6a1
    else:
        reg_p6a12 = colored(p6a12, 'red')
    global reg_all_12, reg_12
    reg_all_12 = all_12+' '+reg_p1a12+reg_p2a12+reg_p3a12+reg_p4a12+reg_p5a12+reg_p6a12
    reg_12 = reg_p1a12+reg_p2a12+reg_p3a12+reg_p4a12+reg_p5a12+reg_p6a12
    print(reg_all_12)
    print("")

    # Create position
    if p1a2 == '':
        pos_p1a2 = p1a2
    else:
        pos_p1a2 = all_2+' '
    if p1a3 == '':
        pos_p1a3 = p1a3
    else:
        pos_p1a3 = all_3+' '
    if p1a4 == '':
        pos_p1a4 = p1a4
    else:
        pos_p1a4 = all_4+' '
    if p1a5 == '':
        pos_p1a5 = p1a5
    else:
        pos_p1a5 = all_5+' '
    if p1a6 == '':
        pos_p1a6 = p1a6
    else:
        pos_p1a6 = all_6+' '
    if p1a7 == '':
        pos_p1a7 = p1a7
    else:
        pos_p1a7 = all_7+' '
    if p1a8 == '':
        pos_p1a8 = p1a8
    else:
        pos_p1a8 = all_8+' '
    if p1a9 == '':
        pos_p1a9 = p1a9
    else:
        pos_p1a9 = all_9+' '
    if p1a10 == '':
        pos_p1a10 = p1a10
    else:
        pos_p1a10 = all_10+' '
    if p1a11 == '':
        pos_p1a11 = p1a11
    else:
        pos_p1a11 = all_11+' '
    if p1a12 == '':
        pos_p1a12 = p1a12
    else:
        pos_p1a12 = all_12+' '
    pos1 = pos_p1a2+pos_p1a3+pos_p1a4+pos_p1a5+pos_p1a6+pos_p1a7+pos_p1a8+pos_p1a9+pos_p1a10+pos_p1a11+pos_p1a12

    if p2a2 == '':
        pos_p2a2 = p2a2
    else:
        pos_p2a2 = all_2+' '
    if p2a3 == '':
        pos_p2a3 = p2a3
    else:
        pos_p2a3 = all_3+' '
    if p2a4 == '':
        pos_p2a4 = p2a4
    else:
        pos_p2a4 = all_4+' '
    if p2a5 == '':
        pos_p2a5 = p2a5
    else:
        pos_p2a5 = all_5+' '
    if p2a6 == '':
        pos_p2a6 = p2a6
    else:
        pos_p2a6 = all_6+' '
    if p2a7 == '':
        pos_p2a7 = p2a7
    else:
        pos_p2a7 = all_7+' '
    if p2a8 == '':
        pos_p2a8 = p2a8
    else:
        pos_p2a8 = all_8+' '
    if p2a9 == '':
        pos_p2a9 = p2a9
    else:
        pos_p2a9 = all_9+' '
    if p2a10 == '':
        pos_p2a10 = p2a10
    else:
        pos_p2a10 = all_10+' '
    if p2a11 == '':
        pos_p2a11 = p2a11
    else:
        pos_p2a11 = all_11+' '
    if p2a12 == '':
        pos_p2a12 = p2a12+' '
    else:
        pos_p2a12 = all_12
    pos2 = pos_p2a2+pos_p2a3+pos_p2a4+pos_p2a5+pos_p2a6+pos_p2a7+pos_p2a8+pos_p2a9+pos_p2a10+pos_p2a11+pos_p2a12

    if p3a2 == '':
        pos_p3a2 = p3a2
    else:
        pos_p3a2 = all_2+' '
    if p3a3 == '':
        pos_p3a3 = p3a3
    else:
        pos_p3a3 = all_3+' '
    if p3a4 == '':
        pos_p3a4 = p3a4
    else:
        pos_p3a4 = all_4+' '
    if p3a5 == '':
        pos_p3a5 = p3a5
    else:
        pos_p3a5 = all_5+' '
    if p3a6 == '':
        pos_p3a6 = p3a6
    else:
        pos_p3a6 = all_6+' '
    if p3a7 == '':
        pos_p3a7 = p3a7
    else:
        pos_p3a7 = all_7+' '
    if p3a8 == '':
        pos_p3a8 = p3a8
    else:
        pos_p3a8 = all_8+' '
    if p3a9 == '':
        pos_p3a9 = p3a9
    else:
        pos_p3a9 = all_9+' '
    if p3a10 == '':
        pos_p3a10 = p3a10
    else:
        pos_p3a10 = all_10+' '
    if p3a11 == '':
        pos_p3a11 = p3a11
    else:
        pos_p3a11 = all_11+' '
    if p3a12 == '':
        pos_p3a12 = p3a12
    else:
        pos_p3a12 = all_12+' '
    pos3 = pos_p3a2+pos_p3a3+pos_p3a4+pos_p3a5+pos_p3a6+pos_p3a7+pos_p3a8+pos_p3a9+pos_p3a10+pos_p3a11+pos_p3a12

    if p4a2 == '':
        pos_p4a2 = p4a2
    else:
        pos_p4a2 = all_2+' '
    if p4a3 == '':
        pos_p4a3 = p4a3
    else:
        pos_p4a3 = all_3+' '
    if p4a4 == '':
        pos_p4a4 = p4a4
    else:
        pos_p4a4 = all_4+' '
    if p4a5 == '':
        pos_p4a5 = p4a5
    else:
        pos_p4a5 = all_5+' '
    if p4a6 == '':
        pos_p4a6 = p4a6
    else:
        pos_p4a6 = all_6+' '
    if p4a7 == '':
        pos_p4a7 = p4a7
    else:
        pos_p4a7 = all_7+' '
    if p4a8 == '':
        pos_p4a8 = p4a8
    else:
        pos_p4a8 = all_8+' '
    if p4a9 == '':
        pos_p4a9 = p4a9
    else:
        pos_p4a9 = all_9+' '
    if p4a10 == '':
        pos_p4a10 = p4a10
    else:
        pos_p4a10 = all_10+' '
    if p4a11 == '':
        pos_p4a11 = p4a11
    else:
        pos_p4a11 = all_11+' '
    if p4a12 == '':
        pos_p4a12 = p4a12
    else:
        pos_p4a12 = all_12+' '
    pos4 = pos_p4a2+pos_p4a3+pos_p4a4+pos_p4a5+pos_p4a6+pos_p4a7+pos_p4a8+pos_p4a9+pos_p4a10+pos_p4a11+pos_p4a12

    if p5a2 == '':
        pos_p5a2 = p5a2
    else:
        pos_p5a2 = all_2+' '
    if p5a3 == '':
        pos_p5a3 = p5a3
    else:
        pos_p5a3 = all_3+' '
    if p5a4 == '':
        pos_p5a4 = p5a4
    else:
        pos_p5a4 = all_4+' '
    if p5a5 == '':
        pos_p5a5 = p5a5
    else:
        pos_p5a5 = all_5+' '
    if p5a6 == '':
        pos_p5a6 = p5a6
    else:
        pos_p5a6 = all_6+' '
    if p5a7 == '':
        pos_p5a7 = p5a7
    else:
        pos_p5a7 = all_7+' '
    if p5a8 == '':
        pos_p5a8 = p5a8
    else:
        pos_p5a8 = all_8+' '
    if p5a9 == '':
        pos_p5a9 = p5a9
    else:
        pos_p5a9 = all_9+' '
    if p5a10 == '':
        pos_p5a10 = p5a10
    else:
        pos_p5a10 = all_10+' '
    if p5a11 == '':
        pos_p5a11 = p5a11
    else:
        pos_p5a11 = all_11+' '
    if p5a12 == '':
        pos_p5a12 = p5a12
    else:
        pos_p5a12 = all_12+' '
    pos5 = pos_p5a2+pos_p5a3+pos_p5a4+pos_p5a5+pos_p5a6+pos_p5a7+pos_p5a8+pos_p5a9+pos_p5a10+pos_p5a11+pos_p5a12

    if p6a2 == '':
        pos_p6a2 = p6a2
    else:
        pos_p6a2 = all_2+' '
    if p6a3 == '':
        pos_p6a3 = p6a3
    else:
        pos_p6a3 = all_3+' '
    if p6a4 == '':
        pos_p6a4 = p6a4
    else:
        pos_p6a4 = all_4+' '
    if p6a5 == '':
        pos_p6a5 = p6a5
    else:
        pos_p6a5 = all_5+' '
    if p6a6 == '':
        pos_p6a6 = p6a6
    else:
        pos_p6a6 = all_6+' '
    if p6a7 == '':
        pos_p6a7 = p6a7
    else:
        pos_p6a7 = all_7+' '
    if p6a8 == '':
        pos_p6a8 = p6a8
    else:
        pos_p6a8 = all_8+' '
    if p6a9 == '':
        pos_p6a9 = p6a9
    else:
         pos_p6a9 = all_9+' '
    if p6a10 == '':
        pos_p6a10 = p6a10
    else:
        pos_p6a10 = all_10+' '
    if p6a11 == '':
        pos_p6a11 = p6a11
    else:
        pos_p6a11 = all_11+' '
    if p6a12 == '':
        pos_p6a12 = p6a12
    else:
        pos_p6a12 = all_12+' '
    pos6 = pos_p6a2+pos_p6a3+pos_p6a4+pos_p6a5+pos_p6a6+pos_p6a7+pos_p6a8+pos_p6a9+pos_p6a10+pos_p6a11+pos_p6a12

    tbl2 = Texttable()
    tbl2.add_rows([['Position', 'Allele']
        , [pos_1, pos1]
        , [pos_2, pos2]
        , [pos_3, pos3]
        , [pos_4, pos4]
        , [pos_5, pos5]
        , [pos_6, pos6]
        ])
    print(tbl2.draw())
    print("")

    # Create allele
    if p1a2 == '':
        all_p1a2 = p1a2
    else:
        all_p1a2 = pos_1+' '
    if p2a2 == '':
        all_p2a2 = p2a2
    else:
        all_p2a2 = pos_2+' '
    if p3a2 == '':
        all_p3a2 = p3a2
    else:
        all_p3a2 = pos_3+' '
    if p4a2 == '':
        all_p4a2 = p4a2
    else:
        all_p4a2 = pos_4+' '
    if p5a2 == '':
        all_p5a2 = p5a2
    else:
        all_p5a2 = pos_5+' '
    if p6a2 == '':
        all_p6a2 = p6a2
    else:
        all_p6a2 = pos_6+' '
    all2 = all_p1a2+all_p2a2+all_p3a2+all_p4a2+all_p5a2+all_p6a2

    if p1a3 == '':
        all_p1a3 = p1a3
    else:
        all_p1a3 = pos_1+' '
    if p2a3 == '':
        all_p2a3 = p2a3
    else:
        all_p2a3 = pos_2+' '
    if p3a3 == '':
        all_p3a3 = p3a3
    else:
        all_p3a3 = pos_3+' '
    if p4a3 == '':
        all_p4a3 = p4a3
    else:
        all_p4a3 = pos_4+' '
    if p5a3 == '':
        all_p5a3 = p5a3
    else:
        all_p5a3 = pos_5+' '
    if p6a3 == '':
        all_p6a3 = p6a3
    else:
        all_p6a3 = pos_6+' '
    all3 = all_p1a3+all_p2a3+all_p3a3+all_p4a3+all_p5a3+all_p6a3

    if p1a4 == '':
        all_p1a4 = p1a4
    else:
        all_p1a4 = pos_1+' '
    if p2a4 == '':
        all_p2a4 = p2a4
    else:
        all_p2a4 = pos_2+' '
    if p3a4 == '':
        all_p3a4 = p3a4
    else:
        all_p3a4 = pos_3+' '
    if p4a4 == '':
        all_p4a4 = p4a4
    else:
        all_p4a4 = pos_4+' '
    if p5a4 == '':
        all_p5a4 = p5a4
    else:
        all_p5a4 = pos_5+' '
    if p6a4 == '':
        all_p6a4 = p6a4
    else:
        all_p6a4 = pos_6+' '
    all4 = all_p1a4+all_p2a4+all_p3a4+all_p4a4+all_p5a4+all_p6a4

    if p1a5 == '':
        all_p1a5 = p1a5
    else:
        all_p1a5 = pos_1+' '
    if p2a5 == '':
        all_p2a5 = p2a5
    else:
        all_p2a5 = pos_2+' '
    if p3a5 == '':
        all_p3a5 = p3a5
    else:
        all_p3a5 = pos_3+' '
    if p4a5 == '':
        all_p4a5 = p4a5
    else:
        all_p4a5 = pos_4+' '
    if p5a5 == '':
        all_p5a5 = p5a5
    else:
        all_p5a5 = pos_5+' '
    if p6a5 == '':
        all_p6a5 = p6a5
    else:
        all_p6a5 = pos_6+' '
    all5 = all_p1a5+all_p2a5+all_p3a5+all_p4a5+all_p5a5+all_p6a5

    if p1a6 == '':
        all_p1a6 = p1a6
    else:
        all_p1a6 = pos_1+' '
    if p2a6 == '':
        all_p2a6 = p2a6
    else:
        all_p2a6 = pos_2+' '
    if p3a6 == '':
        all_p3a6 = p3a6
    else:
        all_p3a6 = pos_3+' '
    if p4a6 == '':
        all_p4a6 = p4a6
    else:
        all_p4a6 = pos_4+' '
    if p5a6 == '':
        all_p5a6 = p5a6
    else:
        all_p5a6 = pos_5+' '
    if p6a6 == '':
        all_p6a6 = p6a6
    else:
        all_p6a6 = pos_6+' '
    all6 = all_p1a6+all_p2a6+all_p3a6+all_p4a6+all_p5a6+all_p6a6

    if p1a7 == '':
        all_p1a7 = p1a7
    else:
        all_p1a7 = pos_1+' '
    if p2a7 == '':
        all_p2a7 = p2a7
    else:
        all_p2a7 = pos_2+' '
    if p3a7 == '':
        all_p3a7 = p3a7
    else:
        all_p3a7 = pos_3+' '
    if p4a7 == '':
        all_p4a7 = p4a7
    else:
        all_p4a7 = pos_4+' '
    if p5a7 == '':
        all_p5a7 = p5a7
    else:
        all_p5a7 = pos_5+' '
    if p6a7 == '':
        all_p6a7 = p6a7
    else:
        all_p6a7 = pos_6+' '
    all7 = all_p1a7+all_p2a7+all_p3a7+all_p4a7+all_p5a7+all_p6a7

    if p1a8 == '':
        all_p1a8 = p1a8
    else:
        all_p1a8 = pos_1+' '
    if p2a8 == '':
        all_p2a8 = p2a8
    else:
        all_p2a8 = pos_2+' '
    if p3a8 == '':
        all_p3a8 = p3a8
    else:
        all_p3a8 = pos_3+' '
    if p4a8 == '':
        all_p4a8 = p4a8
    else:
        all_p4a8 = pos_4+' '
    if p5a8 == '':
        all_p5a8 = p5a8
    else:
        all_p5a8 = pos_5+' '
    if p6a8 == '':
        all_p6a8 = p6a8
    else:
        all_p6a8 = pos_6+' '
    all8 = all_p1a8+all_p2a8+all_p3a8+all_p4a8+all_p5a8+all_p6a8

    if p1a9 == '':
        all_p1a9 = p1a9
    else:
        all_p1a9 = pos_1+' '
    if p2a9 == '':
        all_p2a9 = p2a9
    else:
        all_p2a9 = pos_2+' '
    if p3a9 == '':
        all_p3a9 = p3a9
    else:
        all_p3a9 = pos_3+' '
    if p4a9 == '':
        all_p4a9 = p4a9
    else:
        all_p4a9 = pos_4+' '
    if p5a9 == '':
        all_p5a9 = p5a9
    else:
        all_p5a9 = pos_5+' '
    if p6a9 == '':
        all_p6a9 = p6a9
    else:
        all_p6a9 = pos_6+' '
    all9 = all_p1a9+all_p2a9+all_p3a9+all_p4a9+all_p5a9+all_p6a9

    if p1a10 == '':
        all_p1a10 = p1a10
    else:
        all_p1a10 = pos_1+' '
    if p2a10 == '':
        all_p2a10 = p2a10
    else:
        all_p2a10 = pos_2+' '
    if p3a10 == '':
        all_p3a10 = p3a10
    else:
        all_p3a10 = pos_3+' '
    if p4a10 == '':
        all_p4a10 = p4a10
    else:
        all_p4a10 = pos_4+' '
    if p5a10 == '':
        all_p5a10 = p5a10
    else:
        all_p5a10 = pos_5+' '
    if p6a10 == '':
        all_p6a10 = p6a10
    else:
        all_p6a10 = pos_6+' '
    all10 = all_p1a10+all_p2a10+all_p3a10+all_p4a10+all_p5a10+all_p6a10

    if p1a11 == '':
        all_p1a11 = p1a11
    else:
        all_p1a11 = pos_1+' '
    if p2a11 == '':
        all_p2a11 = p2a11
    else:
        all_p2a11 = pos_2+' '
    if p3a11 == '':
        all_p3a11 = p3a11
    else:
        all_p3a11 = pos_3+' '
    if p4a11 == '':
        all_p4a11 = p4a11
    else:
        all_p4a11 = pos_4+' '
    if p5a11 == '':
        all_p5a11 = p5a11
    else:
        all_p5a11 = pos_5+' '
    if p6a11 == '':
        all_p6a11 = p6a11
    else:
        all_p6a11 = pos_6+' '
    all11 = all_p1a11+all_p2a11+all_p3a11+all_p4a11+all_p5a11+all_p6a11

    if p1a12 == '':
        all_p1a12 = p1a12
    else:
        all_p1a12 = pos_1+' '
    if p2a12 == '':
        all_p2a12 = p2a12
    else:
        all_p2a12 = pos_2+' '
    if p3a12 == '':
        all_p3a12 = p3a12
    else:
        all_p3a12 = pos_3+' '
    if p4a12 == '':
        all_p4a12 = p4a12
    else:
        all_p4a12 = pos_4+' '
    if p5a12 == '':
        all_p5a12 = p5a12
    else:
        all_p5a12 = pos_5+' '
    if p6a12 == '':
        all_p6a12 = p6a12
    else:
        all_p6a12 = pos_6+' '
    all12 = all_p1a12+all_p2a12+all_p3a12+all_p4a12+all_p5a12+all_p6a12

    tbl3 = Texttable()
    tbl3.add_rows([['Allele', 'Position']
        , [all_2, all2]
        , [all_3, all3]
        , [all_4, all4]
        , [all_5, all5]
        , [all_6, all6]
        , [all_7, all7]
        , [all_8, all8]
        , [all_9, all9]
        , [all_10, all10]
        , [all_11, all11]
        , [all_12, all12]
        ])
    print(tbl3.draw())
    print("")

# Create table ##################################################
def create_table():
    while True:

        # Module name
        print("Module 1: Create table")
        print("")

        print("Please enter Allele and Position.")
        print(">>> Enter A for A.")
        print(">>> Enter T for T.")
        print(">>> Enter C for C.")
        print(">>> Enter G for G.")
        print(">>> Enter S for Space.")
        print("")

        global p1a1
        input_p1a1 = ''
        while input_p1a1 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p1a1 = str(input("Allele 1 Position 1: "))
        	if input_p1a1 == 'A':
        		p1a1 = input_p1a1
        	elif input_p1a1 == 'T':
        		p1a1 = input_p1a1
        	elif input_p1a1 == 'C':
        		p1a1 = input_p1a1
        	elif input_p1a1 == 'G':
        		p1a1 = input_p1a1
        	elif input_p1a1 == 'S':
        		p1a1 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p2a1
        input_p2a1 = ''
        while input_p2a1 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p2a1 = str(input("Allele 1 Position 2: "))
        	if input_p2a1 == 'A':
        		p2a1 = input_p2a1
        	elif input_p2a1 == 'T':
        		p2a1 = input_p2a1
        	elif input_p2a1 == 'C':
        		p2a1 = input_p2a1
        	elif input_p2a1 == 'G':
        		p2a1 = input_p2a1
        	elif input_p2a1 == 'S':
        		p2a1 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p3a1
        input_p3a1 = ''
        while input_p3a1 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p3a1 = str(input("Allele 1 Position 3: "))
        	if input_p3a1 == 'A':
        		p3a1 = input_p3a1
        	elif input_p3a1 == 'T':
        		p3a1 = input_p3a1
        	elif input_p3a1 == 'C':
        		p3a1 = input_p3a1
        	elif input_p3a1 == 'G':
        		p3a1 = input_p3a1
        	elif input_p3a1 == 'S':
        		p3a1 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p4a1
        input_p4a1 = ''
        while input_p4a1 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p4a1 = str(input("Allele 1 Position 4: "))
        	if input_p4a1 == 'A':
        		p4a1 = input_p4a1
        	elif input_p4a1 == 'T':
        		p4a1 = input_p4a1
        	elif input_p4a1 == 'C':
        		p4a1 = input_p4a1
        	elif input_p4a1 == 'G':
        		p4a1 = input_p4a1
        	elif input_p4a1 == 'S':
        		p4a1 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p5a1
        input_p5a1 = ''
        while input_p5a1 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p5a1 = str(input("Allele 1 Position 5: "))
        	if input_p5a1 == 'A':
        		p5a1 = input_p5a1
        	elif input_p5a1 == 'T':
        		p5a1 = input_p5a1
        	elif input_p5a1 == 'C':
        		p5a1 = input_p5a1
        	elif input_p5a1 == 'G':
        		p5a1 = input_p5a1
        	elif input_p5a1 == 'S':
        		p5a1 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p6a1
        input_p6a1 = ''
        while input_p6a1 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p6a1 = str(input("Allele 1 Position 6: "))
        	if input_p6a1 == 'A':
        		p6a1 = input_p6a1
        	elif input_p6a1 == 'T':
        		p6a1 = input_p6a1
        	elif input_p6a1 == 'C':
        		p6a1 = input_p6a1
        	elif input_p6a1 == 'G':
        		p6a1 = input_p6a1
        	elif input_p6a1 == 'S':
        		p6a1 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p1a2
        input_p1a2 = ''
        while input_p1a2 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p1a2 = str(input("Allele 2 Position 1: "))
        	if input_p1a2 == 'A':
        		p1a2 = input_p1a2
        	elif input_p1a2 == 'T':
        		p1a2 = input_p1a2
        	elif input_p1a2 == 'C':
        		p1a2 = input_p1a2
        	elif input_p1a2 == 'G':
        		p1a2 = input_p1a2
        	elif input_p1a2 == 'S':
        		p1a2 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p2a2
        input_p2a2 = ''
        while input_p2a2 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p2a2 = str(input("Allele 2 Position 2: "))
        	if input_p2a2 == 'A':
        		p2a2 = input_p2a2
        	elif input_p2a2 == 'T':
        		p2a2 = input_p2a2
        	elif input_p2a2 == 'C':
        		p2a2 = input_p2a2
        	elif input_p2a2 == 'G':
        		p2a2 = input_p2a2
        	elif input_p2a2 == 'S':
        		p2a2 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p3a2
        input_p3a2 = ''
        while input_p3a2 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p3a2 = str(input("Allele 2 Position 3: "))
        	if input_p3a2 == 'A':
        		p3a2 = input_p3a2
        	elif input_p3a2 == 'T':
        		p3a2 = input_p3a2
        	elif input_p3a2 == 'C':
        		p3a2 = input_p3a2
        	elif input_p3a2 == 'G':
        		p3a2 = input_p3a2
        	elif input_p3a2 == 'S':
        		p3a2 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p4a2
        input_p4a2 = ''
        while input_p4a2 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p4a2 = str(input("Allele 2 Position 4: "))
        	if input_p4a2 == 'A':
        		p4a2 = input_p4a2
        	elif input_p4a2 == 'T':
        		p4a2 = input_p4a2
        	elif input_p4a2 == 'C':
        		p4a2 = input_p4a2
        	elif input_p4a2 == 'G':
        		p4a2 = input_p4a2
        	elif input_p4a2 == 'S':
        		p4a2 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p5a2
        input_p5a2 = ''
        while input_p5a2 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p5a2 = str(input("Allele 2 Position 5: "))
        	if input_p5a2 == 'A':
        		p5a2 = input_p5a2
        	elif input_p5a2 == 'T':
        		p5a2 = input_p5a2
        	elif input_p5a2 == 'C':
        		p5a2 = input_p5a2
        	elif input_p5a2 == 'G':
        		p5a2 = input_p5a2
        	elif input_p5a2 == 'S':
        		p5a2 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p6a2
        input_p6a2 = ''
        while input_p6a2 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p6a2 = str(input("Allele 2 Position 6: "))
        	if input_p6a2 == 'A':
        		p6a2 = input_p6a2
        	elif input_p6a2 == 'T':
        		p6a2 = input_p6a2
        	elif input_p6a2 == 'C':
        		p6a2 = input_p6a2
        	elif input_p6a2 == 'G':
        		p6a2 = input_p6a2
        	elif input_p6a2 == 'S':
        		p6a2 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p1a3
        input_p1a3 = ''
        while input_p1a3 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p1a3 = str(input("Allele 3 Position 1: "))
        	if input_p1a3 == 'A':
        		p1a3 = input_p1a3
        	elif input_p1a3 == 'T':
        		p1a3 = input_p1a3
        	elif input_p1a3 == 'C':
        		p1a3 = input_p1a3
        	elif input_p1a3 == 'G':
        		p1a3 = input_p1a3
        	elif input_p1a3 == 'S':
        		p1a3 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p2a3
        input_p2a3 = ''
        while input_p2a3 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p2a3 = str(input("Allele 3 Position 2: "))
        	if input_p2a3 == 'A':
        		p2a3 = input_p2a3
        	elif input_p2a3 == 'T':
        		p2a3 = input_p2a3
        	elif input_p2a3 == 'C':
        		p2a3 = input_p2a3
        	elif input_p2a3 == 'G':
        		p2a3 = input_p2a3
        	elif input_p2a3 == 'S':
        		p2a3 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p3a3
        input_p3a3 = ''
        while input_p3a3 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p3a3 = str(input("Allele 3 Position 3: "))
        	if input_p3a3 == 'A':
        		p3a3 = input_p3a3
        	elif input_p3a3 == 'T':
        		p3a3 = input_p3a3
        	elif input_p3a3 == 'C':
        		p3a3 = input_p3a3
        	elif input_p3a3 == 'G':
        		p3a3 = input_p3a3
        	elif input_p3a3 == 'S':
        		p3a3 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p4a3
        input_p4a3 = ''
        while input_p4a3 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p4a3 = str(input("Allele 3 Position 4: "))
        	if input_p4a3 == 'A':
        		p4a3 = input_p4a3
        	elif input_p4a3 == 'T':
        		p4a3 = input_p4a3
        	elif input_p4a3 == 'C':
        		p4a3 = input_p4a3
        	elif input_p4a3 == 'G':
        		p4a3 = input_p4a3
        	elif input_p4a3 == 'S':
        		p4a3 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p5a3
        input_p5a3 = ''
        while input_p5a3 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p5a3 = str(input("Allele 3 Position 5: "))
        	if input_p5a3 == 'A':
        		p5a3 = input_p5a3
        	elif input_p5a3 == 'T':
        		p5a3 = input_p5a3
        	elif input_p5a3 == 'C':
        		p5a3 = input_p5a3
        	elif input_p5a3 == 'G':
        		p5a3 = input_p5a3
        	elif input_p5a3 == 'S':
        		p5a3 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p6a3
        input_p6a3 = ''
        while input_p6a3 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p6a3 = str(input("Allele 3 Position 6: "))
        	if input_p6a3 == 'A':
        		p6a3 = input_p6a3
        	elif input_p6a3 == 'T':
        		p6a3 = input_p6a3
        	elif input_p6a3 == 'C':
        		p6a3 = input_p6a3
        	elif input_p6a3 == 'G':
        		p6a3 = input_p6a3
        	elif input_p6a3 == 'S':
        		p6a3 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p1a4
        input_p1a4 = ''
        while input_p1a4 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p1a4 = str(input("Allele 4 Position 1: "))
        	if input_p1a4 == 'A':
        		p1a4 = input_p1a4
        	elif input_p1a4 == 'T':
        		p1a4 = input_p1a4
        	elif input_p1a4 == 'C':
        		p1a4 = input_p1a4
        	elif input_p1a4 == 'G':
        		p1a4 = input_p1a4
        	elif input_p1a4 == 'S':
        		p1a4 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p2a4
        input_p2a4 = ''
        while input_p2a4 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p2a4 = str(input("Allele 4 Position 2: "))
        	if input_p2a4 == 'A':
        		p2a4 = input_p2a4
        	elif input_p2a4 == 'T':
        		p2a4 = input_p2a4
        	elif input_p2a4 == 'C':
        		p2a4 = input_p2a4
        	elif input_p2a4 == 'G':
        		p2a4 = input_p2a4
        	elif input_p2a4 == 'S':
        		p2a4 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p3a4
        input_p3a4 = ''
        while input_p3a4 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p3a4 = str(input("Allele 4 Position 3: "))
        	if input_p3a4 == 'A':
        		p3a4 = input_p3a4
        	elif input_p3a4 == 'T':
        		p3a4 = input_p3a4
        	elif input_p3a4 == 'C':
        		p3a4 = input_p3a4
        	elif input_p3a4 == 'G':
        		p3a4 = input_p3a4
        	elif input_p3a4 == 'S':
        		p3a4 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p4a4
        input_p4a4 = ''
        while input_p4a4 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p4a4 = str(input("Allele 4 Position 4: "))
        	if input_p4a4 == 'A':
        		p4a4 = input_p4a4
        	elif input_p4a4 == 'T':
        		p4a4 = input_p4a4
        	elif input_p4a4 == 'C':
        		p4a4 = input_p4a4
        	elif input_p4a4 == 'G':
        		p4a4 = input_p4a4
        	elif input_p4a4 == 'S':
        		p4a4 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p5a4
        input_p5a4 = ''
        while input_p5a4 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p5a4 = str(input("Allele 4 Position 5: "))
        	if input_p5a4 == 'A':
        		p5a4 = input_p5a4
        	elif input_p5a4 == 'T':
        		p5a4 = input_p5a4
        	elif input_p5a4 == 'C':
        		p5a4 = input_p5a4
        	elif input_p5a4 == 'G':
        		p5a4 = input_p5a4
        	elif input_p5a4 == 'S':
        		p5a4 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p6a4
        input_p6a4 = ''
        while input_p6a4 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p6a4 = str(input("Allele 4 Position 6: "))
        	if input_p6a4 == 'A':
        		p6a4 = input_p6a4
        	elif input_p6a4 == 'T':
        		p6a4 = input_p6a4
        	elif input_p6a4 == 'C':
        		p6a4 = input_p6a4
        	elif input_p6a4 == 'G':
        		p6a4 = input_p6a4
        	elif input_p6a4 == 'S':
        		p6a4 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p1a5
        input_p1a5 = ''
        while input_p1a5 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p1a5 = str(input("Allele 5 Position 1: "))
        	if input_p1a5 == 'A':
        		p1a5 = input_p1a5
        	elif input_p1a5 == 'T':
        		p1a5 = input_p1a5
        	elif input_p1a5 == 'C':
        		p1a5 = input_p1a5
        	elif input_p1a5 == 'G':
        		p1a5 = input_p1a5
        	elif input_p1a5 == 'S':
        		p1a5 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p2a5
        input_p2a5 = ''
        while input_p2a5 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p2a5 = str(input("Allele 5 Position 2: "))
        	if input_p2a5 == 'A':
        		p2a5 = input_p2a5
        	elif input_p2a5 == 'T':
        		p2a5 = input_p2a5
        	elif input_p2a5 == 'C':
        		p2a5 = input_p2a5
        	elif input_p2a5 == 'G':
        		p2a5 = input_p2a5
        	elif input_p2a5 == 'S':
        		p2a5 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p3a5
        input_p3a5 = ''
        while input_p3a5 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p3a5 = str(input("Allele 5 Position 3: "))
        	if input_p3a5 == 'A':
        		p3a5 = input_p3a5
        	elif input_p3a5 == 'T':
        		p3a5 = input_p3a5
        	elif input_p3a5 == 'C':
        		p3a5 = input_p3a5
        	elif input_p3a5 == 'G':
        		p3a5 = input_p3a5
        	elif input_p3a5 == 'S':
        		p3a5 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p4a5
        input_p4a5 = ''
        while input_p4a5 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p4a5 = str(input("Allele 5 Position 4: "))
        	if input_p4a5 == 'A':
        		p4a5 = input_p4a5
        	elif input_p4a5 == 'T':
        		p4a5 = input_p4a5
        	elif input_p4a5 == 'C':
        		p4a5 = input_p4a5
        	elif input_p4a5 == 'G':
        		p4a5 = input_p4a5
        	elif input_p4a5 == 'S':
        		p4a5 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p5a5
        input_p5a5 = ''
        while input_p5a5 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p5a5 = str(input("Allele 5 Position 5: "))
        	if input_p5a5 == 'A':
        		p5a5 = input_p5a5
        	elif input_p5a5 == 'T':
        		p5a5 = input_p5a5
        	elif input_p5a5 == 'C':
        		p5a5 = input_p5a5
        	elif input_p5a5 == 'G':
        		p5a5 = input_p5a5
        	elif input_p5a5 == 'S':
        		p5a5 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p6a5
        input_p6a5 = ''
        while input_p6a5 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p6a5 = str(input("Allele 5 Position 6: "))
        	if input_p6a5 == 'A':
        		p6a5 = input_p6a5
        	elif input_p6a5 == 'T':
        		p6a5 = input_p6a5
        	elif input_p6a5 == 'C':
        		p6a5 = input_p6a5
        	elif input_p6a5 == 'G':
        		p6a5 = input_p6a5
        	elif input_p6a5 == 'S':
        		p6a5 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p1a6
        input_p1a6 = ''
        while input_p1a6 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p1a6 = str(input("Allele 6 Position 1: "))
        	if input_p1a6 == 'A':
        		p1a6 = input_p1a6
        	elif input_p1a6 == 'T':
        		p1a6 = input_p1a6
        	elif input_p1a6 == 'C':
        		p1a6 = input_p1a6
        	elif input_p1a6 == 'G':
        		p1a6 = input_p1a6
        	elif input_p1a6 == 'S':
        		p1a6 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p2a6
        input_p2a6 = ''
        while input_p2a6 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p2a6 = str(input("Allele 6 Position 2: "))
        	if input_p2a6 == 'A':
        		p2a6 = input_p2a6
        	elif input_p2a6 == 'T':
        		p2a6 = input_p2a6
        	elif input_p2a6 == 'C':
        		p2a6 = input_p2a6
        	elif input_p2a6 == 'G':
        		p2a6 = input_p2a6
        	elif input_p2a6 == 'S':
        		p2a6 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

       	global p3a6
        input_p3a6 = ''
        while input_p3a6 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p3a6 = str(input("Allele 6 Position 3: "))
        	if input_p3a6 == 'A':
        		p3a6 = input_p3a6
        	elif input_p3a6 == 'T':
        		p3a6 = input_p3a6
        	elif input_p3a6 == 'C':
        		p3a6 = input_p3a6
        	elif input_p3a6 == 'G':
        		p3a6 = input_p3a6
        	elif input_p3a6 == 'S':
        		p3a6 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p4a6
        input_p4a6 = ''
        while input_p4a6 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p4a6 = str(input("Allele 6 Position 4: "))
        	if input_p4a6 == 'A':
        		p4a6 = input_p4a6
        	elif input_p4a6 == 'T':
        		p4a6 = input_p4a6
        	elif input_p4a6 == 'C':
        		p4a6 = input_p4a6
        	elif input_p4a6 == 'G':
        		p4a6 = input_p4a6
        	elif input_p4a6 == 'S':
        		p4a6 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p5a6
        input_p5a6 = ''
        while input_p5a6 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p5a6 = str(input("Allele 6 Position 5: "))
        	if input_p5a6 == 'A':
        		p5a6 = input_p5a6
        	elif input_p5a6 == 'T':
        		p5a6 = input_p5a6
        	elif input_p5a6 == 'C':
        		p5a6 = input_p5a6
        	elif input_p5a6 == 'G':
        		p5a6 = input_p5a6
        	elif input_p5a6 == 'S':
        		p5a6 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p6a6
        input_p6a6 = ''
        while input_p6a6 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p6a6 = str(input("Allele 6 Position 6: "))
        	if input_p6a6 == 'A':
        		p6a6 = input_p6a6
        	elif input_p6a6 == 'T':
        		p6a6 = input_p6a6
        	elif input_p6a6 == 'C':
        		p6a6 = input_p6a6
        	elif input_p6a6 == 'G':
        		p6a6 = input_p6a6
        	elif input_p6a6 == 'S':
        		p6a6 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p1a7
        input_p1a7 = ''
        while input_p1a7 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p1a7 = str(input("Allele 7 Position 1: "))
        	if input_p1a7 == 'A':
        		p1a7 = input_p1a7
        	elif input_p1a7 == 'T':
        		p1a7 = input_p1a7
        	elif input_p1a7 == 'C':
        		p1a7 = input_p1a7
        	elif input_p1a7 == 'G':
        		p1a7 = input_p1a7
        	elif input_p1a7 == 'S':
        		p1a7 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p2a7
        input_p2a7 = ''
        while input_p2a7 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p2a7 = str(input("Allele 7 Position 2: "))
        	if input_p2a7 == 'A':
        		p2a7 = input_p2a7
        	elif input_p2a7 == 'T':
        		p2a7 = input_p2a7
        	elif input_p2a7 == 'C':
        		p2a7 = input_p2a7
        	elif input_p2a7 == 'G':
        		p2a7 = input_p2a7
        	elif input_p2a7 == 'S':
        		p2a7 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p3a7
        input_p3a7 = ''
        while input_p3a7 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p3a7 = str(input("Allele 7 Position 3: "))
        	if input_p3a7 == 'A':
        		p3a7 = input_p3a7
        	elif input_p3a7 == 'T':
        		p3a7 = input_p3a7
        	elif input_p3a7 == 'C':
        		p3a7 = input_p3a7
        	elif input_p3a7 == 'G':
        		p3a7 = input_p3a7
        	elif input_p3a7 == 'S':
        		p3a7 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p4a7
        input_p4a7 = ''
        while input_p4a7 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p4a7 = str(input("Allele 7 Position 4: "))
        	if input_p4a7 == 'A':
        		p4a7 = input_p4a7
        	elif input_p4a7 == 'T':
        		p4a7 = input_p4a7
        	elif input_p4a7 == 'C':
        		p4a7 = input_p4a7
        	elif input_p4a7 == 'G':
        		p4a7 = input_p4a7
        	elif input_p4a7 == 'S':
        		p4a7 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p5a7
        input_p5a7 = ''
        while input_p5a7 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p5a7 = str(input("Allele 7 Position 5: "))
        	if input_p5a7 == 'A':
        		p5a7 = input_p5a7
        	elif input_p5a7 == 'T':
        		p5a7 = input_p5a7
        	elif input_p5a7 == 'C':
        		p5a7 = input_p5a7
        	elif input_p5a7 == 'G':
        		p5a7 = input_p5a7
        	elif input_p5a7 == 'S':
        		p5a7 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p6a7
        input_p6a7 = ''
        while input_p6a7 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p6a7 = str(input("Allele 7 Position 6: "))
        	if input_p6a7 == 'A':
        		p6a7 = input_p6a7
        	elif input_p6a7 == 'T':
        		p6a7 = input_p6a7
        	elif input_p6a7 == 'C':
        		p6a7 = input_p6a7
        	elif input_p6a7 == 'G':
        		p6a7 = input_p6a7
        	elif input_p6a7 == 'S':
        		p6a7 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p1a8
        input_p1a8 = ''
        while input_p1a8 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p1a8 = str(input("Allele 8 Position 1: "))
        	if input_p1a8 == 'A':
        		p1a8 = input_p1a8
        	elif input_p1a8 == 'T':
        		p1a8 = input_p1a8
        	elif input_p1a8 == 'C':
        		p1a8 = input_p1a8
        	elif input_p1a8 == 'G':
        		p1a8 = input_p1a8
        	elif input_p1a8 == 'S':
        		p1a8 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p2a8
        input_p2a8 = ''
        while input_p2a8 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p2a8 = str(input("Allele 8 Position 2: "))
        	if input_p2a8 == 'A':
        		p2a8 = input_p2a8
        	elif input_p2a8 == 'T':
        		p2a8 = input_p2a8
        	elif input_p2a8 == 'C':
        		p2a8 = input_p2a8
        	elif input_p2a8 == 'G':
        		p2a8 = input_p2a8
        	elif input_p2a8 == 'S':
        		p2a8 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p3a8
        input_p3a8 = ''
        while input_p3a8 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p3a8 = str(input("Allele 8 Position 3: "))
        	if input_p3a8 == 'A':
        		p3a8 = input_p3a8
        	elif input_p3a8 == 'T':
        		p3a8 = input_p3a8
        	elif input_p3a8 == 'C':
        		p3a8 = input_p3a8
        	elif input_p3a8 == 'G':
        		p3a8 = input_p3a8
        	elif input_p3a8 == 'S':
        		p3a8 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p4a8
        input_p4a8 = ''
        while input_p4a8 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p4a8 = str(input("Allele 8 Position 4: "))
        	if input_p4a8 == 'A':
        		p4a8 = input_p4a8
        	elif input_p4a8 == 'T':
        		p4a8 = input_p4a8
        	elif input_p4a8 == 'C':
        		p4a8 = input_p4a8
        	elif input_p4a8 == 'G':
        		p4a8 = input_p4a8
        	elif input_p4a8 == 'S':
        		p4a8 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p5a8
        input_p5a8 = ''
        while input_p5a8 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p5a8 = str(input("Allele 8 Position 5: "))
        	if input_p5a8 == 'A':
        		p5a8 = input_p5a8
        	elif input_p5a8 == 'T':
        		p5a8 = input_p5a8
        	elif input_p5a8 == 'C':
        		p5a8 = input_p5a8
        	elif input_p5a8 == 'G':
        		p5a8 = input_p5a8
        	elif input_p5a8 == 'S':
        		p5a8 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p6a8
        input_p6a8 = ''
        while input_p6a8 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p6a8 = str(input("Allele 8 Position 6: "))
        	if input_p6a8 == 'A':
        		p6a8 = input_p6a8
        	elif input_p6a8 == 'T':
        		p6a8 = input_p6a8
        	elif input_p6a8 == 'C':
        		p6a8 = input_p6a8
        	elif input_p6a8 == 'G':
        		p6a8 = input_p6a8
        	elif input_p6a8 == 'S':
        		p6a8 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p1a9
        input_p1a9 = ''
        while input_p1a9 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p1a9 = str(input("Allele 9 Position 1: "))
        	if input_p1a9 == 'A':
        		p1a9 = input_p1a9
        	elif input_p1a9 == 'T':
        		p1a9 = input_p1a9
        	elif input_p1a9 == 'C':
        		p1a9 = input_p1a9
        	elif input_p1a9 == 'G':
        		p1a9 = input_p1a9
        	elif input_p1a9 == 'S':
        		p1a9 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p2a9
        input_p2a9 = ''
        while input_p2a9 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p2a9 = str(input("Allele 9 Position 2: "))
        	if input_p2a9 == 'A':
        		p2a9 = input_p2a9
        	elif input_p2a9 == 'T':
        		p2a9 = input_p2a9
        	elif input_p2a9 == 'C':
        		p2a9 = input_p2a9
        	elif input_p2a9 == 'G':
        		p2a9 = input_p2a9
        	elif input_p2a9 == 'S':
        		p2a9 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p3a9
        input_p3a9 = ''
        while input_p3a9 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p3a9 = str(input("Allele 9 Position 3: "))
        	if input_p3a9 == 'A':
        		p3a9 = input_p3a9
        	elif input_p3a9 == 'T':
        		p3a9 = input_p3a9
        	elif input_p3a9 == 'C':
        		p3a9 = input_p3a9
        	elif input_p3a9 == 'G':
        		p3a9 = input_p3a9
        	elif input_p3a9 == 'S':
        		p3a9 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p4a9
        input_p4a9 = ''
        while input_p4a9 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p4a9 = str(input("Allele 9 Position 4: "))
        	if input_p4a9 == 'A':
        		p4a9 = input_p4a9
        	elif input_p4a9 == 'T':
        		p4a9 = input_p4a9
        	elif input_p4a9 == 'C':
        		p4a9 = input_p4a9
        	elif input_p4a9 == 'G':
        		p4a9 = input_p4a9
        	elif input_p4a9 == 'S':
        		p4a9 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p5a9
        input_p5a9 = ''
        while input_p5a9 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p5a9 = str(input("Allele 9 Position 5: "))
        	if input_p5a9 == 'A':
        		p5a9 = input_p5a9
        	elif input_p5a9 == 'T':
        		p5a9 = input_p5a9
        	elif input_p5a9 == 'C':
        		p5a9 = input_p5a9
        	elif input_p5a9 == 'G':
        		p5a9 = input_p5a9
        	elif input_p5a9 == 'S':
        		p5a9 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p6a9
        input_p6a9 = ''
        while input_p6a9 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p6a9 = str(input("Allele 9 Position 6: "))
        	if input_p6a9 == 'A':
        		p6a9 = input_p6a9
        	elif input_p6a9 == 'T':
        		p6a9 = input_p6a9
        	elif input_p6a9 == 'C':
        		p6a9 = input_p6a9
        	elif input_p6a9 == 'G':
        		p6a9 = input_p6a9
        	elif input_p6a9 == 'S':
        		p6a9 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p1a10
        input_p1a10 = ''
        while input_p1a10 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p1a10 = str(input("Allele 10 Position 1: "))
        	if input_p1a10 == 'A':
        		p1a10 = input_p1a10
        	elif input_p1a10 == 'T':
        		p1a10 = input_p1a10
        	elif input_p1a10 == 'C':
        		p1a10 = input_p1a10
        	elif input_p1a10 == 'G':
        		p1a10 = input_p1a10
        	elif input_p1a10 == 'S':
        		p1a10 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p2a10
        input_p2a10 = ''
        while input_p2a10 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p2a10 = str(input("Allele 10 Position 2: "))
        	if input_p2a10 == 'A':
        		p2a10 = input_p2a10
        	elif input_p2a10 == 'T':
        		p2a10 = input_p2a10
        	elif input_p2a10 == 'C':
        		p2a10 = input_p2a10
        	elif input_p2a10 == 'G':
        		p2a10 = input_p2a10
        	elif input_p2a10 == 'S':
        		p2a10 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p3a10
        input_p3a10 = ''
        while input_p3a10 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p3a10 = str(input("Allele 10 Position 3: "))
        	if input_p3a10 == 'A':
        		p3a10 = input_p3a10
        	elif input_p3a10 == 'T':
        		p3a10 = input_p3a10
        	elif input_p3a10 == 'C':
        		p3a10 = input_p3a10
        	elif input_p3a10 == 'G':
        		p3a10 = input_p3a10
        	elif input_p3a10 == 'S':
        		p3a10 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p4a10
        input_p4a10 = ''
        while input_p4a10 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p4a10 = str(input("Allele 10 Position 4: "))
        	if input_p4a10 == 'A':
        		p4a10 = input_p4a10
        	elif input_p4a10 == 'T':
        		p4a10 = input_p4a10
        	elif input_p4a10 == 'C':
        		p4a10 = input_p4a10
        	elif input_p4a10 == 'G':
        		p4a10 = input_p4a10
        	elif input_p4a10 == 'S':
        		p4a10 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p5a10
        input_p5a10 = ''
        while input_p5a10 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p5a10 = str(input("Allele 10 Position 5: "))
        	if input_p5a10 == 'A':
        		p5a10 = input_p5a10
        	elif input_p5a10 == 'T':
        		p5a10 = input_p5a10
        	elif input_p5a10 == 'C':
        		p5a10 = input_p5a10
        	elif input_p5a10 == 'G':
        		p5a10 = input_p5a10
        	elif input_p5a10 == 'S':
        		p5a10 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p6a10
        input_p6a10 = ''
        while input_p6a10 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p6a10 = str(input("Allele 10 Position 6: "))
        	if input_p6a10 == 'A':
        		p6a10 = input_p6a10
        	elif input_p6a10 == 'T':
        		p6a10 = input_p6a10
        	elif input_p6a10 == 'C':
        		p6a10 = input_p6a10
        	elif input_p6a10 == 'G':
        		p6a10 = input_p6a10
        	elif input_p6a10 == 'S':
        		p6a10 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p1a11
        input_p1a11 = ''
        while input_p1a11 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p1a11 = str(input("Allele 11 Position 1: "))
        	if input_p1a11 == 'A':
        		p1a11 = input_p1a11
        	elif input_p1a11 == 'T':
        		p1a11 = input_p1a11
        	elif input_p1a11 == 'C':
        		p1a11 = input_p1a11
        	elif input_p1a11 == 'G':
        		p1a11 = input_p1a11
        	elif input_p1a11 == 'S':
        		p1a11 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p2a11
        input_p2a11 = ''
        while input_p2a11 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p2a11 = str(input("Allele 11 Position 2: "))
        	if input_p2a11 == 'A':
        		p2a11 = input_p2a11
        	elif input_p2a11 == 'T':
        		p2a11 = input_p2a11
        	elif input_p2a11 == 'C':
        		p2a11 = input_p2a11
        	elif input_p2a11 == 'G':
        		p2a11 = input_p2a11
        	elif input_p2a11 == 'S':
        		p2a11 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p3a11
        input_p3a11 = ''
        while input_p3a11 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p3a11 = str(input("Allele 11 Position 3: "))
        	if input_p3a11 == 'A':
        		p3a11 = input_p3a11
        	elif input_p3a11 == 'T':
        		p3a11 = input_p3a11
        	elif input_p3a11 == 'C':
        		p3a11 = input_p3a11
        	elif input_p3a11 == 'G':
        		p3a11 = input_p3a11
        	elif input_p3a11 == 'S':
        		p3a11 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p4a11
        input_p4a11 = ''
        while input_p4a11 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p4a11 = str(input("Allele 11 Position 4: "))
        	if input_p4a11 == 'A':
        		p4a11 = input_p4a11
        	elif input_p4a11 == 'T':
        		p4a11 = input_p4a11
        	elif input_p4a11 == 'C':
        		p4a11 = input_p4a11
        	elif input_p4a11 == 'G':
        		p4a11 = input_p4a11
        	elif input_p4a11 == 'S':
        		p4a11 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p5a11
        input_p5a11 = ''
        while input_p5a11 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p5a11 = str(input("Allele 11 Position 5: "))
        	if input_p5a11 == 'A':
        		p5a11 = input_p5a11
        	elif input_p5a11 == 'T':
        		p5a11 = input_p5a11
        	elif input_p5a11 == 'C':
        		p5a11 = input_p5a11
        	elif input_p5a11 == 'G':
        		p5a11 = input_p5a11
        	elif input_p5a11 == 'S':
        		p5a11 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p6a11
        input_p6a11 = ''
        while input_p6a11 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p6a11 = str(input("Allele 11 Position 6: "))
        	if input_p6a11 == 'A':
        		p6a11 = input_p6a11
        	elif input_p6a11 == 'T':
        		p6a11 = input_p6a11
        	elif input_p6a11 == 'C':
        		p6a11 = input_p6a11
        	elif input_p6a11 == 'G':
        		p6a11 = input_p6a11
        	elif input_p6a11 == 'S':
        		p6a11 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p1a12
        input_p1a12 = ''
        while input_p1a12 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p1a12 = str(input("Allele 12 Position 1: "))
        	if input_p1a12 == 'A':
        		p1a12 = input_p1a12
        	elif input_p1a12 == 'T':
        		p1a12 = input_p1a12
        	elif input_p1a12 == 'C':
        		p1a12 = input_p1a12
        	elif input_p1a12 == 'G':
        		p1a12 = input_p1a12
        	elif input_p1a12 == 'S':
        		p1a12 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p2a12
        input_p2a12 = ''
        while input_p2a12 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p2a12 = str(input("Allele 12 Position 2: "))
        	if input_p2a12 == 'A':
        		p2a12 = input_p2a12
        	elif input_p2a12 == 'T':
        		p2a12 = input_p2a12
        	elif input_p2a12 == 'C':
        		p2a12 = input_p2a12
        	elif input_p2a12 == 'G':
        		p2a12 = input_p2a12
        	elif input_p2a12 == 'S':
        		p2a12 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p3a12
        input_p3a12 = ''
        while input_p3a12 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p3a12 = str(input("Allele 12 Position 3: "))
        	if input_p3a12 == 'A':
        		p3a12 = input_p3a12
        	elif input_p3a12 == 'T':
        		p3a12 = input_p3a12
        	elif input_p3a12 == 'C':
        		p3a12 = input_p3a12
        	elif input_p3a12 == 'G':
        		p3a12 = input_p3a12
        	elif input_p3a12 == 'S':
        		p3a12 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p4a12
        input_p4a12 = ''
        while input_p4a12 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p4a12 = str(input("Allele 12 Position 4: "))
        	if input_p4a12 == 'A':
        		p4a12 = input_p4a12
        	elif input_p4a12 == 'T':
        		p4a12 = input_p4a12
        	elif input_p4a12 == 'C':
        		p4a12 = input_p4a12
        	elif input_p4a12 == 'G':
        		p4a12 = input_p4a12
        	elif input_p4a12 == 'S':
        		p4a12 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p5a12
        input_p5a12 = ''
        while input_p5a12 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p5a12 = str(input("Allele 12 Position 5: "))
        	if input_p5a12 == 'A':
        		p5a12 = input_p5a12
        	elif input_p5a12 == 'T':
        		p5a12 = input_p5a12
        	elif input_p5a12 == 'C':
        		p5a12 = input_p5a12
        	elif input_p5a12 == 'G':
        		p5a12 = input_p5a12
        	elif input_p5a12 == 'S':
        		p5a12 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        global p6a12
        input_p6a12 = ''
        while input_p6a12 not in {'A', 'T', 'C', 'G', 'S'}:
        	input_p6a12 = str(input("Allele 12 Position 6: "))
        	if input_p6a12 == 'A':
        		p6a12 = input_p6a12
        	elif input_p6a12 == 'T':
        		p6a12 = input_p6a12
        	elif input_p6a12 == 'C':
        		p6a12 = input_p6a12
        	elif input_p6a12 == 'G':
        		p6a12 = input_p6a12
        	elif input_p6a12 == 'S':
        		p6a12 = ''
        	else:
        		print("")
        		print("Not correct!")
        		print("")

        # Create table
        ans = ''
        while ans not in {'y', 'n'}:
        	print("")
        	print("Create table ...?")
        	ans = str(input("y/n: "))
        	if ans == 'y':
        		print("")
	        	tbl1 = Texttable()
	        	tbl1.add_rows([[' ', pos_1, pos_2, pos_3, pos_4, pos_5, pos_16]
	        		, [all_1, p1a1, p2a1, p3a1, p4a1, p5a1, p6a1]
	        		, [all_2, p1a2, p2a2, p3a2, p4a2, p5a2, p6a2]
	        		, [all_3, p1a3, p2a3, p3a3, p4a3, p5a3, p6a3]
	        		, [all_4, p1a4, p2a4, p3a4, p4a4, p5a4, p6a4]
	        		, [all_5, p1a5, p2a5, p3a5, p4a5, p5a5, p6a5]
	        		, [all_6, p1a6, p2a6, p3a6, p4a6, p5a6, p6a6]
	        		, [all_7, p1a7, p2a7, p3a7, p4a7, p5a7, p6a7]
	        		, [all_8, p1a8, p2a8, p3a8, p4a8, p5a8, p6a8]
	        		, [all_9, p1a9, p2a9, p3a9, p4a9, p5a9, p6a9]
	        		, [all_10, p1a10, p2a10, p3a10, p4a10, p5a10, p6a10]
	        		, [all_11, p1a11, p2a11, p3a11, p4a11, p5a11, p6a11]
	        		, [all_12, p1a12, p2a12, p3a12, p4a12, p5a12, p6a12]
	        		])
	        	print(tbl1.draw())
        	elif ans == 'n':
        		break
        		print("")
        	else:
        		print("Not correct!")

    	# Create region
        create_region()

        # Module complete
        print("Module 1: Create table complete")
        print("")
        break

# Load table ##################################################
def load_table():
    while True:

        # Module name
        print("Module 2: Load table")
        print("")

        def_p1a1 = 'T'
        def_p2a1 = 'A'
        def_p3a1 = 'G'
        def_p4a1 = 'C'
        def_p5a1 = 'T'
        def_p6a1 = 'G'

        def_p1a2 = ''
        def_p2a2 = 'G'
        def_p3a2 = ''
        def_p4a2 = ''
        def_p5a2 = ''
        def_p6a2 = ''

        def_p1a3 = ''
        def_p2a3 = ''
        def_p3a3 = 'C'
        def_p4a3 = ''
        def_p5a3 = ''
        def_p6a3 = ''

        def_p1a4 = ''
        def_p2a4 = ''
        def_p3a4 = ''
        def_p4a4 = 'T'
        def_p5a4 = ''
        def_p6a4 = ''

        def_p1a5 = ''
        def_p2a5 = ''
        def_p3a5 = ''
        def_p4a5 = ''
        def_p5a5 = 'A'
        def_p6a5 = ''

        def_p1a6 = ''
        def_p2a6 = ''
        def_p3a6 = ''
        def_p4a6 = ''
        def_p5a6 = ''
        def_p6a6 = 'C'

        def_p1a7 = 'A'
        def_p2a7 = ''
        def_p3a7 = 'C'
        def_p4a7 = ''
        def_p5a7 = ''
        def_p6a7 = ''

        def_p1a8 = ''
        def_p2a8 = 'G'
        def_p3a8 = ''
        def_p4a8 = ''
        def_p5a8 = ''
        def_p6a8 = 'C'

        def_p1a9 = ''
        def_p2a9 = ''
        def_p3a9 = 'C'
        def_p4a9 = ''
        def_p5a9 = 'A'
        def_p6a9 = ''

        def_p1a10 = ''
        def_p2a10 = 'G'
        def_p3a10 = ''
        def_p4a10 = 'T'
        def_p5a10 = ''
        def_p6a10 = ''

        def_p1a11 = 'A'
        def_p2a11 = ''
        def_p3a11 = ''
        def_p4a11 = 'T'
        def_p5a11 = ''
        def_p6a11 = ''

        def_p1a12 = ''
        def_p2a12 = ''
        def_p3a12 = ''
        def_p4a12 = ''
        def_p5a12 = 'A'
        def_p6a12 = 'C'

        # Load table
        ans = ''
        while ans not in {'y', 'n'}:
        	print("\tLoad table ...?")
        	ans = str(input("\ty/n: "))
        	if ans == 'y':

        		global p1a1, p2a1, p3a1, p4a1, p5a1, p6a1
        		p1a1 = def_p1a1
        		p2a1 = def_p2a1
        		p3a1 = def_p3a1
        		p4a1 = def_p4a1
        		p5a1 = def_p5a1
        		p6a1 = def_p6a1

        		global p1a2, p2a2, p3a2, p4a2, p5a2, p6a2
        		p1a2 = def_p1a2
        		p2a2 = def_p2a2
        		p3a2 = def_p3a2
        		p4a2 = def_p4a2
        		p5a2 = def_p5a2
        		p6a2 = def_p6a2

        		global p1a3, p2a3, p3a3, p4a3, p5a3, p6a3
        		p1a3 = def_p1a3
        		p2a3 = def_p2a3
        		p3a3 = def_p3a3
        		p4a3 = def_p4a3
        		p5a3 = def_p5a3
        		p6a3 = def_p6a3

        		global p1a4, p2a4, p3a4, p4a4, p5a4, p6a4
        		p1a4 = def_p1a4
        		p2a4 = def_p2a4
        		p3a4 = def_p3a4
        		p4a4 = def_p4a4
        		p5a4 = def_p5a4
        		p6a4 = def_p6a4

        		global p1a5, p2a5, p3a5, p4a5, p5a5, p6a5
        		p1a5 = def_p1a5
        		p2a5 = def_p2a5
        		p3a5 = def_p3a5
        		p4a5 = def_p4a5
        		p5a5 = def_p5a5
        		p6a5 = def_p6a5

        		global p1a6, p2a6, p3a6, p4a6, p5a6, p6a6
        		p1a6 = def_p1a6
        		p2a6 = def_p2a6
        		p3a6 = def_p3a6
        		p4a6 = def_p4a6
        		p5a6 = def_p5a6
        		p6a6 = def_p6a6

        		global p1a7, p2a7, p3a7, p4a7, p5a7, p6a7
        		p1a7 = def_p1a7
        		p2a7 = def_p2a7
        		p3a7 = def_p3a7
        		p4a7 = def_p4a7
        		p5a7 = def_p5a7
        		p6a7 = def_p6a7

        		global p1a8, p2a8, p3a8, p4a8, p5a8, p6a8
        		p1a8 = def_p1a8
        		p2a8 = def_p2a8
        		p3a8 = def_p3a8
        		p4a8 = def_p4a8
        		p5a8 = def_p5a8
        		p6a8 = def_p6a8

        		global p1a9, p2a9, p3a9, p4a9, p5a9, p6a9
        		p1a9 = def_p1a9
        		p2a9 = def_p2a9
        		p3a9 = def_p3a9
        		p4a9 = def_p4a9
        		p5a9 = def_p5a9
        		p6a9 = def_p6a9

        		global p1a10, p2a10, p3a10, p4a10, p5a10, p6a10
        		p1a10 = def_p1a10
        		p2a10 = def_p2a10
        		p3a10 = def_p3a10
        		p4a10 = def_p4a10
        		p5a10 = def_p5a10
        		p6a10 = def_p6a10

        		global p1a11, p2a11, p3a11, p4a11, p5a11, p6a11
        		p1a11 = def_p1a11
        		p2a11 = def_p2a11
        		p3a11 = def_p3a11
        		p4a11 = def_p4a11
        		p5a11 = def_p5a11
        		p6a11 = def_p6a11

        		global p1a12, p2a12, p3a12, p4a12, p5a12, p6a12
        		p1a12 = def_p1a12
        		p2a12 = def_p2a12
        		p3a12 = def_p3a12
        		p4a12 = def_p4a12
        		p5a12 = def_p5a12
        		p6a12 = def_p6a12

        		print("")
	        	tbl1 = Texttable()
	        	tbl1.add_rows([[' ', pos_1, pos_2, pos_3, pos_4, pos_5, pos_6]
	        		, [all_1, p1a1, p2a1, p3a1, p4a1, p5a1, p6a1]
	        		, [all_2, p1a2, p2a2, p3a2, p4a2, p5a2, p6a2]
	        		, [all_3, p1a3, p2a3, p3a3, p4a3, p5a3, p6a3]
	        		, [all_4, p1a4, p2a4, p3a4, p4a4, p5a4, p6a4]
	        		, [all_5, p1a5, p2a5, p3a5, p4a5, p5a5, p6a5]
	        		, [all_6, p1a6, p2a6, p3a6, p4a6, p5a6, p6a6]
	        		, [all_7, p1a7, p2a7, p3a7, p4a7, p5a7, p6a7]
	        		, [all_8, p1a8, p2a8, p3a8, p4a8, p5a8, p6a8]
	        		, [all_9, p1a9, p2a9, p3a9, p4a9, p5a9, p6a9]
	        		, [all_10, p1a10, p2a10, p3a10, p4a10, p5a10, p6a10]
	        		, [all_11, p1a11, p2a11, p3a11, p4a11, p5a11, p6a11]
	        		, [all_12, p1a12, p2a12, p3a12, p4a12, p5a12, p6a12]
	        		])
	        	print(tbl1.draw())
	        	print("")

	        	# Create region
		        create_region()

        	elif ans == 'n':
        		break
        		print("")
        	else:
        		print("Not correct!")
        		print("")

        # Module complete
        print("Module 2: Load table complete")
        print("")
        break

# Example simulate VCF file ##################################################
def example_simulate_vcf_file():
    while True:

        # Module Name
        print("Module 3: Example simulate VCF file")
        print("")

        print("\tAnswers for test cases with missing positions.")
        print("")

        print("\tExample.")
        txt_1 = '\tAnswers of diplotype is: '+all_2+'/'+all_3+', '+'missing positions '+pos_2+', '+pos_4+', '+pos_5
        print(txt_1)
        txt_2 = '\t'+reg_2+' '+all_2+' '+'(Region 1)'
        print(txt_2)
        txt_3 = '\t'+reg_3+' '+all_3+' '+'(Region 2)'
        print(txt_3)
        print("")

        print("\tFind the possible star allele.")
        txt_4 = '\tEx. answers for region 1 = '+all_2+' '+'(related position'+' '+pos_2+')'
        print(txt_4)
        print("")

        print("\t1). Find all combination of missing positions.")
        x = pos_2+pos_4+pos_5
        for i in range(len(x)):
        	for l in combinations(x, i+1):
        		print('\t\t'+', '.join(l))
        print("")

        print("\t2). Find positions relate answers.")
        print('\t\t'+pos_2)
        print("")

        print("\t3). Insert positions relate answers into all combination in 1).")
        print('\t\t'+'Position relate answers = position'+pos_2)
        x = pos_2+pos_4+pos_5
        for i in range(len(x)):
        	for l in combinations(x, i+1):
        		s = '\t\t'+'2, '+(', '.join(l))
        		print(s)
        print("")

        print("\t4). From 1). select combination include positions relate answer. And remove positions relate answer.")
        print('\t\t'+'Combination include position relate answer = ')
        print('\t\t'+pos_2)
        x = pos_4+pos_5
        for i in range(len(x)):
            for l in combinations(x, i+1):
                s = '\t\t'+'2, '+(', '.join(l))
                print(s)
        print('\t\t'+'Remove position '+pos_2)

        x = pos_4+pos_5
        for i in range(len(x)):
            for l in combinations(x, i+1):
                s = '\t\t'+(', '.join(l))
                print(s)
        print("")

        print("\t5). Merge all combination from 2). 3). 4).")
        print('\t\t'+pos_2)
        print("")
        x = pos_2+pos_4+pos_5
        for i in range(len(x)):
        	for l in combinations(x, i+1):
        		s = '\t\t'+'2, '+(', '.join(l))
        		print(s)
        print("")
        x = pos_4+pos_5
        for i in range(len(x)):
            for l in combinations(x, i+1):
                s = '\t\t'+(', '.join(l))
                print(s)
        print("")

        print("\t6). Each combination from 5. remove duplicate positions.")
        print('\t\t'+pos_2)
        print("")
        x = pos_2+pos_4+pos_5
        for i in range(len(x)):
        	for l in combinations(x, i+1):
        		if '2' in l:
        			s = '\t\t'+(', '.join(l))
        		else:
        			s = '\t\t'+'2, '+(', '.join(l))
        		print(s)
        print("")
        x = pos_4+pos_5
        for i in range(len(x)):
            for l in combinations(x, i+1):
                s = '\t\t'+(', '.join(l))
                print(s)
        print("")

        print("\t7). Remove duplicate combination from 6).")
        print('\t\t'+pos_2)
        x = pos_4+pos_5
        for i in range(len(x)):
        	for l in combinations(x, i+1):
        		s = '\t\t'+'2, '+(', '.join(l))
        		print(s)
        x = pos_4+pos_5
        for i in range(len(x)):
            for l in combinations(x, i+1):
                s = '\t\t'+(', '.join(l))
                print(s)
        print("")

        print("\t8). Find star allele relate positions from 7).")

        sta_allele_2 = '2'
        sta_allele_3 = '3'
        sta_allele_4 = '4'
        sta_allele_5 = '5'
        sta_allele_6 = '6'
        sta_allele_7 = '1, 3'
        sta_allele_8 = '2, 6'
        sta_allele_9 = '3, 5'
        sta_allele_10 = '2, 4'
        sta_allele_11 = '1, 4'
        sta_allele_12 = '5, 6'

        sta_allele_2_list = ['2']
        sta_allele_3_list = ['3']
        sta_allele_4_list = ['4']
        sta_allele_5_list = ['5']
        sta_allele_6_list = ['6']
        sta_allele_7_list = ['1, 3', '3, 1']
        sta_allele_8_list = ['2, 6', '6, 2']
        sta_allele_9_list = ['3, 5', '5, 3']
        sta_allele_10_list = ['2, 4', '4, 2']
        sta_allele_11_list = ['1, 4', '4, 1']
        sta_allele_12_list = ['5, 6', '5, 5']

        all_pos_all_1 = 'a'
        all_pos_all_2 = ''
        all_pos_all_3 = ''
        all_pos_all_4 = ''
        all_pos_all_5 = ''
        all_pos_all_6 = ''
        all_pos_all_7 = ''
        all_pos_all_8 = ''
        all_pos_all_9 = ''
        all_pos_all_10 = ''
        all_pos_all_11 = ''
        all_pos_all_12 = ''

        com_1 = pos_2
        global pos_alleles_com_1_reg_1
        if com_1 in sta_allele_2_list:
            str2 = sta_allele_2+' '+'\t-->'+' '+all_2
            pos_alleles_com_1_reg_1 = all_2+' '
            all_pos_all_2 = 'b'
            print('\t\t'+str2)
        elif com_1 in sta_allele_3_list:
            str3 = sta_allele_3+' '+'\t-->'+' '+all_3
            pos_alleles_com_1_reg_1 = all_3+' '
            all_pos_all_3 = 'c'
            print('\t\t'+str3)
        elif com_1 in sta_allele_4_list:
            str4 = sta_allele_4+' '+'\t-->'+' '+all_4
            pos_alleles_com_1_reg_1 = all_4+' '
            all_pos_all_4 = 'd'
            print('\t\t'+str4)
        elif com_1 in sta_allele_5_list:
            str5 = sta_allele5+' '+'\t-->'+' '+all_5
            pos_alleles_com_1_reg_1 = all_5+' '
            all_pos_all_5 = 'e'
            print('\t\t'+str5)
        elif com_1 in sta_allele_6_list:
            str6 = sta_allele_6+' '+'\t-->'+' '+all_6
            pos_alleles_com_1_reg_1 = all_6+' '
            all_pos_all_6 = 'f'
            print('\t\t'+str6)
        elif com_1 in sta_allele_7_list:
            str7 = sta_allele_7+' '+'\t-->'+' '+all_7
            pos_alleles_com_1_reg_1 = all_7+' '
            all_pos_all_7 = 'g'
            print('\t\t'+str7)
        elif com_1 in sta_allele_8_list:
            str8 = sta_allele_8+' '+'\t-->'+' '+all_8
            pos_alleles_com_1_reg_1 = all_8+' '
            all_pos_all_8 = 'h'
            print('\t\t'+str8)
        elif com_1 in sta_allele_9_list:
            str9 = sta_allele_9+' '+'\t-->'+' '+all_9
            pos_alleles_com_1_reg_1 = all_9+' '
            all_pos_all_9 = 'i'
            print('\t\t'+str9)
        elif com_1 in sta_allele_10_list:
            str10 = sta_allele_10+' '+'\t-->'+' '+all_10
            pos_alleles_com_1_reg_1 = all_10+' '
            all_pos_all_10 = 'j'
            print('\t\t'+str10)
        elif com_1 in sta_allele_11_list:
            str11 = sta_allele_11+' '+'\t-->'+' '+all_11
            pos_alleles_com_1_reg_1 = all_11+' '
            all_pos_all_11 = 'k'
            print('\t\t'+str11)
        elif com_1 in sta_allele_12_list:
            str12 = sta_allele_12+' '+'\t-->'+' '+all_12
            pos_alleles_com_1_reg_1 = all_12+' '
            all_pos_all_12 = 'l'
            print('\t\t'+str12)
        else:
            str_none ='None'
            print('\t\t'+str_none)
        all_pos_alleles_com_1_reg_1 = pos_alleles_com_1_reg_1

        x = pos_4+pos_5
        for i in range(len(x)):
            for l in combinations(x, i+1):
                s2 = '2, '+(', '.join(l))
                com_2 = s2
                global pos_alleles_2_com_2_reg_1, pos_alleles_3_com_2_reg_1, pos_alleles_4_com_2_reg_1, pos_alleles_5_com_2_reg_1, pos_alleles_6_com_2_reg_1, pos_alleles_7_com_2_reg_1, pos_alleles_8_com_2_reg_1, pos_alleles_9_com_2_reg_1, pos_alleles_10_com_2_reg_1, pos_alleles_11_com_2_reg_1, pos_alleles_12_com_2_reg_1, pos_alleles_none_com_2_reg_1
                if com_2 in sta_allele_2_list:
                    str2 = com_2+' '+'\t-->'+' '+all_2
                    pos_alleles_2_com_2_reg_1 = all_2+' '
                    all_pos_all_2 = 'b'
                    print('\t\t'+str2)
                elif com_2 in sta_allele_3_list:
                    str3 = com_2+' '+'\t-->'+' '+all_3
                    pos_alleles_3_com_2_reg_1 = all_3+' '
                    all_pos_all_3 = 'c'
                    print('\t\t'+str3)
                elif com_2 in sta_allele_4_list:
                    str4 = com_2+' '+'\t-->'+' '+all_4
                    pos_alleles_4_com_2_reg_1 = all_4+' '
                    all_pos_all_4 = 'd'
                    print('\t\t'+str4)
                elif com_2 in sta_allele_5_list:
                    str5 = com_2+' '+'\t-->'+' '+all_5
                    pos_alleles_5_com_2_reg_1 = all_5+' '
                    all_pos_all_5 = 'e'
                    print('\t\t'+str5)
                elif com_2 in sta_allele_6_list:
                    str6 = com_2+' '+'\t-->'+' '+all_6
                    pos_alleles_6_com_2_reg_1 = all_6+' '
                    all_pos_all_6 = 'f'
                    print('\t\t'+str6)
                elif com_2 in sta_allele_7_list:
                    str7 = com_2+' '+'\t-->'+' '+all_7
                    pos_alleles_7_com_2_reg_1 = all_7+' '
                    all_pos_all_7 = 'g'
                    print('\t\t'+str7)
                elif com_2 in sta_allele_8_list:
                    str8 = com_2+' '+'\t-->'+' '+all_8
                    pos_alleles_8_com_2_reg_1 = all_8+' '
                    all_pos_all_8 = 'h'
                    print('\t\t'+str8)
                elif com_2 in sta_allele_9_list:
                    str9 = com_2+' '+'\t-->'+' '+all_9
                    pos_alleles_9_com_2_reg_1 = all_9+' '
                    all_pos_all_9 = 'i'
                    print('\t\t'+str9)
                elif com_2 in sta_allele_10_list:
                    str10 = com_2+' '+'\t-->'+' '+all_10
                    pos_alleles_10_com_2_reg_1 = all_10+' '
                    all_pos_all_10 = 'j'
                    print('\t\t'+str10)
                elif com_2 in sta_allele_11_list:
                    str11 = com_2+' '+'\t-->'+' '+all_11
                    pos_alleles_11_com_2_reg_1 = all_11+' '
                    all_pos_all_11 = 'k'
                    print('\t\t'+str11)
                elif com_2 in sta_allele_12_list:
                    str12 = com_2+' '+'\t-->'+' '+all_12
                    pos_alleles_12_com_2_reg_1 = all_12+' '
                    all_pos_all_12 = 'l'
                    print('\t\t'+str12)
                else:
                    str_none = com_2+' '+'\t-->'+' '+'None'
                    pos_alleles_none_com_2_reg_1 = ''
                    print('\t\t'+str_none)

        all_pos_alleles_com_2_reg_1 = pos_alleles_2_com_2_reg_1+pos_alleles_3_com_2_reg_1+pos_alleles_4_com_2_reg_1+pos_alleles_5_com_2_reg_1+pos_alleles_6_com_2_reg_1+pos_alleles_7_com_2_reg_1+pos_alleles_8_com_2_reg_1+pos_alleles_9_com_2_reg_1+pos_alleles_10_com_2_reg_1+pos_alleles_11_com_2_reg_1+pos_alleles_12_com_2_reg_1+pos_alleles_none_com_2_reg_1

        x = pos_4+pos_5
        for i in range(len(x)):
            for l in combinations(x, i+1):
                s3 = ', '.join(l)
                com_3 = s3
                global pos_alleles_2_com_3_reg_1, pos_alleles_3_com_3_reg_1, pos_alleles_4_com_3_reg_1, pos_alleles_5_com_3_reg_1, pos_alleles_6_com_3_reg_1, pos_alleles_7_com_3_reg_1, pos_alleles_8_com_3_reg_1, pos_alleles_9_com_3_reg_1, pos_alleles_10_com_3_reg_1, pos_alleles_11_com_3_reg_1, pos_alleles_12_com_3_reg_1, pos_alleles_none_com_3_reg_1
                if com_3 in sta_allele_2_list:
                    str2 = com_3+' '+'\t-->'+' '+all_2
                    pos_alleles_2_com_3_reg_1 = all_2+' '
                    all_pos_all_2 = 'b'
                    print('\t\t'+str2)
                elif com_3 in sta_allele_3_list:
                    str3 = com_3+' '+'\t-->'+' '+all_3
                    pos_alleles_3_com_3_reg_1 = all_3+' '
                    all_pos_all_3 = 'c'
                    print('\t\t'+str3)
                elif com_3 in sta_allele_4_list:
                    str4 = com_3+' '+'\t-->'+' '+all_4
                    pos_alleles_4_com_3_reg_1 = all_4+' '
                    all_pos_all_4 = 'd'
                    print('\t\t'+str4)
                elif com_3 in sta_allele_5_list:
                    str5 = com_3+' '+'\t-->'+' '+all_5
                    pos_alleles_5_com_3_reg_1 = all_5+' '
                    all_pos_all_5 = 'e'
                    print('\t\t'+str5)
                elif com_3 in sta_allele_6_list:
                    str6 = com_3+' '+'\t-->'+' '+all_6
                    pos_alleles_6_com_3_reg_1 = all_6+' '
                    all_pos_all_6 = 'f'
                    print('\t\t'+str6)
                elif com_3 in sta_allele_7_list:
                    str7 = com_3+' '+'\t-->'+' '+all_7
                    pos_alleles_7_com_3_reg_1 = all_7+' '
                    all_pos_all_7 = 'g'
                    print('\t\t'+str7)
                elif com_3 in sta_allele_8_list:
                    str8 = com_3+' '+'\t-->'+' '+all_8
                    pos_alleles_8_com_3_reg_1 = all_8+' '
                    all_pos_all_8 = 'h'
                    print('\t\t'+str8)
                elif com_3 in sta_allele_9_list:
                    str9 = com_3+' '+'\t-->'+' '+all_9
                    pos_alleles_9_com_3_reg_1 = all_9+' '
                    all_pos_all_9 = 'i'
                    print('\t\t'+str9)
                elif com_3 in sta_allele_10_list:
                    str10 = com_3+' '+'\t-->'+' '+all_10
                    pos_alleles_10_com_3_reg_1 = all_10+' '
                    all_pos_all_10 = 'j'
                    print('\t\t'+str10)
                elif com_3 in sta_allele_11_list:
                    str11 = com_3+' '+'\t-->'+' '+all_11
                    pos_alleles_11_com_3_reg_1 = all_11+' '
                    all_pos_all_11 = 'k'
                    print('\t\t'+str11)
                elif com_3 in sta_allele_12_list:
                    str12 = com_3+' '+'\t-->'+' '+all_12
                    pos_alleles_12_com_3_reg_1 = all_12+' '
                    all_pos_all_12 = 'l'
                    print('\t\t'+str12)
                else:
                    str_none = com_3+' '+'\t-->'+' '+'None'
                    pos_alleles_none_com_3_reg_1 = ''
                    print('\t\t'+str_none)
        
        print("")

        all_pos_alleles_com_3_reg_1 = pos_alleles_2_com_3_reg_1+pos_alleles_3_com_3_reg_1+pos_alleles_4_com_3_reg_1+pos_alleles_5_com_3_reg_1+pos_alleles_6_com_3_reg_1+pos_alleles_7_com_3_reg_1+pos_alleles_8_com_3_reg_1+pos_alleles_9_com_3_reg_1+pos_alleles_10_com_3_reg_1+pos_alleles_11_com_3_reg_1+pos_alleles_12_com_3_reg_1+pos_alleles_none_com_3_reg_1
        txt1  = all_pos_all_2+all_pos_all_3+all_pos_all_4+all_pos_all_5+all_pos_all_6+all_pos_all_7+all_pos_all_8+all_pos_all_9+all_pos_all_10+all_pos_all_11+all_pos_all_12

        cc = []

        for i in txt1:
            if i == 'a':
                i = '*1'
            elif i == 'b':
                i = '*2'
            elif i == 'c':
                i = '*3'
            elif i == 'd':
                i = '*4'
            elif i == 'e':
                i = '*5'
            elif i == 'f':
                i = '*6'
            elif i == 'g':
                i = '*7'
            elif i == 'h':
                i = '*8'
            elif i == 'i':
                i = '*9'
            elif i == 'j':
                i = '*10'
            elif i == 'k':
                i = '*11'
            elif i == 'l':
                i = '*12'
            cc.append(i)

        txt2  = '\t\tAll possible alleles for answer of missing in region 1. = '
        print(txt2+', '.join(cc))

        print("")

        print("\t9). If all position relate answer is missing, add *1.")
        txt1 = ''
        all_pos_miss = pos_2+', '+pos_4+', '+pos_5
        if pos_2 in all_pos_miss:
            txt1 = all_pos_all_1+all_pos_all_2+all_pos_all_3+all_pos_all_4+all_pos_all_5+all_pos_all_6+all_pos_all_7+all_pos_all_8+all_pos_all_9+all_pos_all_10+all_pos_all_11+all_pos_all_12
            txt2 = '\t\tAll position relate answer = '+pos_2+', '+'that in list missing position '+pos_2+', '+pos_4+', '+pos_5
            
            dd = []

            for i in txt1:
                if i == 'a':
                    i = '*1'
                elif i == 'b':
                    i = '*2'
                elif i == 'c':
                    i = '*3'
                elif i == 'd':
                    i = '*4'
                elif i == 'e':
                    i = '*5'
                elif i == 'f':
                    i = '*6'
                elif i == 'g':
                    i = '*7'
                elif i == 'h':
                    i = '*8'
                elif i == 'i':
                    i = '*9'
                elif i == 'j':
                    i = '*10'
                elif i == 'k':
                    i = '*11'
                elif i == 'l':
                    i = '*12'
                dd.append(i)

            txt3  = '\t\tSo, all possible alleles for answer of missing in region 1. = '
            print(txt2)
            print("")
            print(txt3+', '.join(dd))
        else:
            txt1 = all_pos_all_2+all_pos_all_3+all_pos_all_4+all_pos_all_5+all_pos_all_6+all_pos_all_7+all_pos_all_8+all_pos_all_9+all_pos_all_10+all_pos_all_11+all_pos_all_12
            txt2 = '\t\tAll position relate answer = '+pos_2+', '+'that not in list missing position '+pos_2+', '+pos_4+', '+pos_5
            
            dd = []

            for i in txt1:
                if i == 'a':
                    i = '*1'
                elif i == 'b':
                    i = '*2'
                elif i == 'c':
                    i = '*3'
                elif i == 'd':
                    i = '*4'
                elif i == 'e':
                    i = '*5'
                elif i == 'f':
                    i = '*6'
                elif i == 'g':
                    i = '*7'
                elif i == 'h':
                    i = '*8'
                elif i == 'i':
                    i = '*9'
                elif i == 'j':
                    i = '*10'
                elif i == 'k':
                    i = '*11'
                elif i == 'l':
                    i = '*12'
                dd.append(i)

            txt3 = '\t\tSo, all possible alleles for answer of missing in region 1. = '
            print(txt2)
            print("")
            print(txt3+', '.join(dd))

        print("")

        print("\t10). Repeat from 2). to 9). for region 2.")
        all_pos_all_3 = 'c'
        all_pos_all_9 = 'i'
        txt6 = all_pos_all_3+all_pos_all_9
        
        ee = []

        for i in txt6:
            if i == 'a':
                i = '*1'
            elif i == 'b':
                i = '*2'
            elif i == 'c':
                i = '*3'
            elif i == 'd':
                i = '*4'
            elif i == 'e':
                i = '*5'
            elif i == 'f':
                i = '*6'
            elif i == 'g':
                i = '*7'
            elif i == 'h':
                i = '*8'
            elif i == 'i':
                i = '*9'
            elif i == 'j':
                i = '*10'
            elif i == 'k':
                i = '*11'
            elif i == 'l':
                i = '*12'
            ee.append(i)

        txt7 = '\t\tSo, all possible alleles for answer of missing in region 2. = '
        print(txt7+', '.join(ee))
        print("")

        print("\t11). From 9). and 10). find all diplotypes combination of allele in each region.")
        print('\t\t'+reg_2+' '+', '.join(dd))
        print('\t\t'+reg_3+' '+', '.join(ee))
        print("")
        result = []
        for i in txt1:
            if i == 'a':
                i = '*1'
            elif i == 'b':
                i = '*2'
            elif i == 'c':
                i = '*3'
            elif i == 'd':
                i = '*4'
            elif i == 'e':
                i = '*5'
            elif i == 'f':
                i = '*6'
            elif i == 'g':
                i = '*7'
            elif i == 'h':
                i = '*8'
            elif i == 'i':
                i = '*9'
            elif i == 'j':
                i = '*10'
            elif i == 'k':
                i = '*11'
            elif i == 'l':
                i = '*12'
            for j in txt6:
                if j == 'a':
                    j = '*1'
                elif j == 'b':
                    j = '*2'
                elif j == 'c':
                    j = '*3'
                elif j == 'd':
                    j = '*4'
                elif j == 'e':
                    j = '*5'
                elif j == 'f':
                    j = '*6'
                elif j == 'g':
                    j = '*7'
                elif j == 'h':
                    j = '*8'
                elif j == 'i':
                    j = '*9'
                elif j == 'j':
                    j = '*10'
                elif j == 'k':
                    j = '*11'
                elif j == 'l':
                    j = '*12'
                result.append((i+'/'+ j))
        txt8 = '\t\tAll possible diplotypes = '
        values = ', '.join(str(v) for v in result)
        print(txt8+' '+values)

        print("")
        # Module complete
        print("Module 3: Example simulate VCF file complete")
        print("")
        break

# Example test cases ##################################################
def example_test_cases():
    while True:

        # Module name
        print("Module 4: Example test cases")
        print("")

        print("\tExample case for answer of sample diplotype "+all_2+'/'+all_3)
        print("")
        print('\t'+reg_2+' '+all_2+' '+'(Region 1)')
        print('\t'+reg_3+' '+all_3+' '+'(Region 2)')
        print("")

        print("\tMissing "+pos_1+" Position")
        print("")
        print("\t\t1). missing position "+pos_2)
        print('\t\t\tRegion 1 '+' '+all_1+', '+all_2)
        print('\t\t\tRegion 2 '+' '+all_3)
        print("\t\t2). missing position "+pos_3)
        print('\t\t\tRegion 1 '+' '+all_2)
        print('\t\t\tRegion 2 '+' '+all_1+', '+all_3)
        print("\t\t3). missing position "+pos_4)
        print('\t\t\tRegion 1 '+' '+all_2+', '+all_10)
        print('\t\t\tRegion 2 '+' '+all_3)
        print("\t\t4). missing position "+pos_5)
        print('\t\t\tRegion 1 '+' '+all_2)
        print('\t\t\tRegion 2 '+' '+all_3+', '+all_9)
        print("\t\t5). missing position "+pos_6)
        print('\t\t\tRegion 1 '+' '+all_2+', '+all_8)
        print('\t\t\tRegion 2 '+' '+all_3)
        print("")

        print("\tMissing "+pos_2+" Position")
        print("")
        print("\t\t1). missing position "+pos_2+', '+pos_3)
        print('\t\t\tRegion 1 '+' '+all_1+', '+all_2+', '+all_3)
        print('\t\t\tRegion 2 '+' '+all_1+', '+all_2+', '+all_3)
        print("\t\t2). missing position "+pos_2+', '+pos_4)
        print('\t\t\tRegion 1 '+' '+all_1+', '+all_2+', '+all_4+', '+all_10)
        print('\t\t\tRegion 2 '+' '+all_3)
        print("\t\t3). missing position "+pos_3+', '+pos_5)
        print('\t\t\tRegion 1 '+' '+all_2)
        print('\t\t\tRegion 2 '+' '+all_1+', '+all_3+', '+all_5+', '+all_9)
        print("\t\t4). missing position "+pos_1+', '+pos_4)
        print('\t\t\tRegion 1 '+' '+all_2+', '+all_10)
        print('\t\t\tRegion 2 '+' '+all_3+', '+all_7)
        print("")

        print("\tMissing "+pos_3+" Position")
        print("")
        print("\t\t1). missing position "+pos_2+', '+pos_3+', '+pos_4)
        print('\t\t\tRegion 1 '+' '+all_1+', '+all_2+', '+all_3+', '+all_4+', '+all_10)
        print('\t\t\tRegion 2 '+' '+all_1+', '+all_2+', '+all_3+', '+all_4+', '+all_10)
        print("\t\t2). missing position "+pos_2+', '+pos_4+', '+pos_5)
        print('\t\t\tRegion 1 '+' '+all_1+', '+all_2+', '+all_4+', '+all_5+', '+all_10)
        print('\t\t\tRegion 2 '+' '+all_3+', '+all_9)
        print("\t\t3). missing position "+pos_3+', '+pos_4+', '+pos_5)
        print('\t\t\tRegion 1 '+' '+all_2+', '+all_10)
        print('\t\t\tRegion 2 '+' '+all_1+', '+all_3+', '+all_4+', '+all_5+', '+all_9)
        print("\t\t4). missing position "+pos_4+', '+pos_5+', '+pos_6)
        print('\t\t\tRegion 1 '+' '+all_2+', '+all_8+', '+all_10)
        print('\t\t\tRegion 2 '+' '+all_3+', '+all_9)
        print("")

        print("\tMissing "+pos_4+" Position")
        print("")
        print("\t\t1). missing position "+pos_2+', '+pos_3+', '+pos_4+', '+pos_5)
        print('\t\t\tRegion 1 '+' '+all_1+', '+all_2+', '+all_3+', '+all_4+', '+all_5+', '+all_9+', '+all_10)
        print('\t\t\tRegion 2 '+' '+all_1+', '+all_2+', '+all_3+', '+all_4+', '+all_5+', '+all_9+', '+all_10)
        print("")

        # Module complete
        print("Module 4: Example test cases complete")
        print("")
        break

# Find answers ##################################################
def find_answers():
    while True:

        # Module name
        print("Module 5: Find answers cases")
        print("")

        input_ans_dip_2 = '*2'
        input_ans_dip_3 = '*3'

        print("\tPlease enter answers of diplotype.")
        print("")
        print("\tInput answers of diplotype region 1: "+input_ans_dip_2)
        print("\tInput answers of diplotype region 2: "+input_ans_dip_3)

        input_ans_dip_reg_1 = input_ans_dip_2
        input_ans_dip_reg_2 = input_ans_dip_3

        print("")

        ##### Start: Input missing position #####

        print("\tPlease enter missing position.")
        print("")

        global mis_pos_1
        input_mis_pos_1 = ''
        while input_mis_pos_1 not in {'y', 'n'}:
            input_mis_pos_1 = str(input("\tMissing position 1 ...? y/n: "))
            if input_mis_pos_1 == 'y':
                mis_pos_1 = '1'
            elif input_mis_pos_1 == 'n':
                mis_pos_1 = ''
            else:
                print("")
                print("\tNot correct!")
                print("")

        global mis_pos_2
        input_mis_pos_2 = ''
        while input_mis_pos_2 not in {'y', 'n'}:
            input_mis_pos_2 = str(input("\tMissing position 2 ...? y/n: "))
            if input_mis_pos_2 == 'y':
                mis_pos_2 = '2'
            elif input_mis_pos_2 == 'n':
                mis_pos_2 = ''
            else:
                print("")
                print("\tNot correct!")
                print("")

        global mis_pos_3
        input_mis_pos_3 = ''
        while input_mis_pos_3 not in {'y', 'n'}:
            input_mis_pos_3 = str(input("\tMissing position 3 ...? y/n: "))
            if input_mis_pos_3 == 'y':
                mis_pos_3 = '3'
            elif input_mis_pos_3 == 'n':
                mis_pos_3 = ''
            else:
                print("")
                print("\tNot correct!")
                print("")

        global mis_pos_4
        input_mis_pos_4 = ''
        while input_mis_pos_4 not in {'y', 'n'}:
            input_mis_pos_4 = str(input("\tMissing position 4 ...? y/n: "))
            if input_mis_pos_4 == 'y':
                mis_pos_4 = '4'
            elif input_mis_pos_4 == 'n':
                mis_pos_4 = ''
            else:
                print("")
                print("\tNot correct!")
                print("")

        global mis_pos_5
        input_mis_pos_5 = ''
        while input_mis_pos_5 not in {'y', 'n'}:
            input_mis_pos_5 = str(input("\tMissing position 5 ...? y/n: "))
            if input_mis_pos_5 == 'y':
                mis_pos_5 = '5'
            elif input_mis_pos_5 == 'n':
                mis_pos_5 = ''
            else:
                print("")
                print("\tNot correct!")
                print("")

        global mis_pos_6
        input_mis_pos_6 = ''
        while input_mis_pos_6 not in {'y', 'n'}:
            input_mis_pos_6 = str(input("\tMissing position 6 ...? y/n: "))
            if input_mis_pos_6 == 'y':
                mis_pos_6 = '6'
            elif input_mis_pos_6 == 'n':
                mis_pos_6 = ''
            else:
                print("")
                print("\tNot correct!")
                print("")

        input_mis_pos = mis_pos_1+mis_pos_2+mis_pos_3+mis_pos_4+mis_pos_5+mis_pos_6
        input_mis_pos_list = []
        for i in input_mis_pos:
            input_mis_pos_list.append(i)

        print("")
        txt1  = '\tInput missing position: '
        print(txt1+', '.join(input_mis_pos_list))
        print("")

        input_calulate = ''
        while input_calulate not in {'y'}:
            input_calulate = str(input("\tCalculate possible diplotype? y: "))
            if input_calulate == 'y':
                None
            else:
                print("")
                print("\tNot correct!")
                print("")

        ##### End: Input missing position #####

        print("")

        print("\t1). Find all combination of missing positions.")
        x = input_mis_pos
        for i in range(len(x)):
            for l in combinations(x, i+1):
                print('\t\t'+', '.join(l))
        print("")

        print("\t2/1). Find positions relate answers.")
        if input_ans_dip_reg_1 == '*2':
            pos_rel_ans = '2'
        elif input_ans_dip_reg_1 == '*3':
            pos_rel_ans = '3'
        print('\t\t'+pos_rel_ans)
        print("")

        print("\t3/1). Insert positions relate answers into all combination in 1).")
        x = input_mis_pos
        for i in range(len(x)):
            for l in combinations(x, i+1):
                s = '\t\t'+pos_rel_ans+', '+(', '.join(l))
                print(s)
        print("")

        print("\t4/1). From 1). select combination include positions relate answer. And remove positions relate answer.")
        if pos_rel_ans in input_mis_pos:
            x = input_mis_pos.replace(pos_rel_ans,'')
            for i in range(len(x)):
                for l in combinations(x, i+1):
                    print('\t\t'+', '.join(l))

        else:
            #x = 'None'
            #print('\t\t'+x)
            None
        print("")

        print("\t5/1). Merge all combination from 2/1). 3/1). 4/1).")
        if input_ans_dip_reg_1 == '*2':
            pos_rel_ans = '2'
        elif input_ans_dip_reg_1 == '*3':
            pos_rel_ans = '3'
        print('\t\t'+pos_rel_ans)
        print("")

        x = input_mis_pos
        for i in range(len(x)):
            for l in combinations(x, i+1):
                s = '\t\t'+pos_rel_ans+', '+(', '.join(l))
                print(s)
        print("")

        if pos_rel_ans in input_mis_pos:
            x = input_mis_pos.replace(pos_rel_ans,'')
            for i in range(len(x)):
                for l in combinations(x, i+1):
                    print('\t\t'+', '.join(l))
        else:
            #x = 'None'
            #print('\t\t'+x)
            None
        print("")

        print("\t6/1). Each combination from 5/1). remove duplicate positions.")
        if input_ans_dip_reg_1 == '*2':
            pos_rel_ans = '2'
        elif input_ans_dip_reg_1 == '*3':
            pos_rel_ans = '3'
        print('\t\t'+pos_rel_ans)
        print("")

        x = input_mis_pos
        for i in range(len(x)):
            for l in combinations(x, i+1):
                if pos_rel_ans in l:
                    s = '\t\t'+(', '.join(l))
                else:
                    s = '\t\t'+pos_rel_ans+', '+(', '.join(l))
                print(s)
        print("")

        if pos_rel_ans in input_mis_pos:
            x = input_mis_pos.replace(pos_rel_ans,'')
            for i in range(len(x)):
                for l in combinations(x, i+1):
                    print('\t\t'+', '.join(l))
        else:
            #x = 'None'
            #print('\t\t'+x)
            None
        print("")

        print("\t7/1). Remove duplicate combination from 6/1).")
        if input_ans_dip_reg_1 == '*2':
            pos_rel_ans = '2'
        elif input_ans_dip_reg_1 == '*3':
            pos_rel_ans = '3'
        print('\t\t'+pos_rel_ans)

        if pos_rel_ans in input_mis_pos:
            x = input_mis_pos.replace(pos_rel_ans,'')
            for i in range(len(x)):
                for l in combinations(x, i+1):
                    s = '\t\t'+pos_rel_ans+', '+(', '.join(l))
                    print(s)
        else:
            x = input_mis_pos.replace(pos_rel_ans,'')
            for i in range(len(x)):
                for l in combinations(x, i+1):
                    s = '\t\t'+pos_rel_ans+', '+(', '.join(l))
                    print(s)

        if pos_rel_ans in input_mis_pos:
            x = input_mis_pos.replace(pos_rel_ans,'')
            for i in range(len(x)):
                for l in combinations(x, i+1):
                    print('\t\t'+', '.join(l))
        else:
            #x = 'None'
            #print('\t\t'+x)
            None
        print("")

        print("\t8/1). Find star allele relate positions from 7/1).")
        sta_allele_2 = '2'
        sta_allele_3 = '3'
        sta_allele_4 = '4'
        sta_allele_5 = '5'
        sta_allele_6 = '6'
        sta_allele_7 = '1, 3'
        sta_allele_8 = '2, 6'
        sta_allele_9 = '3, 5'
        sta_allele_10 = '2, 4'
        sta_allele_11 = '1, 4'
        sta_allele_12 = '5, 6'

        sta_allele_2_list = ['2']
        sta_allele_3_list = ['3']
        sta_allele_4_list = ['4']
        sta_allele_5_list = ['5']
        sta_allele_6_list = ['6']
        sta_allele_7_list = ['1, 3', '3, 1']
        sta_allele_8_list = ['2, 6', '6, 2']
        sta_allele_9_list = ['3, 5', '5, 3']
        sta_allele_10_list = ['2, 4', '4, 2']
        sta_allele_11_list = ['1, 4', '4, 1']
        sta_allele_12_list = ['5, 6', '5, 5']

        all_pos_all_1 = 'a'
        all_pos_all_2 = ''
        all_pos_all_3 = ''
        all_pos_all_4 = ''
        all_pos_all_5 = ''
        all_pos_all_6 = ''
        all_pos_all_7 = ''
        all_pos_all_8 = ''
        all_pos_all_9 = ''
        all_pos_all_10 = ''
        all_pos_all_11 = ''
        all_pos_all_12 = ''

        if input_ans_dip_reg_1 == '*2':
            pos_rel_ans = '2'
        elif input_ans_dip_reg_1 == '*3':
            pos_rel_ans = '3'
        com_1 = pos_rel_ans
        global pos_alleles_com_1_reg_1, pos_alleles_none_com_1_reg_1
        if com_1 in sta_allele_2_list:
            str2 = sta_allele_2+' '+'\t-->'+' '+all_2
            pos_alleles_com_1_reg_1 = all_2+' '
            all_pos_all_2 = 'b'
            print('\t\t'+str2)
        elif com_1 in sta_allele_3_list:
            str3 = sta_allele_3+' '+'\t-->'+' '+all_3
            pos_alleles_com_1_reg_1 = all_3+' '
            all_pos_all_3 = 'c'
            print('\t\t'+str3)
        elif com_1 in sta_allele_4_list:
            str4 = sta_allele_4+' '+'\t-->'+' '+all_4
            pos_alleles_com_1_reg_1 = all_4+' '
            all_pos_all_4 = 'd'
            print('\t\t'+str4)
        elif com_1 in sta_allele_5_list:
            str5 = sta_allele5+' '+'\t-->'+' '+all_5
            pos_alleles_com_1_reg_1 = all_5+' '
            all_pos_all_5 = 'e'
            print('\t\t'+str5)
        elif com_1 in sta_allele_6_list:
            str6 = sta_allele_6+' '+'\t-->'+' '+all_6
            pos_alleles_com_1_reg_1 = all_6+' '
            all_pos_all_6 = 'f'
            print('\t\t'+str6)
        elif com_1 in sta_allele_7_list:
            str7 = sta_allele_7+' '+'\t-->'+' '+all_7
            pos_alleles_com_1_reg_1 = all_7+' '
            all_pos_all_7 = 'g'
            print('\t\t'+str7)
        elif com_1 in sta_allele_8_list:
            str8 = sta_allele_8+' '+'\t-->'+' '+all_8
            pos_alleles_com_1_reg_1 = all_8+' '
            all_pos_all_8 = 'h'
            print('\t\t'+str8)
        elif com_1 in sta_allele_9_list:
            str9 = sta_allele_9+' '+'\t-->'+' '+all_9
            pos_alleles_com_1_reg_1 = all_9+' '
            all_pos_all_9 = 'i'
            print('\t\t'+str9)
        elif com_1 in sta_allele_10_list:
            str10 = sta_allele_10+' '+'\t-->'+' '+all_10
            pos_alleles_com_1_reg_1 = all_10+' '
            all_pos_all_10 = 'j'
            print('\t\t'+str10)
        elif com_1 in sta_allele_11_list:
            str11 = sta_allele_11+' '+'\t-->'+' '+all_11
            pos_alleles_com_1_reg_1 = all_11+' '
            all_pos_all_11 = 'k'
            print('\t\t'+str11)
        elif com_1 in sta_allele_12_list:
            str12 = sta_allele_12+' '+'\t-->'+' '+all_12
            pos_alleles_com_1_reg_1 = all_12+' '
            all_pos_all_12 = 'l'
            print('\t\t'+str12)
        else:
            str_none = com_1+' '+'\t-->'+' '+'None'
            pos_alleles_none_com_1_reg_1 = ''
            print('\t\t'+str_none)
        all_pos_alleles_com_1_reg_1 = pos_alleles_com_1_reg_1+pos_alleles_none_com_1_reg_1

        if pos_rel_ans in input_mis_pos:
            x = input_mis_pos.replace(pos_rel_ans,'')
            for i in range(len(x)):
                for l in combinations(x, i+1):
                    s = pos_rel_ans+', '+(', '.join(l))
                    com_2 = s
                    global pos_alleles_2_com_2_reg_1, pos_alleles_3_com_2_reg_1, pos_alleles_4_com_2_reg_1, pos_alleles_5_com_2_reg_1, pos_alleles_6_com_2_reg_1, pos_alleles_7_com_2_reg_1, pos_alleles_8_com_2_reg_1, pos_alleles_9_com_2_reg_1, pos_alleles_10_com_2_reg_1, pos_alleles_11_com_2_reg_1, pos_alleles_12_com_2_reg_1, pos_alleles_none_com_2_reg_1
                    if com_2 in sta_allele_2_list:
                        str2 = com_2+' '+'\t-->'+' '+all_2
                        pos_alleles_2_com_2_reg_1 = all_2+' '
                        all_pos_all_2 = 'b'
                        print('\t\t'+str2)
                    elif com_2 in sta_allele_3_list:
                        str3 = com_2+' '+'\t-->'+' '+all_3
                        pos_alleles_3_com_2_reg_1 = all_3+' '
                        all_pos_all_3 = 'c'
                        print('\t\t'+str3)
                    elif com_2 in sta_allele_4_list:
                        str4 = com_2+' '+'\t-->'+' '+all_4
                        pos_alleles_4_com_2_reg_1 = all_4+' '
                        all_pos_all_4 = 'd'
                        print('\t\t'+str4)
                    elif com_2 in sta_allele_5_list:
                        str5 = com_2+' '+'\t-->'+' '+all_5
                        pos_alleles_5_com_2_reg_1 = all_5+' '
                        all_pos_all_5 = 'e'
                        print('\t\t'+str5)
                    elif com_2 in sta_allele_6_list:
                        str6 = com_2+' '+'\t-->'+' '+all_6
                        pos_alleles_6_com_2_reg_1 = all_6+' '
                        all_pos_all_6 = 'f'
                        print('\t\t'+str6)
                    elif com_2 in sta_allele_7_list:
                        str7 = com_2+' '+'\t-->'+' '+all_7
                        pos_alleles_7_com_2_reg_1 = all_7+' '
                        all_pos_all_7 = 'g'
                        print('\t\t'+str7)
                    elif com_2 in sta_allele_8_list:
                        str8 = com_2+' '+'\t-->'+' '+all_8
                        pos_alleles_8_com_2_reg_1 = all_8+' '
                        all_pos_all_8 = 'h'
                        print('\t\t'+str8)
                    elif com_2 in sta_allele_9_list:
                        str9 = com_2+' '+'\t-->'+' '+all_9
                        pos_alleles_9_com_2_reg_1 = all_9+' '
                        all_pos_all_9 = 'i'
                        print('\t\t'+str9)
                    elif com_2 in sta_allele_10_list:
                        str10 = com_2+' '+'\t-->'+' '+all_10
                        pos_alleles_10_com_2_reg_1 = all_10+' '
                        all_pos_all_10 = 'j'
                        print('\t\t'+str10)
                    elif com_2 in sta_allele_11_list:
                        str11 = com_2+' '+'\t-->'+' '+all_11
                        pos_alleles_11_com_2_reg_1 = all_11+' '
                        all_pos_all_11 = 'k'
                        print('\t\t'+str11)
                    elif com_2 in sta_allele_12_list:
                        str12 = com_2+' '+'\t-->'+' '+all_12
                        pos_alleles_12_com_2_reg_1 = all_12+' '
                        all_pos_all_12 = 'l'
                        print('\t\t'+str12)
                    else:
                        str_none = com_2+' '+'\t-->'+' '+'None'
                        pos_alleles_none_com_2_reg_1 = ''
                        print('\t\t'+str_none)
        else:
            x = input_mis_pos.replace(pos_rel_ans,'')
            for i in range(len(x)):
                for l in combinations(x, i+1):
                    s = pos_rel_ans+', '+(', '.join(l))
                    com_2 = s
                    if com_2 in sta_allele_2_list:
                        str2 = com_2+' '+'\t-->'+' '+all_2
                        pos_alleles_2_com_2_reg_1 = all_2+' '
                        all_pos_all_2 = 'b'
                        print('\t\t'+str2)
                    elif com_2 in sta_allele_3_list:
                        str3 = com_2+' '+'\t-->'+' '+all_3
                        pos_alleles_3_com_2_reg_1 = all_3+' '
                        all_pos_all_3 = 'c'
                        print('\t\t'+str3)
                    elif com_2 in sta_allele_4_list:
                        str4 = com_2+' '+'\t-->'+' '+all_4
                        pos_alleles_4_com_2_reg_1 = all_4+' '
                        all_pos_all_4 = 'd'
                        print('\t\t'+str4)
                    elif com_2 in sta_allele_5_list:
                        str5 = com_2+' '+'\t-->'+' '+all_5
                        pos_alleles_5_com_2_reg_1 = all_5+' '
                        all_pos_all_5 = 'e'
                        print('\t\t'+str5)
                    elif com_2 in sta_allele_6_list:
                        str6 = com_2+' '+'\t-->'+' '+all_6
                        pos_alleles_6_com_2_reg_1 = all_6+' '
                        all_pos_all_6 = 'f'
                        print('\t\t'+str6)
                    elif com_2 in sta_allele_7_list:
                        str7 = com_2+' '+'\t-->'+' '+all_7
                        pos_alleles_7_com_2_reg_1 = all_7+' '
                        all_pos_all_7 = 'g'
                        print('\t\t'+str7)
                    elif com_2 in sta_allele_8_list:
                        str8 = com_2+' '+'\t-->'+' '+all_8
                        pos_alleles_8_com_2_reg_1 = all_8+' '
                        all_pos_all_8 = 'h'
                        print('\t\t'+str8)
                    elif com_2 in sta_allele_9_list:
                        str9 = com_2+' '+'\t-->'+' '+all_9
                        pos_alleles_9_com_2_reg_1 = all_9+' '
                        all_pos_all_9 = 'i'
                        print('\t\t'+str9)
                    elif com_2 in sta_allele_10_list:
                        str10 = com_2+' '+'\t-->'+' '+all_10
                        pos_alleles_10_com_2_reg_1 = all_10+' '
                        all_pos_all_10 = 'j'
                        print('\t\t'+str10)
                    elif com_2 in sta_allele_11_list:
                        str11 = com_2+' '+'\t-->'+' '+all_11
                        pos_alleles_11_com_2_reg_1 = all_11+' '
                        all_pos_all_11 = 'k'
                        print('\t\t'+str11)
                    elif com_2 in sta_allele_12_list:
                        str12 = com_2+' '+'\t-->'+' '+all_12
                        pos_alleles_12_com_2_reg_1 = all_12+' '
                        all_pos_all_12 = 'l'
                        print('\t\t'+str12)
                    else:
                        str_none = com_2+' '+'\t-->'+' '+'None'
                        pos_alleles_none_com_2_reg_1 = ''
                        print('\t\t'+str_none)
        all_pos_alleles_com_2_reg_1 = pos_alleles_2_com_2_reg_1+pos_alleles_3_com_2_reg_1+pos_alleles_4_com_2_reg_1+pos_alleles_5_com_2_reg_1+pos_alleles_6_com_2_reg_1+pos_alleles_7_com_2_reg_1+pos_alleles_8_com_2_reg_1+pos_alleles_9_com_2_reg_1+pos_alleles_10_com_2_reg_1+pos_alleles_11_com_2_reg_1+pos_alleles_12_com_2_reg_1+pos_alleles_none_com_2_reg_1

        if pos_rel_ans in input_mis_pos:
            x = input_mis_pos.replace(pos_rel_ans,'')
            for i in range(len(x)):
                for l in combinations(x, i+1):
                    s = ', '.join(l)
                    com_3 = s
                    global pos_alleles_2_com_3_reg_1, pos_alleles_3_com_3_reg_1, pos_alleles_4_com_3_reg_1, pos_alleles_5_com_3_reg_1, pos_alleles_6_com_3_reg_1, pos_alleles_7_com_3_reg_1, pos_alleles_8_com_3_reg_1, pos_alleles_9_com_3_reg_1, pos_alleles_10_com_3_reg_1, pos_alleles_11_com_3_reg_1, pos_alleles_12_com_3_reg_1, pos_alleles_none_com_3_reg_1
                    if com_3 in sta_allele_2_list:
                        str2 = com_3+' '+'\t-->'+' '+all_2
                        pos_alleles_2_com_3_reg_1 = all_2+' '
                        all_pos_all_2 = 'b'
                        print('\t\t'+str2)
                    elif com_3 in sta_allele_3_list:
                        str3 = com_3+' '+'\t-->'+' '+all_3
                        pos_alleles_3_com_3_reg_1 = all_3+' '
                        all_pos_all_3 = 'c'
                        print('\t\t'+str3)
                    elif com_3 in sta_allele_4_list:
                        str4 = com_3+' '+'\t-->'+' '+all_4
                        pos_alleles_4_com_3_reg_1 = all_4+' '
                        all_pos_all_4 = 'd'
                        print('\t\t'+str4)
                    elif com_3 in sta_allele_5_list:
                        str5 = com_3+' '+'\t-->'+' '+all_5
                        pos_alleles_5_com_3_reg_1 = all_5+' '
                        all_pos_all_5 = 'e'
                        print('\t\t'+str5)
                    elif com_3 in sta_allele_6_list:
                        str6 = com_3+' '+'\t-->'+' '+all_6
                        pos_alleles_6_com_3_reg_1 = all_6+' '
                        all_pos_all_6 = 'f'
                        print('\t\t'+str6)
                    elif com_3 in sta_allele_7_list:
                        str7 = com_3+' '+'\t-->'+' '+all_7
                        pos_alleles_7_com_3_reg_1 = all_7+' '
                        all_pos_all_7 = 'g'
                        print('\t\t'+str7)
                    elif com_3 in sta_allele_8_list:
                        str8 = com_3+' '+'\t-->'+' '+all_8
                        pos_alleles_8_com_3_reg_1 = all_8+' '
                        all_pos_all_8 = 'h'
                        print('\t\t'+str8)
                    elif com_3 in sta_allele_9_list:
                        str9 = com_3+' '+'\t-->'+' '+all_9
                        pos_alleles_9_com_3_reg_1 = all_9+' '
                        all_pos_all_9 = 'i'
                        print('\t\t'+str9)
                    elif com_3 in sta_allele_10_list:
                        str10 = com_3+' '+'\t-->'+' '+all_10
                        pos_alleles_10_com_3_reg_1 = all_10+' '
                        all_pos_all_10 = 'j'
                        print('\t\t'+str10)
                    elif com_3 in sta_allele_11_list:
                        str11 = com_3+' '+'\t-->'+' '+all_11
                        pos_alleles_11_com_3_reg_1 = all_11+' '
                        all_pos_all_11 = 'k'
                        print('\t\t'+str11)
                    elif com_3 in sta_allele_12_list:
                        str12 = com_3+' '+'\t-->'+' '+all_12
                        pos_alleles_12_com_3_reg_1 = all_12+' '
                        all_pos_all_12 = 'l'
                        print('\t\t'+str12)
                    else:
                        str_none = com_3+' '+'\t-->'+' '+'None'
                        pos_alleles_none_com_3_reg_1 = ''
                        print('\t\t'+str_none)
        else:
            #x = 'None'
            #print('\t\t'+x)
            None
        print("")
        all_pos_alleles_com_3_reg_1 = pos_alleles_2_com_3_reg_1+pos_alleles_3_com_3_reg_1+pos_alleles_4_com_3_reg_1+pos_alleles_5_com_3_reg_1+pos_alleles_6_com_3_reg_1+pos_alleles_7_com_3_reg_1+pos_alleles_8_com_3_reg_1+pos_alleles_9_com_3_reg_1+pos_alleles_10_com_3_reg_1+pos_alleles_11_com_3_reg_1+pos_alleles_12_com_3_reg_1+pos_alleles_none_com_3_reg_1
        
        

        print("\t9/1). If all position relate answer is missing, add *1.")
        if pos_rel_ans in input_mis_pos:
            str_all_pos_all_reg_1 = all_pos_all_1+all_pos_all_2+all_pos_all_3+all_pos_all_4+all_pos_all_5+all_pos_all_6+all_pos_all_7+all_pos_all_8+all_pos_all_9+all_pos_all_10+all_pos_all_11+all_pos_all_12
            lis_all_pos_all_reg_1 = []
            for i in str_all_pos_all_reg_1:
                if i == 'a':
                    i = '*1'
                elif i == 'b':
                    i = '*2'
                elif i == 'c':
                    i = '*3'
                elif i == 'd':
                    i = '*4'
                elif i == 'e':
                    i = '*5'
                elif i == 'f':
                    i = '*6'
                elif i == 'g':
                    i = '*7'
                elif i == 'h':
                    i = '*8'
                elif i == 'i':
                    i = '*9'
                elif i == 'j':
                    i = '*10'
                elif i == 'k':
                    i = '*11'
                elif i == 'l':
                    i = '*12'
                lis_all_pos_all_reg_1.append(i)

            txt1 = '\t\tAll position relate answer = '+pos_rel_ans+', '+'that in list missing position '+', '.join(input_mis_pos_list)
            txt2  = '\t\tSo, all possible alleles for answer of missing in region 1. = '
            print(txt1)
            print("")
            print(txt2+', '.join(lis_all_pos_all_reg_1))
        else:
            str_all_pos_all_reg_1 = all_pos_all_2+all_pos_all_3+all_pos_all_4+all_pos_all_5+all_pos_all_6+all_pos_all_7+all_pos_all_8+all_pos_all_9+all_pos_all_10+all_pos_all_11+all_pos_all_12
            lis_all_pos_all_reg_1 = []
            for i in str_all_pos_all_reg_1:
                if i == 'a':
                    i = '*1'
                elif i == 'b':
                    i = '*2'
                elif i == 'c':
                    i = '*3'
                elif i == 'd':
                    i = '*4'
                elif i == 'e':
                    i = '*5'
                elif i == 'f':
                    i = '*6'
                elif i == 'g':
                    i = '*7'
                elif i == 'h':
                    i = '*8'
                elif i == 'i':
                    i = '*9'
                elif i == 'j':
                    i = '*10'
                elif i == 'k':
                    i = '*11'
                elif i == 'l':
                    i = '*12'
                lis_all_pos_all_reg_1.append(i)

            txt1 = '\t\tAll position relate answer = '+pos_rel_ans+', '+'that not in list missing position '+', '.join(input_mis_pos_list)
            txt2 = '\t\tSo, all possible alleles for answer of missing in region 1. = '
            print(txt1)
            print("")
            print(txt2+', '.join(lis_all_pos_all_reg_1))
        print("")

        print("\t10). Repeat from 2). to 9). for region 2.")
        print("")

        print("\t2/2). Find positions relate answers.")
        if input_ans_dip_reg_2 == '*2':
            pos_rel_ans = '2'
        elif input_ans_dip_reg_2 == '*3':
            pos_rel_ans = '3'
        print('\t\t'+pos_rel_ans)
        print("")

        print("\t3/2). Insert positions relate answers into all combination in 1).")
        x = input_mis_pos
        for i in range(len(x)):
            for l in combinations(x, i+1):
                s = '\t\t'+pos_rel_ans+', '+(', '.join(l))
                print(s)
        print("")

        print("\t4/2). From 1). select combination include positions relate answer. And remove positions relate answer.")
        if pos_rel_ans in input_mis_pos:
            x = input_mis_pos.replace(pos_rel_ans,'')
            for i in range(len(x)):
                for l in combinations(x, i+1):
                    print('\t\t'+', '.join(l))
        else:
            #x = 'None'
            #print('\t\t'+x)
            None
        print("")

        print("\t5/2). Merge all combination from 2/2). 3/2). 4/2).")
        if input_ans_dip_reg_2 == '*2':
            pos_rel_ans = '2'
        elif input_ans_dip_reg_2 == '*3':
            pos_rel_ans = '3'
        print('\t\t'+pos_rel_ans)
        print("")

        x = input_mis_pos
        for i in range(len(x)):
            for l in combinations(x, i+1):
                s = '\t\t'+pos_rel_ans+', '+(', '.join(l))
                print(s)
        print("")

        if pos_rel_ans in input_mis_pos:
            x = input_mis_pos.replace(pos_rel_ans,'')
            for i in range(len(x)):
                for l in combinations(x, i+1):
                    print('\t\t'+', '.join(l))
        else:
            #x = 'None'
            #print('\t\t'+x)
            None
        print("")

        print("\t6/2). Each combination from 5/2). remove duplicate positions.")
        if input_ans_dip_reg_2 == '*2':
            pos_rel_ans = '2'
        elif input_ans_dip_reg_2 == '*3':
            pos_rel_ans = '3'
        print('\t\t'+pos_rel_ans)
        print("")

        x = input_mis_pos
        for i in range(len(x)):
            for l in combinations(x, i+1):
                if pos_rel_ans in l:
                    s = '\t\t'+(', '.join(l))
                else:
                    s = '\t\t'+pos_rel_ans+', '+(', '.join(l))
                print(s)
        print("")

        if pos_rel_ans in input_mis_pos:
            x = input_mis_pos.replace(pos_rel_ans,'')
            for i in range(len(x)):
                for l in combinations(x, i+1):
                    print('\t\t'+', '.join(l))
        else:
            #x = 'None'
            #print('\t\t'+x)
            None
        print("")

        print("\t7/2). Remove duplicate combination from 6/2).")
        if input_ans_dip_reg_2 == '*2':
            pos_rel_ans = '2'
        elif input_ans_dip_reg_2 == '*3':
            pos_rel_ans = '3'
        print('\t\t'+pos_rel_ans)

        if pos_rel_ans in input_mis_pos:
            x = input_mis_pos.replace(pos_rel_ans,'')
            for i in range(len(x)):
                for l in combinations(x, i+1):
                    s = '\t\t'+pos_rel_ans+', '+(', '.join(l))
                    print(s)
        else:
            x = input_mis_pos.replace(pos_rel_ans,'')
            for i in range(len(x)):
                for l in combinations(x, i+1):
                    s = '\t\t'+pos_rel_ans+', '+(', '.join(l))
                    print(s)

        if pos_rel_ans in input_mis_pos:
            x = input_mis_pos.replace(pos_rel_ans,'')
            for i in range(len(x)):
                for l in combinations(x, i+1):
                    print('\t\t'+', '.join(l))
        else:
            #x = 'None'
            #print('\t\t'+x)
            None
        print("")

        print("\t8/2). Find star allele relate positions from 7/2).")

        sta_allele_2 = '2'
        sta_allele_3 = '3'
        sta_allele_4 = '4'
        sta_allele_5 = '5'
        sta_allele_6 = '6'
        sta_allele_7 = '1, 3'
        sta_allele_8 = '2, 6'
        sta_allele_9 = '3, 5'
        sta_allele_10 = '2, 4'
        sta_allele_11 = '1, 4'
        sta_allele_12 = '5, 6'

        sta_allele_2_list = ['2']
        sta_allele_3_list = ['3']
        sta_allele_4_list = ['4']
        sta_allele_5_list = ['5']
        sta_allele_6_list = ['6']
        sta_allele_7_list = ['1, 3', '3, 1']
        sta_allele_8_list = ['2, 6', '6, 2']
        sta_allele_9_list = ['3, 5', '5, 3']
        sta_allele_10_list = ['2, 4', '4, 2']
        sta_allele_11_list = ['1, 4', '4, 1']
        sta_allele_12_list = ['5, 6', '5, 5']

        all_pos_all_1 = 'a'
        all_pos_all_2 = ''
        all_pos_all_3 = ''
        all_pos_all_4 = ''
        all_pos_all_5 = ''
        all_pos_all_6 = ''
        all_pos_all_7 = ''
        all_pos_all_8 = ''
        all_pos_all_9 = ''
        all_pos_all_10 = ''
        all_pos_all_11 = ''
        all_pos_all_12 = ''

        if input_ans_dip_reg_2 == '*2':
            pos_rel_ans = '2'
        elif input_ans_dip_reg_2 == '*3':
            pos_rel_ans = '3'
        com_1 = pos_rel_ans

        global pos_alleles_com_1_reg_2, pos_alleles_none_com_1_reg_2
        if com_1 in sta_allele_2_list:
            str2 = sta_allele_2+' '+'\t-->'+' '+all_2
            pos_alleles_com_1_reg_2 = all_2+' '
            all_pos_all_2 = 'b'
            print('\t\t'+str2)
        elif com_1 in sta_allele_3_list:
            str3 = sta_allele_3+' '+'\t-->'+' '+all_3
            pos_alleles_com_1_reg_2 = all_3+' '
            all_pos_all_3 = 'c'
            print('\t\t'+str3)
        elif com_1 in sta_allele_4_list:
            str4 = sta_allele_4+' '+'\t-->'+' '+all_4
            pos_alleles_com_1_reg_2 = all_4+' '
            all_pos_all_4 = 'd'
            print('\t\t'+str4)
        elif com_1 in sta_allele_5_list:
            str5 = sta_allele5+' '+'\t-->'+' '+all_5
            pos_alleles_com_1_reg_2 = all_5+' '
            all_pos_all_5 = 'e'
            print('\t\t'+str5)
        elif com_1 in sta_allele_6_list:
            str6 = sta_allele_6+' '+'\t-->'+' '+all_6
            pos_alleles_com_1_reg_2 = all_6+' '
            all_pos_all_6 = 'f'
            print('\t\t'+str6)
        elif com_1 in sta_allele_7_list:
            str7 = sta_allele_7+' '+'\t-->'+' '+all_7
            pos_alleles_com_1_reg_2 = all_7+' '
            all_pos_all_7 = 'g'
            print('\t\t'+str7)
        elif com_1 in sta_allele_8_list:
            str8 = sta_allele_8+' '+'\t-->'+' '+all_8
            pos_alleles_com_1_reg_2 = all_8+' '
            all_pos_all_8 = 'h'
            print('\t\t'+str8)
        elif com_1 in sta_allele_9_list:
            str9 = sta_allele_9+' '+'\t-->'+' '+all_9
            pos_alleles_com_1_reg_2 = all_9+' '
            all_pos_all_9 = 'i'
            print('\t\t'+str9)
        elif com_1 in sta_allele_10_list:
            str10 = sta_allele_10+' '+'\t-->'+' '+all_10
            pos_alleles_com_1_reg_2 = all_10+' '
            all_pos_all_10 = 'j'
            print('\t\t'+str10)
        elif com_1 in sta_allele_11_list:
            str11 = sta_allele_11+' '+'\t-->'+' '+all_11
            pos_alleles_com_1_reg_2 = all_11+' '
            all_pos_all_11 = 'k'
            print('\t\t'+str11)
        elif com_1 in sta_allele_12_list:
            str12 = sta_allele_12+' '+'\t-->'+' '+all_12
            pos_alleles_com_1_reg_2 = all_12+' '
            all_pos_all_12 = 'l'
            print('\t\t'+str12)
        else:
            str_none = com_1+' '+'\t-->'+' '+'None'
            pos_alleles_none_com_1_reg_2 = ''
            print('\t\t'+str_none)
        all_pos_alleles_com_1_reg_2 = pos_alleles_com_1_reg_2+pos_alleles_none_com_1_reg_2

        if pos_rel_ans in input_mis_pos:
            x = input_mis_pos.replace(pos_rel_ans,'')
            for i in range(len(x)):
                for l in combinations(x, i+1):
                    s = pos_rel_ans+', '+(', '.join(l))
                    com_2 = s
                    global pos_alleles_2_com_2_reg_2, pos_alleles_3_com_2_reg_2, pos_alleles_4_com_2_reg_2, pos_alleles_5_com_2_reg_2, pos_alleles_6_com_2_reg_2, pos_alleles_7_com_2_reg_2, pos_alleles_8_com_2_reg_2, pos_alleles_9_com_2_reg_2, pos_alleles_10_com_2_reg_2, pos_alleles_11_com_2_reg_2, pos_alleles_12_com_2_reg_2, pos_alleles_none_com_2_reg_2
                    if com_2 in sta_allele_2_list:
                        str2 = com_2+' '+'\t-->'+' '+all_2
                        pos_alleles_2_com_2_reg_2 = all_2+' '
                        all_pos_all_2 = 'b'
                        print('\t\t'+str2)
                    elif com_2 in sta_allele_3_list:
                        str3 = com_2+' '+'\t-->'+' '+all_3
                        pos_alleles_3_com_2_reg_2 = all_3+' '
                        all_pos_all_3 = 'c'
                        print('\t\t'+str3)
                    elif com_2 in sta_allele_4_list:
                        str4 = com_2+' '+'\t-->'+' '+all_4
                        pos_alleles_4_com_2_reg_2 = all_4+' '
                        all_pos_all_4 = 'd'
                        print('\t\t'+str4)
                    elif com_2 in sta_allele_5_list:
                        str5 = com_2+' '+'\t-->'+' '+all_5
                        pos_alleles_5_com_2_reg_2 = all_5+' '
                        all_pos_all_5 = 'e'
                        print('\t\t'+str5)
                    elif com_2 in sta_allele_6_list:
                        str6 = com_2+' '+'\t-->'+' '+all_6
                        pos_alleles_6_com_2_reg_2 = all_6+' '
                        all_pos_all_6 = 'f'
                        print('\t\t'+str6)
                    elif com_2 in sta_allele_7_list:
                        str7 = com_2+' '+'\t-->'+' '+all_7
                        pos_alleles_7_com_2_reg_2 = all_7+' '
                        all_pos_all_7 = 'g'
                        print('\t\t'+str7)
                    elif com_2 in sta_allele_8_list:
                        str8 = com_2+' '+'\t-->'+' '+all_8
                        pos_alleles_8_com_2_reg_2 = all_8+' '
                        all_pos_all_8 = 'h'
                        print('\t\t'+str8)
                    elif com_2 in sta_allele_9_list:
                        str9 = com_2+' '+'\t-->'+' '+all_9
                        pos_alleles_9_com_2_reg_2 = all_9+' '
                        all_pos_all_9 = 'i'
                        print('\t\t'+str9)
                    elif com_2 in sta_allele_10_list:
                        str10 = com_2+' '+'\t-->'+' '+all_10
                        pos_alleles_10_com_2_reg_2 = all_10+' '
                        all_pos_all_10 = 'j'
                        print('\t\t'+str10)
                    elif com_2 in sta_allele_11_list:
                        str11 = com_2+' '+'\t-->'+' '+all_11
                        pos_alleles_11_com_2_reg_2 = all_11+' '
                        all_pos_all_11 = 'k'
                        print('\t\t'+str11)
                    elif com_2 in sta_allele_12_list:
                        str12 = com_2+' '+'\t-->'+' '+all_12
                        pos_alleles_12_com_2_reg_2 = all_12+' '
                        all_pos_all_12 = 'l'
                        print('\t\t'+str12)
                    else:
                        str_none = com_2+' '+'\t-->'+' '+'None'
                        pos_alleles_none_com_2_reg_2 = ''
                        print('\t\t'+str_none)
        else:
            x = input_mis_pos.replace(pos_rel_ans,'')
            for i in range(len(x)):
                for l in combinations(x, i+1):
                    s = pos_rel_ans+', '+(', '.join(l))
                    com_2 = s
                    if com_2 in sta_allele_2_list:
                        str2 = com_2+' '+'\t-->'+' '+all_2
                        pos_alleles_2_com_2_reg_2 = all_2+' '
                        all_pos_all_2 = 'b'
                        print('\t\t'+str2)
                    elif com_2 in sta_allele_3_list:
                        str3 = com_2+' '+'\t-->'+' '+all_3
                        pos_alleles_3_com_2_reg_2 = all_3+' '
                        all_pos_all_3 = 'c'
                        print('\t\t'+str3)
                    elif com_2 in sta_allele_4_list:
                        str4 = com_2+' '+'\t-->'+' '+all_4
                        pos_alleles_4_com_2_reg_2 = all_4+' '
                        all_pos_all_4 = 'd'
                        print('\t\t'+str4)
                    elif com_2 in sta_allele_5_list:
                        str5 = com_2+' '+'\t-->'+' '+all_5
                        pos_alleles_5_com_2_reg_2 = all_5+' '
                        all_pos_all_5 = 'e'
                        print('\t\t'+str5)
                    elif com_2 in sta_allele_6_list:
                        str6 = com_2+' '+'\t-->'+' '+all_6
                        pos_alleles_6_com_2_reg_2 = all_6+' '
                        all_pos_all_6 = 'f'
                        print('\t\t'+str6)
                    elif com_2 in sta_allele_7_list:
                        str7 = com_2+' '+'\t-->'+' '+all_7
                        pos_alleles_7_com_2_reg_2 = all_7+' '
                        all_pos_all_7 = 'g'
                        print('\t\t'+str7)
                    elif com_2 in sta_allele_8_list:
                        str8 = com_2+' '+'\t-->'+' '+all_8
                        pos_alleles_8_com_2_reg_2 = all_8+' '
                        all_pos_all_8 = 'h'
                        print('\t\t'+str8)
                    elif com_2 in sta_allele_9_list:
                        str9 = com_2+' '+'\t-->'+' '+all_9
                        pos_alleles_9_com_2_reg_2 = all_9+' '
                        all_pos_all_9 = 'i'
                        print('\t\t'+str9)
                    elif com_2 in sta_allele_10_list:
                        str10 = com_2+' '+'\t-->'+' '+all_10
                        pos_alleles_10_com_2_reg_2 = all_10+' '
                        all_pos_all_10 = 'j'
                        print('\t\t'+str10)
                    elif com_2 in sta_allele_11_list:
                        str11 = com_2+' '+'\t-->'+' '+all_11
                        pos_alleles_11_com_2_reg_2 = all_11+' '
                        all_pos_all_11 = 'k'
                        print('\t\t'+str11)
                    elif com_2 in sta_allele_12_list:
                        str12 = com_2+' '+'\t-->'+' '+all_12
                        pos_alleles_12_com_2_reg_2 = all_12+' '
                        all_pos_all_12 = 'l'
                        print('\t\t'+str12)
                    else:
                        str_none = com_2+' '+'\t-->'+' '+'None'
                        pos_alleles_none_com_2_reg_2 = ''
                        print('\t\t'+str_none)
        all_pos_alleles_com_2_reg_2 = pos_alleles_2_com_2_reg_2+pos_alleles_3_com_2_reg_2+pos_alleles_4_com_2_reg_2+pos_alleles_5_com_2_reg_2+pos_alleles_6_com_2_reg_2+pos_alleles_7_com_2_reg_2+pos_alleles_8_com_2_reg_2+pos_alleles_9_com_2_reg_2+pos_alleles_10_com_2_reg_2+pos_alleles_11_com_2_reg_2+pos_alleles_12_com_2_reg_2+pos_alleles_none_com_2_reg_2

        if pos_rel_ans in input_mis_pos:
            x = input_mis_pos.replace(pos_rel_ans,'')
            for i in range(len(x)):
                for l in combinations(x, i+1):
                    s = ', '.join(l)
                    com_3 = s
                    global pos_alleles_2_com_3_reg_2, pos_alleles_3_com_3_reg_2, pos_alleles_4_com_3_reg_2, pos_alleles_5_com_3_reg_2, pos_alleles_6_com_3_reg_2, pos_alleles_7_com_3_reg_2, pos_alleles_8_com_3_reg_2, pos_alleles_9_com_3_reg_2, pos_alleles_10_com_3_reg_2, pos_alleles_11_com_3_reg_2, pos_alleles_12_com_3_reg_2, pos_alleles_none_com_3_reg_2
                    if com_3 in sta_allele_2_list:
                        str2 = com_3+' '+'\t-->'+' '+all_2
                        pos_alleles_2_com_3_reg_2 = all_2+' '
                        all_pos_all_2 = 'b'
                        print('\t\t'+str2)
                    elif com_3 in sta_allele_3_list:
                        str3 = com_3+' '+'\t-->'+' '+all_3
                        pos_alleles_3_com_3_reg_2 = all_3+' '
                        all_pos_all_3 = 'c'
                        print('\t\t'+str3)
                    elif com_3 in sta_allele_4_list:
                        str4 = com_3+' '+'\t-->'+' '+all_4
                        pos_alleles_4_com_3_reg_2 = all_4+' '
                        all_pos_all_4 = 'd'
                        print('\t\t'+str4)
                    elif com_3 in sta_allele_5_list:
                        str5 = com_3+' '+'\t-->'+' '+all_5
                        pos_alleles_5_com_3_reg_2 = all_5+' '
                        all_pos_all_5 = 'e'
                        print('\t\t'+str5)
                    elif com_3 in sta_allele_6_list:
                        str6 = com_3+' '+'\t-->'+' '+all_6
                        pos_alleles_6_com_3_reg_2 = all_6+' '
                        all_pos_all_6 = 'f'
                        print('\t\t'+str6)
                    elif com_3 in sta_allele_7_list:
                        str7 = com_3+' '+'\t-->'+' '+all_7
                        pos_alleles_7_com_3_reg_2 = all_7+' '
                        all_pos_all_7 = 'g'
                        print('\t\t'+str7)
                    elif com_3 in sta_allele_8_list:
                        str8 = com_3+' '+'\t-->'+' '+all_8
                        pos_alleles_8_com_3_reg_2 = all_8+' '
                        all_pos_all_8 = 'h'
                        print('\t\t'+str8)
                    elif com_3 in sta_allele_9_list:
                        str9 = com_3+' '+'\t-->'+' '+all_9
                        pos_alleles_9_com_3_reg_2 = all_9+' '
                        all_pos_all_9 = 'i'
                        print('\t\t'+str9)
                    elif com_3 in sta_allele_10_list:
                        str10 = com_3+' '+'\t-->'+' '+all_10
                        pos_alleles_10_com_3_reg_2 = all_10+' '
                        all_pos_all_10 = 'j'
                        print('\t\t'+str10)
                    elif com_3 in sta_allele_11_list:
                        str11 = com_3+' '+'\t-->'+' '+all_11
                        pos_alleles_11_com_3_reg_2 = all_11+' '
                        all_pos_all_11 = 'k'
                        print('\t\t'+str11)
                    elif com_3 in sta_allele_12_list:
                        str12 = com_3+' '+'\t-->'+' '+all_12
                        pos_alleles_12_com_3_reg_2 = all_12+' '
                        all_pos_all_12 = 'l'
                        print('\t\t'+str12)
                    else:
                        str_none = com_3+' '+'\t-->'+' '+'None'
                        pos_alleles_none_com_3_reg_2 = ''
                        print('\t\t'+str_none)
        else:
            #x = 'None'
            #print('\t\t'+x)
            None
        print("")
        all_pos_alleles_com_3_reg_2 = pos_alleles_2_com_3_reg_2+pos_alleles_3_com_3_reg_2+pos_alleles_4_com_3_reg_2+pos_alleles_5_com_3_reg_2+pos_alleles_6_com_3_reg_2+pos_alleles_7_com_3_reg_2+pos_alleles_8_com_3_reg_2+pos_alleles_9_com_3_reg_2+pos_alleles_10_com_3_reg_2+pos_alleles_11_com_3_reg_2+pos_alleles_12_com_3_reg_2+pos_alleles_none_com_3_reg_2
        
        str_all_pos_all  = all_pos_all_2+all_pos_all_3+all_pos_all_4+all_pos_all_5+all_pos_all_6+all_pos_all_7+all_pos_all_8+all_pos_all_9+all_pos_all_10+all_pos_all_11+all_pos_all_12
        lis_all_pos_all = []
        for i in str_all_pos_all:
            if i == 'a':
                i = '*1'
            elif i == 'b':
                i = '*2'
            elif i == 'c':
                i = '*3'
            elif i == 'd':
                i = '*4'
            elif i == 'e':
                i = '*5'
            elif i == 'f':
                i = '*6'
            elif i == 'g':
                i = '*7'
            elif i == 'h':
                i = '*8'
            elif i == 'i':
                i = '*9'
            elif i == 'j':
                i = '*10'
            elif i == 'k':
                i = '*11'
            elif i == 'l':
                i = '*12'
            lis_all_pos_all.append(i)
        txt1  = '\t\tAll possible alleles for answer of missing in region 1. = '
        print(txt1+', '.join(lis_all_pos_all))
        print("")

        print("\t9/2). If all position relate answer is missing, add *1.")
        if pos_rel_ans in input_mis_pos:
            str_all_pos_all_reg_2 = all_pos_all_1+all_pos_all_2+all_pos_all_3+all_pos_all_4+all_pos_all_5+all_pos_all_6+all_pos_all_7+all_pos_all_8+all_pos_all_9+all_pos_all_10+all_pos_all_11+all_pos_all_12
            lis_all_pos_all_reg_2 = []
            for i in str_all_pos_all_reg_2:
                if i == 'a':
                    i = '*1'
                elif i == 'b':
                    i = '*2'
                elif i == 'c':
                    i = '*3'
                elif i == 'd':
                    i = '*4'
                elif i == 'e':
                    i = '*5'
                elif i == 'f':
                    i = '*6'
                elif i == 'g':
                    i = '*7'
                elif i == 'h':
                    i = '*8'
                elif i == 'i':
                    i = '*9'
                elif i == 'j':
                    i = '*10'
                elif i == 'k':
                    i = '*11'
                elif i == 'l':
                    i = '*12'
                lis_all_pos_all_reg_2.append(i)

            txt1 = '\t\tAll position relate answer = '+pos_rel_ans+', '+'that in list missing position '+', '.join(input_mis_pos_list)
            txt2  = '\t\tSo, all possible alleles for answer of missing in region 1. = '
            print(txt1)
            print("")
            print(txt2+', '.join(lis_all_pos_all_reg_2))
        else:
            str_all_pos_all_reg_2 = all_pos_all_2+all_pos_all_3+all_pos_all_4+all_pos_all_5+all_pos_all_6+all_pos_all_7+all_pos_all_8+all_pos_all_9+all_pos_all_10+all_pos_all_11+all_pos_all_12
            lis_all_pos_all_reg_2 = []
            for i in str_all_pos_all_reg_2:
                if i == 'a':
                    i = '*1'
                elif i == 'b':
                    i = '*2'
                elif i == 'c':
                    i = '*3'
                elif i == 'd':
                    i = '*4'
                elif i == 'e':
                    i = '*5'
                elif i == 'f':
                    i = '*6'
                elif i == 'g':
                    i = '*7'
                elif i == 'h':
                    i = '*8'
                elif i == 'i':
                    i = '*9'
                elif i == 'j':
                    i = '*10'
                elif i == 'k':
                    i = '*11'
                elif i == 'l':
                    i = '*12'
                lis_all_pos_all_reg_2.append(i)

            txt1 = '\t\tAll position relate answer = '+pos_rel_ans+', '+'that not in list missing position '+', '.join(input_mis_pos_list)
            txt2 = '\t\tSo, all possible alleles for answer of missing in region 1. = '
            print(txt1)
            print("")
            print(txt2+', '.join(lis_all_pos_all_reg_2))
        print("")

        print("\t10). Already repeat from 2). to 9). for region 1 and region 2.")
        print("")


        print("\t11). From 9/1). and 9/2). find all diplotypes combination of allele in each region.")
        print('\t\t'+reg_2+' '+', '.join(lis_all_pos_all_reg_1))
        print('\t\t'+reg_3+' '+', '.join(lis_all_pos_all_reg_2))
        print("")

        lis_all_pos_all_com = []
        for i in str_all_pos_all_reg_1:
            if i == 'a':
                i = '*1'
            elif i == 'b':
                i = '*2'
            elif i == 'c':
                i = '*3'
            elif i == 'd':
                i = '*4'
            elif i == 'e':
                i = '*5'
            elif i == 'f':
                i = '*6'
            elif i == 'g':
                i = '*7'
            elif i == 'h':
                i = '*8'
            elif i == 'i':
                i = '*9'
            elif i == 'j':
                i = '*10'
            elif i == 'k':
                i = '*11'
            elif i == 'l':
                i = '*12'
            for j in str_all_pos_all_reg_2:
                if j == 'a':
                    j = '*1'
                elif j == 'b':
                    j = '*2'
                elif j == 'c':
                    j = '*3'
                elif j == 'd':
                    j = '*4'
                elif j == 'e':
                    j = '*5'
                elif j == 'f':
                    j = '*6'
                elif j == 'g':
                    j = '*7'
                elif j == 'h':
                    j = '*8'
                elif j == 'i':
                    j = '*9'
                elif j == 'j':
                    j = '*10'
                elif j == 'k':
                    j = '*11'
                elif j == 'l':
                    j = '*12'
                lis_all_pos_all_com.append((i+'/'+ j))

        txt1 = '\t\tAll possible diplotypes = '
        str_all_pos_all_com = ', '.join(str(v) for v in lis_all_pos_all_com)
        print(txt1+' '+str_all_pos_all_com)

        print("")

        # Module complete
        print("Module 5: Find answers complete")
        print("")
        break

# Module Main ##################################################
def main_phased_genotypes_missing_positions():
    while True:
        print("")
        print("\tMenu: Case: Phased genotypes with missing positions")
        print("")
        print("\t\tModule 1. Create table\t\t\t\t\t: Button[1]")
        print("\t\tModule 2. Load table\t\t\t\t\t: Button[2]")
        print("\t\tModule 3. Example simulate VCF file\t\t\t: Button[3]")
        print("\t\tModule 4. Example test cases\t\t\t\t: Button[4]")
        print("\t\tModule 5. Find answers\t\t\t\t\t: Button[5]")
        print("\t\tModule 6. Main Menu \t\t\t\t\t: Button[6]")
        print("\t\tModule 7. Exit. \t\t\t\t\t: Button[7]")
        print("")
        try:
            button = int(input("\t\t\tPlease push your Button. \t\t\t: "))
        except ValueError:
            print("")
            print("\t\t\t\t\t\t\t\t\t: No Module.")
            print("")
            continue
        else:
            if button == 1:
                #print("")
                #create_table()
                print("")
                print("\t\t\t\t\t\t\t\t\t: Coming soon!")
                print("")
            elif button == 2:
                print("")
                load_table()
            elif button == 3:
                print("")
                if p1a1 != '':
                    example_simulate_vcf_file()
                else:
                    print("")
                    print("\t\t\t\t\t\t\t\t\t: Please load table in Module 2.")
                    print("")
            elif button == 4:
                print("")
                if p1a1 != '':
                    example_test_cases()
                else:
                    print("")
                    print("\t\t\t\t\t\t\t\t\t: Please load table in Module 2.")
                    print("")
            elif button == 5:
                print("")
                if p1a1 != '':
                    find_answers()
                else:
                    print("")
                    print("\t\t\t\t\t\t\t\t\t: Please load table in Module 2.")
                    print("")
            elif button == 6:
                print("")
                import main
                main.main()
            elif button == 7:
                print("")
                bar = progressbar.ProgressBar(max_value=progressbar.UnknownLength)
                for i in range(100):
                    time.sleep(0.1)
                    bar.update(i)
                print(" Exit.")
                print("")
                sys.exit(0)
            else:
                print("")
                print("\t\t\t\t\t\t\t\t\t: No Module.")
                print("")

main_phased_genotypes_missing_positions()
