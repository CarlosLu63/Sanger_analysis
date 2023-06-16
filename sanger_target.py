# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 16:00:33 2023

@author: lucarlos
"""

import Bio
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
import sys
import os
import pandas as pd
import numpy as np
from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline
os.chdir("D:\sanger")

###read ab1 data###
info_dict = {}
for record in SeqIO.parse("A3_GBJ2-17_20220325_172857.ab1", "abi"):
    info_dict["seq"] = record.seq
    info_dict["id"] = record.id
    
file = open(info_dict["id"]+".txt","w")
reads = ">A3_GBJ2-17_20220325_172857\n" + info_dict["seq"].lower()
file.writelines(reads)
file.close()

sequence_data = open(info_dict["id"]+"_20220325_172857.txt").read()
result_handle = NCBIWWW.qblast("blastn", "nt", sequence_data)

with open('D:/sanger/results.xml', 'w') as save_file: 
    blast_results = result_handle.read() 
    save_file.write(blast_results)

result_handle = open("D:/sanger/results.xml")
blast_record = NCBIXML.read(result_handle)

#設定比對率
E_VALUE_THRESH = 0.04
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
            print ('****Alignment****')
            print ('sequence:', alignment.title)
            print ('length:', alignment.length)
            print ('e value:', hsp.expect)
            print (hsp.query[0:75] + '...')
            print (hsp.match[0:75] + '...')
            print (hsp.sbjct[0:75] + '...')

