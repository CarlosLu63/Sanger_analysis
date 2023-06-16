from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.Seq import Seq
import sys
import os

path = "/home/carloslu/Sanger/"
sample_path = "/home/carloslu/Sanger/Sample_file/"
temp_file = "/home/carloslu/Sanger/temp/"
output_path = "/home/carloslu/Sanger/Output_file/"
os.chdir(path)

#####Create database#####
#os.system("/opt/ncbi-blast-2.13.0+/bin/makeblastdb -in GJB2_db2.fasta -dbtype nucl")

#####Parse sample abi file#####
info_dict = {}
for record in SeqIO.parse(sample_path + "E1_GJB2-53_20220412_171323.ab1", "abi"):
    info_dict["seq"] = record.seq
    info_dict["id"] = record.id
    
file = open("/home/carloslu/Sanger/temp/" + info_dict["id"] + ".txt","w")
reads = ">" + info_dict["id"] + "\n" + info_dict["seq"].lower()
file.writelines(reads)
file.close()

#####RUN blastn#####
os.system("/opt/ncbi-blast-2.13.0+/bin/blastn -task blastn -outfmt 5 -db " + 
          path + "GJB2_db2.fasta -query " + 
          temp_file + info_dict["id"] + ".txt " + 
          "-out " + output_path + info_dict["id"] + ".xml")

#####Read xml file#####
result_handle = open(output_path + info_dict["id"] + ".xml")
blast_record = NCBIXML.read(result_handle)

#設定比對率
E_VALUE_THRESH = 0.0000004
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

my_seq = Seq('agcaaaccgcccagagtagaagatggattggggcacgctgcagacgatcctggggggtgtgaacaaacactccaccagcattggaaagatctggctcaccgtcctcttcatttttcgcattatgatcctcgttgtggctgcaaaggaggtgtggggagatgagcaggccgactttgtctgcaacaccctgcagccaggctgcaagaacgtgtgctacgatcactacttccccatctcccacatccggctatgggccctgcagctgatcttcgtgtccacgccagcgctcctagtggccatgcacgtggcctaccggagacatgagaagaagaggaagttcatcaagggggagataaagagtgaatttaaggacatcgaggagatcaaaacccagaaggtccgcatcgaaggctccctgtggtggacctacacaagcagcatcttcttccgggtcatcttcgaagccgccttcatgtacgtcttctatgtcatgtacgacggcttctccatgcagcggctggtgaagtgcaacgcctggccttgtcccaacactgtggactgctttgtgtcccggcccacggagaagactgtcttcacagtgttcatgattgcagtgtctggaatttgcatcctgctgaatgtcactgaattgtgttatttgctaattagatattgttctgggaagtcaaaaaagccagtttaacgcattgcccagttgttagattaagaaatagacagcatgagagggatgaggcaacccgtgctcagctgtcaaggctcagtcgctagcatttcccaacacaaagattctgaccttaaatgcaaccatttgaaacccctgtaggcctcaggtgaaactccagatgccacaatggagctctgctcccctaaagcctcaaaacaaaggcctaattctatgcctgtcttaattttctttcacttaagttagttccactgagaccccaggctgttaggggttattggtgtaaggtactttcatattttaaacagaggatatcggcatttgtttctttctctgaggacaagagaaaaaagccaggttccacagaggacacagagaaggtttgggtgtcctcctggggttctttttgccaactttccccacgttaaaggtgaacattggttctttcatttgctttggaagttttaatctctaacagtggacaaagttaccagtgccttaaactctgttacactttttggaagtgaaaactttgtagtatgataggttattttgatgtaaagatgttctggataccattatatgttccccctgtttcagaggctcagattgtaatatgtaaatggtatgtcattcgctactatgatttaatttgaaatatggtcttttggttatgaatactttgcagcacagctgagaggctgtctgttgtattcattgtggtcatagcacctaacaacattgtagcctcaatcgagtgagacagactagaagttcctagtgatggcttatgatagcaaatggcctcatgtcaaatatttagatgtaattttgtgtaagaaatacagactggatgtaccaccaactactacctgtaatgacaggcctgtccaacacatctcccttttccatgactgtggtagccagcatcggaaagaacgctgatttaaagaggtcgcttgggaattttattgacacagtaccatttaatggggaggacaaaatggggcaggggagggagaagtttctgtcgttaaaaacagatttggaaagactggactctaaagtctgttgattaaagatgagctttgtctacttcaaaagtttgtttgcttaccccttcagcctccaattttttaagtgaaaatatagctaataacatgtgaaaagaatagaagctaaggtttagataaatattgagcagatctataggaagattgaacctgaatattgccattatgcttgacatggtttccaaaaaatggtactccacatatttcagtgagggtaagtattttcctgttgtcaagaatagcattgtaaaagcattttgtaataataaagaatagctttaatgatatgcttgtaactaaaataattttgtaatgtatcaaatacatttaaaacattaaaatataatctctataataa')
my_seq2 = Seq('tttatgatgtgtttaaagattgggtgaattactcaggtgaacaagctactttttatcagagaacacctaaaaacacgttcaagagggtttgggaactatacatttaatcctatgacaaactaagttggttctgtcttcacctgttttggtgaggttgtgtaagagttggtgtttgctcaggaagagatttaagcatgcttgcttacccagactcagagaagtctccctgttctgtcctagctagtgattcctgtgttgtgtgcattcgtcttttccagagcaaaccgcccagagtagaagATGGATTGGGGCACGCTGCAGACGATCCTGGGGGGTGTGAACAAACACTCCACCAGCATTGGAAAGATCTGGCTCACCGTCCTCTTCATTTTTCGCATTATGATCCTCGTTGTGGCTGCAAAGGAGGTGTGGGGAGATGAGCAGGCCGACTTTGTCTGCAACACCCTGCAGCCAGGCTGCAAGAACGTGTGCTACGATCACTACTTCCCCATCTCCCACATCCGGCTATGGGCCCTGCAGCTGATCTTCGTGTCCACGCCAGCGCTCCTAGTGGCCATGCACGTGGCCTACCGGAGACATGAGAAGAAGAGGAAGTTCATCAAGGGGGAGATAAAGAGTGAATTTAAGGACATCGAGGAGATCAAAACCCAGAAGGTCCGCATCGAAGGCTCCCTGTGGTGGACCTACACAAGCAGCATCTTCTTCCGGGTCATCTTCGAAGCCGCCTTCATGTACGTCTTCTATGTCATGTACGACGGCTTCTCCATGCAGCGGCTGGTGAAGTGCAACGCCTGGCCTTGTCCCAACACTGTGGACTGCTTTGTGTCCCGGCCCACGGAGAAGACTGTCTTCACAGTGTTCATGATTGCAGTGTCTGGAATTTGCATCCTGCTGAATGTCACTGAATTGTGTTATTTGCTAATTAGATATTGTTCTGGGAAGTCAAAAAAGCCAGTTTAAcgcattgcccagttgttagattaagaaatagacagcatgagagggatgaggcaacccgtgctcagctgtcaaggctcagtcgctagcatttcccaacacaaagattctgaccttaaatgcaaccatttgaaacccctgtaggcctcaggtgaaactccagatgccacaatggagctctgctcccctaaagcctcaaaacaaaggcctaattctatgcctgtcttaattttctttcacttaagttagttccactgagaccccaggctgttaggggttattggtgtaaggtactttcatatt')
my_seq2 = my_seq2.lower()
my_seq_rp = my_seq.reverse_complement()
