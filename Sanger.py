from abifpy import Trace
import sys
import os

path = "/home/carloslu/Sanger/"
sample_path = "/home/carloslu/Sanger/Sample_file/"
temp_file = "/home/carloslu/Sanger/temp/"
output_path = "/home/carloslu/Sanger/Output_file/"
ref_path = "/share/NGS_data/Reference/ref_b37_GRCh37/human_g1k_v37.fasta"
scriptDir = "/home/carloslu/Sanger/code/"

os.chdir(path)
#####abifpy#####
record = Trace(sample_path + 'A3_GBJ2-17_20220325_172857.ab1')
sample_name = record.name

###Trimm sequence and quality value with Q20###
record.seq = record.trim(record.seq, cutoff = 0.01)
record.qual = record.trim(record.qual, cutoff = 0.01)
record.export(out_file = temp_file + sample_name + ".fq", fmt = 'fastq')

#####RUN BWA MEM#####
RGline = "'@RG\\tID:"+ sample_name + "\\tPL:illumina\\tLB:" + sample_name + "\\tSM:" + sample_name + "'"
os.system("/opt/bwa/bwa mem -M " +
          "-R " + RGline + " " +
          ref_path +
          " " + temp_file + sample_name + ".fq" + " " +
          "> " + output_path + sample_name + ".bwa.mem.sam")

#####SAM to BAM#####
samoutput = output_path + sample_name + ".bwa.mem.sam"
bamoutput = output_path + sample_name + ".sort.bam"
bamoutput_nameSrt = output_path + sample_name + ".nameSrt.bam"
bamoutput_mate = output_path + sample_name +  ".nameSrt.mate.bam"
bamoutput_mkdp = output_path + sample_name + ".sort.mkdp.bam"
metricsFile = output_path + "MarkDup_metrics.txt"
statsFile = output_path + "MarkDup_stats.txt"

os.system("samtools sort -@ 2 -m 1G -n -O bam -o " + bamoutput_nameSrt  + " " + samoutput)
os.system("samtools fixmate -@ 2 -O bam -m " + bamoutput_nameSrt  + " " + bamoutput_mate)
os.system("samtools sort -@ 2 -m 1G -O bam -o " + bamoutput  + " " + bamoutput_mate)
os.system("samtools markdup -@ 2 -f " + statsFile + " -d 2500 -O bam " + bamoutput  + " " + bamoutput_mkdp)
os.system("samtools index -@ 2 -b " + bamoutput_mkdp)

os.system("unlink " + samoutput)
os.system("unlink " + bamoutput_nameSrt)
os.system("unlink " + bamoutput_mate)
os.system("unlink " + bamoutput)


##### Variant calling 需要用 Winter 的帳號執行!!!
os.system("cp " + bamoutput_mkdp + " -t " + "/share/NGS_Project/GBJ2-17/Mapping/GBJ2-17")
os.system("cp " + bamoutput_mkdp + ".bai -t " + "/share/NGS_Project/GBJ2-17/Mapping/GBJ2-17")

GATK = '/opt/gatk-4.3.0.0/gatk'
os.system ("perl "  + scriptDir + "runGATK_hardFilter.pl " + sample_name)