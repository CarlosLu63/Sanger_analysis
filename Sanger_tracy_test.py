import os
import pandas as pd
import glob

path = "/home/carloslu/Sanger/"
sample_path = "/home/carloslu/Sanger/Sample_file/GJB2_cases/"
output_path = "/home/carloslu/Sanger/output_file/"
ref_path = "/home/carloslu/Sanger/Reference/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz"
scriptDir = "/home/carloslu/Sanger/code/"
tool_path = "/home/carloslu/Sanger/tracy/bin/tracy"
vc_parameter = " -q 50 -u 50 -p 0.33 "

###position list###
#pos_list = [20763150, 20763210, 20763294, 20763421, 20763486, 20763612, 20763614, 20763686, 20763691, 20763531]
pos_list = ["rs397516878", "rs773528125", "rs80338948", "rs111033204", "rs80338943", "rs750188782", "rs72474224", "rs587783644", "rs1801002", "rs80338939"]

###create reference index must be bgzip file###b
#os.system(tool_path + " index -o " + ref_path + "human_g1k_v37.fasta.fm9 " + ref_path + "human_g1k_v37.fasta.gz")
#os.system("samtools faidx " + ref_path + "human_g1k_v37.fasta.gz")

os.chdir(sample_path)
sanger_check = {}
file_list = glob.glob('*.ab1')

os.chdir(output_path)
for j in file_list:
    print(j) 
    input_filename = j
    output_filename = input_filename.split("_")[1].split("_")[0]
    vcf_path = output_path + output_filename + '/' + output_filename + ".vcf"
    os.mkdir(output_filename)
    ###variant calling###
    os.system(tool_path + " decompose -v -a homo_sapiens_hg19 -r " + ref_path + vc_parameter + " -o " + output_path + output_filename + '/' + output_filename + " "  + sample_path + input_filename)
    os.system("rm -f  " + output_path + output_filename + '/' + output_filename + "*.align1")
    os.system("rm -f  " + output_path + output_filename + '/' + output_filename + "*.align2")
    os.system("rm -f  " + output_path + output_filename + '/' + output_filename + "*.align3")
    os.system("rm -f  " + output_path + output_filename + '/' + output_filename + "*.decomp")
    os.system("rm -f  " + output_path + output_filename + '/' + output_filename + "*.json")
    os.system("rm -f  " + output_path + output_filename + '/' + output_filename + "*.abif")
    ###convert result to vcf###
    os.system("/opt/bcftools/bcftools view " + output_path + output_filename + '/' + output_filename + ".bcf" + " -O v > " + vcf_path)
    ###read vcf file###
    with open(vcf_path, "r") as f:
        lines = f.readlines()
        chrom_index = [i for i, line in enumerate(lines) if line.strip().startswith("#CHROM")]
        data = lines[chrom_index[0]:]  
        header = data[0].strip().split("\t")
        informations = [d.strip().split("\t") for d in data[1:]]
    vcf_out = pd.DataFrame(informations, columns = header)
    sanger_check[output_filename]={}
    for i in pos_list:
        if i in vcf_out["ID"].values.tolist():
            #print(i)
            subdata = pd.DataFrame(vcf_out[vcf_out['ID']==i][["INFO", "sample"]])
            type_ = subdata.iloc[0][0].split("=")[1].split(";")[0]
            genotype = subdata.iloc[0][1].split(":")[0]
            if genotype == '0/0':
                gt = 'Homozygous'
            elif genotype == '0/1':
                gt = 'Heterozygous'
            elif genotype == '1/1':
                gt = 'Homozygous ALT'
            else:
                gt = 'NA'
            #print(type_)
            #print(gt)
            sanger_check[output_filename][i] = [type_, gt]
        else:
            sanger_check[output_filename][i] = 'Normal'

tracy_out = pd.DataFrame.from_dict(sanger_check)
tracy_out.to_csv(r'/home/carloslu/Sanger/output_file/GJB2_Tracy_test.csv', index = True, header = True)

input_filename = "F1_GJB2-17-F_20220429_173428.ab1"
output_filename = input_filename.split("_")[1].split("_")[0]
os.system(tool_path + " basecall -f fastq " + " -o " + output_path + output_filename + '/' + output_filename + ".fastq "  + sample_path + input_filename)
os.system(tool_path + " align -r " + ref_path + " " + sample_path + input_filename)
