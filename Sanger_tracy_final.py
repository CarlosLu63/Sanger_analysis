import os
import pandas as pd
import numpy as np

path = "/home/carloslu/Sanger/"
sample_path = "/home/carloslu/Sanger/Sample_file/GJB2_cases/"
output_path = "/home/carloslu/Sanger/output_file/"
ref_path = "/home/carloslu/Sanger/Reference/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz"
scriptDir = "/home/carloslu/Sanger/code/"
tool_path = "/home/carloslu/Sanger/tracy/bin/tracy"
bcftools_path = "/opt/bcftools/bcftools"
input_filename = "A1_GJB2-1-F_20220308_151446.ab1"
output_filename = input_filename.split("_")[1].split("_")[0]
vcf_path = output_path + output_filename + "/" + output_filename + ".vcf"
vc_parameter = " -q 50 -u 50 -p 0.33 "
os.chdir(path)

###position list###
#pos_list = [20763150, 20763210, 20763294, 20763421, 20763486, 20763612, 20763614, 20763686, 20763691, 20763531]
pos_list = ["rs397516878", "rs773528125", "rs80338948", "rs111033204", "rs80338943", "rs750188782", "rs72474224", "rs587783644", "rs1801002", "rs80338939"]

###create reference index must be bgzip file###
#os.system(tool_path + " index -o " + ref_path + "human_g1k_v37.fasta.fm9 " + ref_path + "human_g1k_v37.fasta.gz")
#os.system("samtools faidx " + ref_path + "human_g1k_v37.fasta.gz")

#############################################
###variant calling only forward or reverse###
#############################################

os.system(tool_path + " decompose -v -a homo_sapiens_hg19 -r " + ref_path + vc_parameter + " -o " + output_path + output_filename + "/" + output_filename + " "  + sample_path + input_filename)
os.system("rm -f " + output_path +  output_filename + "/" + "*.align1")
os.system("rm -f " + output_path +  output_filename + "/" + "*.align2")
os.system("rm -f " + output_path +  output_filename + "/" + "*.align3")
os.system("rm -f " + output_path +  output_filename + "/" + "*.decomp")
os.system("rm -f " + output_path +  output_filename + "/" + "*.json")
os.system("rm -f " + output_path +  output_filename + "/" + "*.abif")

###convert result to vcf###
os.system(bcftools_path + " view " + output_path + output_filename + "/" + output_filename + ".bcf" + " -O v > " + vcf_path)

###read vcf file###
with open(vcf_path, "r") as f:
    lines = f.readlines()
    chrom_index = [i for i, line in enumerate(lines) if line.strip().startswith("#CHROM")]
    data = lines[chrom_index[0]:]  
    header = data[0].strip().split("\t")
    informations = [d.strip().split("\t") for d in data[1:]]

vcf_out = pd.DataFrame(informations, columns = header)

###output all result###
sanger_allres = {}

for n in range(vcf_out.shape[0]):
    sanger_allres[n] = {}
    sanger_allres[n]['CHROM'] = vcf_out.iloc[n][0]
    sanger_allres[n]['POS'] = vcf_out.iloc[n][1]
    sanger_allres[n]['ID'] = vcf_out.iloc[n][2]
    sanger_allres[n]['REF'] = vcf_out.iloc[n][3]
    sanger_allres[n]['ALT'] = vcf_out.iloc[n][4]
    sanger_allres[n]['QUAL'] = vcf_out.iloc[n][5]
    sanger_allres[n]['FILTER'] = vcf_out.iloc[n][6]
    sanger_allres[n]['TYPE'] = vcf_out.iloc[n][7].split("=")[1].split(";")[0]
    genotype = vcf_out.iloc[n][9].split(":")[0]
    if genotype == '0/0':
        sanger_allres[n]['GT'] = 'Homozygous'
    elif genotype == '0/1':
        sanger_allres[n]['GT'] = 'Heterozygous'
    elif genotype == '1/1':
        sanger_allres[n]['GT'] = 'Homozygous ALT'
    else:
        sanger_allres[n]['GT'] = 'NA'
sanger_all = pd.DataFrame.from_dict(sanger_allres).transpose()
sanger_all.to_csv(output_path + output_filename + "/" + output_filename + "_allvariant.csv", index = False, header = True)    

###checking 10 positions###
sanger_check = {}

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
        sanger_check[i] = [type_, gt]
    else:
        sanger_check[i] = 'Normal'

print(sanger_check)



###########################################
###variant calling for forward + reverse###
###########################################

input_file_F = "F1_GJB2-17-F_20220429_173428.ab1"
input_file_R = "A3_GJB2-17_20220325_172857.ab1"
output_filename = input_file_R.split("_")[1].split("_")[0]
vcf_path_FR = output_path + output_filename + "/" + "merged.vcf"

os.system(tool_path + " decompose -o " + output_path + output_filename + "/" + "forward -a homo_sapiens_hg19 -r " + ref_path + " " + sample_path + input_file_F)
os.system(tool_path + " decompose -o " + output_path + output_filename + "/" + "reverse -a homo_sapiens_hg19 -r " + ref_path + " " + sample_path + input_file_R)

os.system(bcftools_path + " norm -O b --write-index -o " + output_path + output_filename + "/" + "forward.norm.bcf -f " + ref_path + " " + output_path + output_filename + "/" + "forward.bcf")
os.system(bcftools_path + " norm -O b --write-index -o " + output_path + output_filename + "/" + "reverse.norm.bcf -f " + ref_path + " " + output_path + output_filename + "/" + "reverse.bcf")

os.system("rm -f " + output_path +  output_filename + "/" + "*.align1")
os.system("rm -f " + output_path +  output_filename + "/" + "*.align2")
os.system("rm -f " + output_path +  output_filename + "/" + "*.align3")
os.system("rm -f " + output_path +  output_filename + "/" + "*.decomp")
os.system("rm -f " + output_path +  output_filename + "/" + "*.json")
os.system("rm -f " + output_path +  output_filename + "/" + "*.abif")

os.system(bcftools_path + " merge --force-samples " + output_path + output_filename + "/" + "forward.norm.bcf" + " " + output_path + output_filename + "/" + "reverse.norm.bcf" + " -O v > " + vcf_path_FR)

###read vcf file###
with open(vcf_path_FR, "r") as f:
    lines = f.readlines()
    chrom_index = [i for i, line in enumerate(lines) if line.strip().startswith("#CHROM")]
    data = lines[chrom_index[0]:]  
    header = data[0].strip().split("\t")
    informations = [d.strip().split("\t") for d in data[1:]]

vcf_out = pd.DataFrame(informations, columns = header)

###output all result###
sanger_allres = {}

for n in range(vcf_out.shape[0]):
    sanger_allres[n] = {}
    sanger_allres[n]['CHROM'] = vcf_out.iloc[n][0]
    sanger_allres[n]['POS'] = vcf_out.iloc[n][1]
    sanger_allres[n]['ID'] = vcf_out.iloc[n][2]
    sanger_allres[n]['REF'] = vcf_out.iloc[n][3]
    sanger_allres[n]['ALT'] = vcf_out.iloc[n][4]
    sanger_allres[n]['QUAL'] = vcf_out.iloc[n][5]
    sanger_allres[n]['FILTER'] = vcf_out.iloc[n][6]
    sanger_allres[n]['TYPE'] = vcf_out.iloc[n][7].split("=")[1].split(";")[0]
    genotype1 = vcf_out.iloc[n][9].split(":")[0]
    genotype2 = vcf_out.iloc[n][10].split(":")[0]
    if genotype1 == '0/0' or genotype2 == '0/0':
        sanger_allres[n]['GT'] = 'Homozygous'
    elif genotype1 == '0/1' or genotype2 == '0/1':
        sanger_allres[n]['GT'] = 'Heterozygous'
    elif genotype1 == '1/1' or genotype2 == '1/1':
        sanger_allres[n]['GT'] = 'Homozygous ALT'
    else:
        sanger_allres[n]['GT'] = 'NA'
sanger_all = pd.DataFrame.from_dict(sanger_allres).transpose()
sanger_all.to_csv(output_path + output_filename + "/" + output_filename + "_allvariant.csv", index = False, header = True)    

###checking 10 positions###
sanger_check = {}

for i in pos_list:
    if i in vcf_out["ID"].values.tolist():
        #print(i)
        subdata = pd.DataFrame(vcf_out[vcf_out['ID']==i][["INFO", "sample", "2:sample"]])
        type_ = subdata.iloc[0][0].split("=")[1].split(";")[0]
        genotype1 = subdata.iloc[0][1].split(":")[0]
        genotype2 = subdata.iloc[0][2].split(":")[0]
        if genotype1 == '0/0' or genotype2 == '0/0':
            gt = 'Homozygous'
        elif genotype1 == '0/1' or genotype2 == '0/1':
            gt = 'Heterozygous'
        elif genotype1 == '1/1' or genotype2 == '1/1':
            gt = 'Homozygous ALT'
        else:
            gt = 'NA'
        #print(type_)
        #print(gt)
        sanger_check[i] = [type_, gt]
    else:
        sanger_check[i] = 'Normal'

print(sanger_check)
