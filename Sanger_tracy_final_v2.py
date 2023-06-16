import os
import pandas as pd
import glob

class Sanger_analysis(object):

    def main_flow(self):

        ###screen ab1 file###
        input_file = glob.glob(self.sample_path + "*.ab1")[0]

        ###setup path name###
        vcf_path, input_filename, output_filename = self.read_ab1(input_file)

        ###apply variant calling###
        self.variant_calling(output_filename, input_filename)

        ###convert bcf to vcf###
        vcf_out = self.bcf_to_vcf(output_filename, vcf_path)
        
        ###output all variant as csv file###
        self.output_all_var(vcf_out, output_filename)
        
        ###check 10 pos###
        sanger_check = self.check_pos(vcf_out)

        ###print dictionary###
        print(sanger_check)


    def __init__(self):
        self.path = "/home/carloslu/Sanger/"
        self.sample_path = "/home/carloslu/Sanger/Sample_file/GJB2_cases/"
        self.output_path = "/home/carloslu/Sanger/output_file/"
        self.ref_path = "/home/carloslu/Sanger/Reference/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz"
        self.scriptDir = "/home/carloslu/Sanger/code/"
        self.tool_path = "/home/carloslu/Sanger/tracy/bin/tracy"
        self.bcftools_path = "/opt/bcftools/bcftools"
        self.vc_parameter = " -q 50 -u 50 -p 0.33 "
        ###position list###
        self.pos_list = pd.read_excel(self.path + 'position_list.xlsx')['pos_list'].to_list()


    def read_ab1(self, input_file):
        input_filename = input_file.split("/")[-1]
        output_filename = input_filename.split("_")[1].split("_")[0]
        vcf_path = self.output_path + output_filename + "/" + output_filename + ".vcf"
        return vcf_path, input_filename, output_filename


    def variant_calling(self, output_filename, input_filename):
        if os.path.exists(self.output_path + output_filename) == False:
            os.mkdir(self.output_path + output_filename)

        os.system(self.tool_path + " decompose -v -a homo_sapiens_hg19 -r " + self.ref_path + self.vc_parameter + " -o " + self.output_path + output_filename + "/" + output_filename + " "  + self.sample_path + input_filename)
        os.system("rm -f " + self.output_path +  output_filename + "/" + "*.align1")
        os.system("rm -f " + self.output_path +  output_filename + "/" + "*.align2")
        os.system("rm -f " + self.output_path +  output_filename + "/" + "*.align3")
        os.system("rm -f " + self.output_path +  output_filename + "/" + "*.decomp")
        os.system("rm -f " + self.output_path +  output_filename + "/" + "*.json")
        os.system("rm -f " + self.output_path +  output_filename + "/" + "*.abif")


    def bcf_to_vcf(self, output_filename, vcf_path):
        ###convert result to vcf###
        os.system(self.bcftools_path + " view " + self.output_path + output_filename + "/" + output_filename + ".bcf" + " -O v > " + vcf_path)

        ###read vcf file###
        with open(vcf_path, "r") as f:
            lines = f.readlines()
            chrom_index = [i for i, line in enumerate(lines) if line.strip().startswith("#CHROM")]
            data = lines[chrom_index[0]:]  
            header = data[0].strip().split("\t")
            informations = [d.strip().split("\t") for d in data[1:]]

        vcf_out = pd.DataFrame(informations, columns = header)
        return vcf_out


    ###output all result###
    def output_all_var(self, vcf_out, output_filename):
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
        sanger_all.to_csv(self.output_path + output_filename + "/" + output_filename + "_allvariant.csv", index = False, header = True)    


    ###checking 10 positions###
    def check_pos(self, vcf_out):
        sanger_check = {}
        for i in self.pos_list:
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
        return sanger_check

def main():	 
    main_obj = Sanger_analysis()
    main_obj.main_flow()

if __name__ == "__main__":
### simple main
    main()
