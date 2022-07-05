#!/beegfs/store4/chenyangkang/miniconda3/bin/python
import pandas as pd
import numpy as np
import os
import subprocess
import time
import sys
sys.stdout = open("pop.log", "a")

config_file = pd.read_csv('auto_input.csv', sep=',')
config_file.columns = ['sample_name','species_or_pop_name','generation_length',\
            'mutation_rate','mutation_rate_lower_bound','mutation_rate_upper_bound',\
                'unknow1','unknow2','reference_genome_path','bam_file_path','vcf_file_path']
# print(config_file)


##### define middle-product vcf output directory
output_mid_file_dict="/beegfs/store4/chenyangkang/22.pop_pipeline/mid_file"   #### Your ideal output directory of raw vcf and other files in the pipeline. Absolute path. For somereason, the vcf and other files will be generated as a middle-product.

#### software path:
samtools_path="/public/software/samtools-1.8/samtools"
bcftools_path="/public/software/bcftools-1.11/bin/bcftools"
psmc_path="/beegfs/store4/chenyangkang/software/psmc/psmc"
vcfutils_pl_path="/public/software/bcftools-1.11/bin/vcfutils.pl" #### usually in your bcftools package
fq2psmcfa_path="/beegfs/store4/chenyangkang/software/psmc/utils/fq2psmcfa" #### usually in your psmc package
psmc_plot_path="/beegfs/store4/chenyangkang/software/psmc/utils/psmc_plot.pl"
perl_path="/beegfs/store4/chenyangkang/miniconda3/bin/perl"
psmc_util_path = "/beegfs/store4/chenyangkang/22.pop_pipeline/psmc/utils"
java_path = "/beegfs/store4/chenyangkang/software/jdk-11.0.2/bin/java"
picard_jar_path="/beegfs/store4/chenyangkang/software/picard.jar"
gatk_path="/beegfs/store4/chenyangkang/software/gatk-4.1.9.0/gatk"

if not os.path.exists(f'{output_mid_file_dict}/04.bam'):
    os.mkdir(f'{output_mid_file_dict}/04.bam')



from generate_script_utils_call_GATKvariants import *
from generate_script_utils_for_psmc import *

def run_GATKvariants_calling_seperate_step():
    for index, line in config_file.iterrows():
        script_str = generate_GATKvariants_calling_seperate_script_str(
            java_path=java_path, 
            picard_jar_path=picard_jar_path, 
            gatk_path=gatk_path,
            input_bam_path=line['bam_file_path'], 
            output_mid_file_dict=output_mid_file_dict, 
            sample_name=line['sample_name'], 
            pop_name=line['species_or_pop_name'], 
            ref_genome=line['reference_genome_path'])
        
        filename = f"test_{line['sample_name']}_{line['species_or_pop_name']}_{index}.GATKvariants_calling_seperate.sh"
        qsub_monitor(filename, script_str, threads_count=1, max_tasks_count=100, delete=False)

    
def run_GATKvariants_calling_combined_step():
    for index, pop_name in enumerate(config_file['species_or_pop_name']):
        sub_data = config_file[config_file['species_or_pop_name']==pop_name]
        script_str = generate_GATKvariants_calling_combine_script_str(
            gatk_path=gatk_path, 
            output_mid_file_dict=output_mid_file_dict, 
            pop_name=pop_name, 
            ref_genome=sub_data['reference_genome_path'].tolist()[0])

        filename = f"test_{pop_name}_{index}.GATKvariants_calling_combine.sh"
        qsub_monitor(filename, script_str, threads_count=1, max_tasks_count=100, delete=False)

###########################
run_sep = True
run_comb = False

if run_sep == True:
    run_GATKvariants_calling_seperate_step()

if run_comb == True:
    run_GATKvariants_calling_combined_step()