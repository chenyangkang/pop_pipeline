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
psmc_util_path = "/beegfs/store4/chenyangkang/software/psmc/utils"



################################ First, run psmc ########################################################
#### psmc parameters
run_psmc = True  #### True or False. For example: do_psmc=True. For example2: do_psmc=False.
run_psmc_pool = True
min_depth=5  #### minimum depth for the sites to be taken into calculation. That's a parameter of "vcfutils.pl vcf2fq -d" in psmc. For example min_depth=5.
sequence_quality=20  #### sequence quality as a parameter of "fq2psmcfa -q". Usually 20. For example, sequence_quality=20.
psmc_N=25    #### maximum number of iterations. Default is 30.
psmc_t=15    #### maximum 2N0 coalescent time. Default is 15.
psmc_r=5    #### initial theta/rho ratio. Default if 4.
psmc_p="4+25*2+4+6"   #### pattern of parameters in psmc, 64*1 is the highest solution. default is 4+5*3+4.



from generate_script_utils_for_psmc import *


def boostrap_psmc_each_sample():
    for index,sample_line in config_file.iterrows():
        print(f"No. {index} PSMC, boostrap, {sample_line['sample_name']}, {sample_line['species_or_pop_name']}")
        sys.stdout.flush()
        script_str_list =  generate_psmc_script_str_boostrap(
            100, 
            output_mid_file_dict, 
            pooled_psmcfa=f"{output_mid_file_dict}/03.psmc/{sample_line['sample_name']}.{sample_line['species_or_pop_name']}.d{min_depth}.split.psmcfa",
            psmc_path=psmc_path, 
            sample_name=sample_line['sample_name'], 
            species_or_pop_name=sample_line['species_or_pop_name'], 
            psmc_N=psmc_N,
            psmc_t=psmc_t,
            psmc_r=psmc_r,
            psmc_p=psmc_p)

        for index2,script in enumerate(script_str_list):
            filename = f"test_{sample_line['sample_name']}_{index2}.boostrap.sh"
            qsub_monitor(filename, script, threads_count=1, max_tasks_count=100, delete=False)



def boostrap_psmc_pooled():
    for index,pop_name in enumerate(config_file['species_or_pop_name'].unique()):
        print(pop_name)
        sub_data = config_file[config_file['species_or_pop_name']==pop_name]
        print(sub_data)
        print(f"No. {index} PSMC, boostrap {pop_name}")
        sys.stdout.flush()
        script_str_list = generate_psmc_script_pooled_str_boostrap(
            100, 
            output_mid_file_dict, 
            pooled_psmcfa=f"{output_mid_file_dict}/03.psmc/{pop_name}.merged.d{min_depth}.split.psmcfa", 
            psmc_path=psmc_path, 
            species_or_pop_name=pop_name, 
            psmc_N=psmc_N, 
            psmc_t=psmc_t, 
            psmc_r=psmc_r, 
            psmc_p=psmc_p)

        for index2,script in enumerate(script_str_list):
            filename = f"test_{pop_name}_{index2}.boostrap.sh"
            qsub_monitor(filename, script, threads_count=1, max_tasks_count=100, delete=False)




################ If PSMC main are done, you can do boostrap #################
if run_psmc == True:
    if run_psmc_pool == True:
        boostrap_psmc_pooled()
    else:
        boostrap_psmc_each_sample()


