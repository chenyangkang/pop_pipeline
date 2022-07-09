#!/beegfs/store4/chenyangkang/miniconda3/bin/python
import pandas as pd
import numpy as np
import os
import subprocess
import time
import sys


import path_config_file
from generate_script_utils_call_GATKvariants import *
from generate_script_utils_for_psmc import *

sys.stdout = open("pop.log", "a")

config_file = path_config_file.config_file
config_file.columns = path_config_file.config_file.columns
output_mid_file_dict = path_config_file.output_mid_file_dict
samtools_path = path_config_file.samtools_path
bcftools_path = path_config_file.bcftools_path
psmc_path=path_config_file.psmc_path
vcfutils_pl_path=path_config_file.vcfutils_pl_path
fq2psmcfa_path=path_config_file.fq2psmcfa_path
psmc_plot_path=path_config_file.psmc_plot_path
perl_path=path_config_file.perl_path
psmc_util_path = path_config_file.psmc_util_path
java_path = path_config_file.java_path
picard_jar_path=path_config_file.picard_jar_path
gatk_path=path_config_file.gatk_path



################################ First, run psmc ########################################################
#### psmc parameters
run_psmc = path_config_file.run_psmc  #### True or False. For example: do_psmc=True. For example2: do_psmc=False.
run_psmc_pool = path_config_file.run_psmc_pool
min_depth=path_config_file.min_depth  #### minimum depth for the sites to be taken into calculation. That's a parameter of "vcfutils.pl vcf2fq -d" in psmc. For example min_depth=5.
sequence_quality=path_config_file.sequence_quality  #### sequence quality as a parameter of "fq2psmcfa -q". Usually 20. For example, sequence_quality=20.
psmc_N=path_config_file.psmc_N  #### maximum number of iterations. Default is 30.
psmc_t=path_config_file.psmc_t    #### maximum 2N0 coalescent time. Default is 15.
psmc_r=path_config_file.psmc_r    #### initial theta/rho ratio. Default if 4.
psmc_p=path_config_file.psmc_p  #### pattern of parameters in psmc, 64*1 is the highest solution. default is 4+5*3+4.



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


