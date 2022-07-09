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

clean = path_config_file.clean
if clean == True:
    import shutil
    shutil.rmtree(output_mid_file_dict)
if not os.path.exists(output_mid_file_dict):
    os.mkdir(output_mid_file_dict)
    print('mid file directory made: ', output_mid_file_dict)


########## Creat directory #######################################
if not os.path.exists(f'{output_mid_file_dict}/01.vcf'):
    os.mkdir(f'{output_mid_file_dict}/01.vcf')
if not os.path.exists(f'{output_mid_file_dict}/02.vcfgz'):
    os.mkdir(f'{output_mid_file_dict}/02.vcfgz')
if not os.path.exists(f'{output_mid_file_dict}/03.psmc'):
    os.mkdir(f'{output_mid_file_dict}/03.psmc')


def run_psmc_each_sample():
    for index,sample_line in config_file.iterrows():
        print(f"No. {index} PSMC, {sample_line['sample_name']}, {sample_line['species_or_pop_name']}")
        sys.stdout.flush()
        script_str = generate_psmc_script_str(
            bam_path=sample_line['bam_file_path'],
            sample_name=sample_line['sample_name'],
            species_or_pop_name=sample_line['species_or_pop_name'],
            reference_genome_path=sample_line['reference_genome_path'],
            mutation_rate=sample_line['mutation_rate'],
            generation_length=sample_line['generation_length'],
            samtools_path=samtools_path,
            bcftools_path=samtools_path,
            psmc_path=psmc_path,
            output_mid_file_dict=output_mid_file_dict,
            vcfutils_pl_path=vcfutils_pl_path,
            min_depth=min_depth,
            fq2psmcfa_path=fq2psmcfa_path,
            sequence_quality=sequence_quality,
            psmc_N=psmc_N,
            psmc_t=psmc_t,
            psmc_r=psmc_r,
            psmc_p=psmc_p,
            perl_path=perl_path,
            psmc_plot_path=psmc_plot_path,
            psmc_util_path=psmc_util_path)

        filename = f"test_{sample_line['sample_name']}_{index}.sh"
        qsub_monitor(filename, script_str, threads_count=1, max_tasks_count=100, delete=False)


def run_psmc_pooled():
    for index,pop_name in enumerate(config_file['species_or_pop_name'].unique()):
        print(pop_name)
        sub_data = config_file[config_file['species_or_pop_name']==pop_name]
        print(sub_data)
        print(f"No. {index} PSMC, {pop_name}")
        sys.stdout.flush()
        script_str = generate_psmc_script_pooled_str(
            bam_path=sub_data['bam_file_path'],
            species_or_pop_name=pop_name,
            reference_genome_path=sub_data['reference_genome_path'].tolist()[0],
            mutation_rate=sub_data['mutation_rate'].tolist()[0],
            generation_length=sub_data['generation_length'].tolist()[0],
            samtools_path=samtools_path,
            bcftools_path=samtools_path,
            psmc_path=psmc_path,
            output_mid_file_dict=output_mid_file_dict,
            vcfutils_pl_path=vcfutils_pl_path,
            min_depth=min_depth,
            fq2psmcfa_path=fq2psmcfa_path,
            sequence_quality=sequence_quality,
            psmc_N=psmc_N,
            psmc_t=psmc_t,
            psmc_r=psmc_r,
            psmc_p=psmc_p,
            perl_path=perl_path,
            psmc_plot_path=psmc_plot_path,
            psmc_util_path=psmc_util_path)

        filename = f"test_{pop_name}_{index}.sh"
        qsub_monitor(filename, script_str, threads_count=1, max_tasks_count=100, delete=False)



if run_psmc == True:
    if run_psmc_pool == True:
        run_psmc_pooled()

    else:
        run_psmc_each_sample()

else:
    pass