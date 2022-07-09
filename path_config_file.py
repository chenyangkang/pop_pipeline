
import pandas as pd
import numpy as np
import os
import subprocess
import time
import sys

###########
clean=False  ### do not remove previous files

#############
config_file = pd.read_csv('auto_input.csv', sep=',')
config_file.columns = ['sample_name','species_or_pop_name','generation_length',\
            'mutation_rate','mutation_rate_lower_bound','mutation_rate_upper_bound',\
                'unknow1','unknow2','reference_genome_path','bam_file_path','vcf_file_path']

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
java_path = "/beegfs/store4/chenyangkang/software/jdk-11.0.2/bin/java"
picard_jar_path="/beegfs/store4/chenyangkang/software/picard.jar"
gatk_path="/beegfs/store4/chenyangkang/software/gatk-4.1.9.0/gatk"



################################ First, run psmc ########################################################
#### psmc parameters
run_psmc = True  #### True or False. For example: do_psmc=True. For example2: do_psmc=False.
run_psmc_pool = True
min_depth=10  #### minimum depth for the sites to be taken into calculation. That's a parameter of "vcfutils.pl vcf2fq -d" in psmc. For example min_depth=5.
sequence_quality=20  #### sequence quality as a parameter of "fq2psmcfa -q". Usually 20. For example, sequence_quality=20.
psmc_N=25    #### maximum number of iterations. Default is 30.
psmc_t=15    #### maximum 2N0 coalescent time. Default is 15.
psmc_r=5    #### initial theta/rho ratio. Default if 4.
psmc_p="4+25*2+4+6"   #### pattern of parameters in psmc, 64*1 is the highest solution. default is 4+5*3+4.

