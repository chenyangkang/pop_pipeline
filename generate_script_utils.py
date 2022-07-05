
def qsub_monitor(filename, str_to_write_into_file, threads_count=1, max_tasks_count=100, delete=False):
    import os, time
    current_task_count = int(os.popen("qstat|grep chenyangkang|wc -l").read().strip())
    while current_task_count >=max_tasks_count:
        time.sleep(5)

    with open(filename,'w') as file:
        file.write(str_to_write_into_file)

    os.popen(f"chmod 755 {filename}",'w')
    os.popen(f"qsub -l nodes=1:ppn={threads_count} {filename}")
    time.sleep(0.5)
    if delete==True:
        os.remove(filename)



def generate_psmc_script_str(bam_path, sample_name, species_or_pop_name, reference_genome_path,\
    mutation_rate, generation_length, samtools_path, bcftools_path, psmc_path, output_mid_file_dict,\
        vcfutils_pl_path, min_depth, fq2psmcfa_path, sequence_quality,\
        psmc_N, psmc_t, psmc_r, psmc_p, perl_path, psmc_plot_path):
    psmc_script = f"""
cd $PBS_O_WORKDIR
### this sample
bam_path={bam_path}
base_name={sample_name}
species_or_pop_name={species_or_pop_name}
this_reference={reference_genome_path}
this_mutationrate={mutation_rate}
this_generation_length={generation_length}
samtools_path={samtools_path}
bcftools_path={bcftools_path}
psmc_path={psmc_path}
output_mid_file_dict={output_mid_file_dict}
vcfutils_pl_path={vcfutils_pl_path}
min_depth={min_depth}
fq2psmcfa_path={fq2psmcfa_path}
sequence_quality={sequence_quality}
psmc_N={psmc_N}
psmc_t={psmc_t}
psmc_r={psmc_r}
psmc_p={psmc_p}
perl_path={perl_path}
psmc_plot_path={psmc_plot_path}
""" +\
"""
######## psmc process
cd ${output_mid_file_dict}/01.vcf
${samtools_path} sort ${bam_path}|bcftools mpileup -Ou -I -f {this_reference} - |bcftools call -c -Ov > ${output_mid_file_dict}/01.vcf/${base_name}.${species_or_pop_name}.vcf
wait
cd ..

${vcfutils_pl_path} vcf2fq -d ${min_depth} ${output_mid_file_dict}/01.vcf/${base_name}.${species_or_pop_name}.vcf| gzip > ${output_mid_file_dict}/03.psmc/${base_name}.${species_or_pop_name}.d${min_depth}.psmc.fq.gz
wait

${fq2psmcfa_path} -q${sequence_quality} ${output_mid_file_dict}/03.psmc/${base_name}.${species_or_pop_name}.d${min_depth}.psmc.fq.gz > ${output_mid_file_dict}/03.psmc/${base_name}.${species_or_pop_name}.d${min_depth}.psmcfa
wait

${psmc_path} -N${psmc_N} -t${psmc_t} -r${psmc_r} -p "${psmc_p}" -o ${output_mid_file_dict}/03.psmc/${base_name}.${species_or_pop_name}.d${min_depth}.psmc ${output_mid_file_dict}/03.psmc/${base_name}.${species_or_pop_name}.d${min_depth}.psmcfa
wait

psmc_title="PSMC plot for sample"${base_name}_${species_or_pop_name}"\n g="${this_generation_length}", mu="${this_mutationrate}", depth="${min_depth}", seuqence quality="${sequence_quality}

${perl_path} ${psmc_plot_path} -u ${this_mutationrate} -g ${this_generation_length} -T \
${psmc_title} \
-P "outside" -f "Helvetica,12" \
-p ${base_name}_${species_or_pop_name} \
${output_mid_file_dict}/03.psmc/${base_name}.${species_or_pop_name}.d${min_depth}.psmc

wait
    """
    return psmc_script
    




def generate_psmc_script_pooled_str(bam_path, species_or_pop_name, reference_genome_path,\
    mutation_rate, generation_length, samtools_path, bcftools_path, psmc_path, output_mid_file_dict,\
        vcfutils_pl_path, min_depth, fq2psmcfa_path, sequence_quality,\
        psmc_N, psmc_t, psmc_r, psmc_p, perl_path, psmc_plot_path):
    import pandas
    bam_list_path = f"{output_mid_file_dict}/01.vcf/{species_or_pop_name}_bamlist.txt"
    bam_path.to_csv(bam_list_path, sep=',',header=False,index=False)
    psmc_script = f"""
cd $PBS_O_WORKDIR
### this sample
bam_list_path={bam_list_path}
species_or_pop_name={species_or_pop_name}
this_reference={reference_genome_path}
this_mutationrate={mutation_rate}
this_generation_length={generation_length}
samtools_path={samtools_path}
bcftools_path={bcftools_path}
psmc_path={psmc_path}
output_mid_file_dict={output_mid_file_dict}
vcfutils_pl_path={vcfutils_pl_path}
min_depth={min_depth}
fq2psmcfa_path={fq2psmcfa_path}
sequence_quality={sequence_quality}
psmc_N={psmc_N}
psmc_t={psmc_t}
psmc_r={psmc_r}
psmc_p={psmc_p}
perl_path={perl_path}
psmc_plot_path={psmc_plot_path}
""" +\
"""
######## psmc process
cd ${output_mid_file_dict}/01.vcf

${samtools_path} merge -b ${bam_list_path} ${species_or_pop_name}.bam.merged

${samtools_path} sort ${species_or_pop_name}.bam.merged|bcftools mpileup -Ou -I -f {this_reference} - |bcftools call -c -Ov > ${output_mid_file_dict}/01.vcf/${species_or_pop_name}.merged.vcf
wait
cd ..

${vcfutils_pl_path} vcf2fq -d ${min_depth} ${output_mid_file_dict}/01.vcf/${species_or_pop_name}.merged.vcf| gzip > ${output_mid_file_dict}/03.psmc/${species_or_pop_name}.merged.d${min_depth}.psmc.fq.gz
wait

${fq2psmcfa_path} -q${sequence_quality} ${output_mid_file_dict}/03.psmc/${species_or_pop_name}.merged.d${min_depth}.psmc.fq.gz > ${output_mid_file_dict}/03.psmc/${species_or_pop_name}.merged.d${min_depth}.psmcfa
wait

${psmc_path} -N${psmc_N} -t${psmc_t} -r${psmc_r} -p "${psmc_p}" -o ${output_mid_file_dict}/03.psmc/${species_or_pop_name}.merged.d${min_depth}.psmc ${output_mid_file_dict}/03.psmc/${species_or_pop_name}.merged.d${min_depth}.psmcfa
wait

psmc_title="PSMC plot for sample"${species_or_pop_name}_merged_"\n g="${this_generation_length}", mu="${this_mutationrate}", depth="${min_depth}", seuqence quality="${sequence_quality}

${perl_path} ${psmc_plot_path} -u ${this_mutationrate} -g ${this_generation_length} -T \
${psmc_title} \
-P "outside" -f "Helvetica,12" \
-p ${species_or_pop_name}_merged \
${output_mid_file_dict}/03.psmc/${species_or_pop_name}_merged.d${min_depth}.psmc

wait
    """
    return psmc_script