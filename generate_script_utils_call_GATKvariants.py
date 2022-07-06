def generate_GATKvariants_calling_seperate_script_str(java_path, picard_jar_path, gatk_path,
    input_bam_path, output_mid_file_dict, sample_name, pop_name, ref_genome):
    script_str=f"""
cd $PBS_O_WORKDIR
cd {output_mid_file_dict}/04.bam

{java_path} -jar {picard_jar_path}  AddOrReplaceReadGroups \
    VALIDATION_STRINGENCY=SILENT \
    I={input_bam_path} \
    O={output_mid_file_dict}/04.bam/{sample_name}.{pop_name}.renamed.bam \
    SORT_ORDER=coordinate \
    RGID={sample_name} \
    RGLB={sample_name} \
    RGPL=illumina \
    RGSM={sample_name} \
    RGPU={sample_name} \
    CREATE_INDEX=True
wait

{java_path} -jar {picard_jar_path} MarkDuplicates REMOVE_DUPLICATES=true  \
    VALIDATION_STRINGENCY=SILENT \
    MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 \
    INPUT={output_mid_file_dict}/04.bam/{sample_name}.{pop_name}.renamed.bam  \
    OUTPUT={output_mid_file_dict}/04.bam/{sample_name}.{pop_name}.renamed.marked.bam \
    METRICS_FILE={output_mid_file_dict}/04.bam/{sample_name}.{pop_name}.renamed.marked.bam.metrics
wait

samtools index {output_mid_file_dict}/04.bam/{sample_name}.{pop_name}.renamed.marked.bam
wait

{gatk_path} HaplotypeCaller \
    -R {ref_genome} \
    -I {output_mid_file_dict}/04.bam/{sample_name}.{pop_name}.renamed.marked.bam \
    -O {output_mid_file_dict}/02.vcfgz/{sample_name}.{pop_name}.g.vcf.gz \
    -ERC GVCF \
    -ploidy 2
    """
    return script_str


def generate_GATKvariants_calling_combine_script_str(gatk_path, output_mid_file_dict, pop_name, ref_genome):
    script_str=f"""
cd $PBS_O_WORKDIR
cd {output_mid_file_dict}/02.vcfgz
ls {output_mid_file_dict}/02.vcfgz/*.{pop_name}.g.vcf.gz > {output_mid_file_dict}/02.vcfgz/{pop_name}.gvcf.list
{gatk_path} CombineGVCFs \
    -R {ref_genome} -V {output_mid_file_dict}/02.vcfgz/{pop_name}.gvcf.list \
    -O {output_mid_file_dict}/02.vcfgz/{pop_name}.g.vcf.gz

{gatk_path} GenotypeGVCFs -R {ref_genome} -V {output_mid_file_dict}/02.vcfgz/{pop_name}.g.vcf.gz -O {output_mid_file_dict}/02.vcfgz/{pop_name}.genotype.vcf.gz

{gatk_path} VariantFiltration \
    -V {output_mid_file_dict}/02.vcfgz/{pop_name}.genotype.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O {output_mid_file_dict}/02.vcfgz/{pop_name}.genotype.filtered.vcf.gz

tabix -p vcf {output_mid_file_dict}/02.vcfgz/{pop_name}.genotype.filtered.vcf.gz
    """
    return script_str