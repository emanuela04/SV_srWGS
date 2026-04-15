nextflow.enable.dsl = 2
/*
 * Pipeline SVs detection for SR WGS sequencing data
 * Autor: Emanuela Iovino emanuela.iovino@aosp.bo.it
 * This pipeline relies on Nexflow and it works using Nexflow version >= 21.10.6.5660
 */


if (params.help) {
    log.info """
    -----------------------------------------------------------------------
    🧬  Structural Variant Detection Pipeline for SR-WGS data
    ----------------------------------------------------------
    This pipeline detects, filters, and merges SVs from short-read
    WGS data using multiple tools and custom processing steps.

    Usage :
    ______________________________________________________________________

    Required:
    --sample_alignments_tsv      PATH              Path to a TSV file with two columns: the sample ID and the full (absolute) path to the CRAM file
    --reference_fasta            PATH              Reference genome to which the reads are aligned.
    --singularity_cache          PATH              Path to singularity cache directory
    --delly_exclude_regions_bed	 PATH		   Path to bad regiosn used from Delly
    --smoove_exclude_regions_bed PATH		   Path to bed regions used from Smoove
    --expansion_hunter_variant_catalog_json        Path to STR variant catalog file: choose between ExpansionHunter_variant_catalog.json (174k polymorphic loci) or variant_catalog_genes.json (30 known pathogenic loci)
     
    --account_name               NAME              BC+user name
    --profile           NAME              slurm

    """.stripIndent()
    exit 0
}


println '''
🧬  Structural Variant Detection Pipeline 🧬
-------------------------------------------
Detects structural variants from SR-WGS data using multiple tools.
'''


data = channel
        .fromPath(params.sample_alignments_tsv, type: "file", checkIfExists: true)
        .splitCsv(sep: "\t", header: ["sample_id", "alignment_file"])
        .map { row -> tuple(row.sample_id, row.alignment_file) }

reference_fasta = file(params.reference_fasta, type: "file", checkIfExists: true)
reference_fasta_fai = file("${reference_fasta}.fai", checkIfExists: true)
delly_exclude_regions_bed = file(params.delly_exclude_regions_bed, type: "file", checkIfExists: true)
smoove_exclude_regions_bed = file(params.smoove_exclude_regions_bed, type: "file", checkIfExists: true)
expansion_hunter_variant_catalog_json = file(params.expansion_hunter_variant_catalog_json, type: "file", checkIfExists: true)
recode_delly_python_script = file(params.recode_delly_python_script, type: "file", checkIfExists: true)
mosdepth_segmental_duplications_bed = file(params.mosdepth_segmental_duplications_bed, type:"file", checkIfExists: true)
samtools_path = params.samtools_path

manta_inv = file(params.manta_inv, type: "file", checkIfExists: true)
repeat_masker = file(params.repeat_masker, type: "file", checkIfExists: true)
seg_duplications = file(params.seg_duplications, type: "file", checkIfExists: true)

// Run mosdepth for getting coverage statistics.

process COVERAGE {
    publishDir "results/$sample_id/", mode: "copy"
    input:
    tuple val(sample_id), val(alignment_file) 
    output:
    path("${sample_id}.mosdepth.summary.txt") 
    path("${sample_id}.mosdepth.global.dist.txt") 
    script:
    """
    mosdepth --threads $task.cpus --no-per-base --fast-mode --fasta $reference_fasta $sample_id $alignment_file
    mosdepth --threads $task.cpus --by $mosdepth_segmental_duplications_bed --use-median --mapq 30 --fasta $reference_fasta ${sample_id}.segdups $alignment_file
    """
    stub:
    """
    touch ${sample_id}.mosdepth.summary.txt
    touch ${sample_id}.mosdepth.global.dist.txt   
    """
}

process cnvpytor {                                                              
    publishDir "results/$sample_id/", mode:"copy"                               
    input:                                                                      
    tuple val(sample_id), val(alignment_file)                                   
    output:                                                                     
    tuple val(sample_id), path("cnvpytor/$sample_id*.tsv")                      
    script:                                                                     
    """                                                                         
    cnvpytor -root ${sample_id}.pytor -chrom echo chr{1..22} -rd $alignment_file -T $reference_fasta
    cnvpytor -root ${sample_id}.pytor -his 1000                                 
                                                                                
    cnvpytor -root ${sample_id}.pytor -partition 1000                           
    cnvpytor -root ${sample_id}.pytor  -call 1000 > ${sample_id}.tsv            
    cut -f 2  ${sample_id}.tsv | cnvpytor -root ${sample_id}.pytor  -genotype 1000
    #perl Hg38/cnvnator2VCF.pl -prefix $sample_id -reference GRCh38 ${sample_id}.tsv $reference_fasta > ${sample_id}_CNV.vcf
    """                                                                         
    stub:                                                                       
    """                                                                         
    #touch ${sample_id}_CNV.vcf  
    touch ${sample_id}.tsv                                                
    """                                                                         
}     

// Run the Delly variant caller.
process DELLY {
    publishDir "results/$sample_id/", mode: "copy"
    input:
    tuple val(sample_id), val(alignment_file)
    output:
    tuple val(sample_id), path("$sample_id-delly.bcf")
    script:
    """
    delly call -g $reference_fasta -x $delly_exclude_regions_bed  -o tmp.bcf $alignment_file -q 20 -s 15 -z 5
    delly call -g $reference_fasta -v tmp.bcf -x $delly_exclude_regions_bed -o $sample_id-delly.bcf  $alignment_file -q 20 -s 15 -z 5
    rm tmp.bcf
    """
    stub:
    """
    touch $sample_id-delly.bcf
    """
}
//// Recode the Delly .bcf file and convert it into a .vcf file.
process RECODE {
    input:
    tuple val(sample_id), path(output_delly)
    output:
    tuple val(sample_id), path("$sample_id-delly-recode.vcf")
    script:
    """
    python $recode_delly_python_script $output_delly ${sample_id}-delly-recode.vcf
    """
    stub:
    """
    touch ${sample_id}-delly-recode.vcf
    """
}


// Run the Manta variant caller.


process MANTA {
    publishDir "results/$sample_id/", mode: "copy"
    input:
    tuple val(sample_id), val(alignment_file)
    output:
    tuple val(sample_id), path("$sample_id-manta.vcf.gz")
    script:
    """
    configManta.py --bam $alignment_file --referenceFasta $reference_fasta --runDir run_folder/ 
    cd run_folder
    python runWorkflow.py
    python  $manta_inv $samtools_path $reference_fasta results/variants/diploidSV.vcf.gz > results/variants/diploidSV_inv.vcf 
    mv results/variants/diploidSV_inv.vcf ../$sample_id-manta.vcf.gz
    """
    stub:
    """
    touch $sample_id-diploidSV.vcf.gz
    """
}

// Run the Smoove variant caller.
process SMOOVE {
    publishDir "results/$sample_id/", mode: "copy"
    input:
    tuple val(sample_id), val(alignment_file)
    output:
    tuple val(sample_id), path("$sample_id-smoove.genotyped.vcf.gz")
    script:
    """
    smoove call --name $sample_id --exclude $smoove_exclude_regions_bed --fasta $reference_fasta  -p $task.cpus --duphold --genotype $alignment_file 
    """
    stub:
    """
    touch $sample_id-smoove.genotyped.vcf.gz
    """
}

// Filter VCF with bcftools (filter on SVLEN and PASS flag).
process FILTER {
    input:
    tuple val(sample_id), path(vcf)
    output: 
    tuple val(sample_id), path("${sample_id}_$vcf-filtered.vcf")
    
    publishDir "results/$sample_id/", mode: "copy"  
    
    script:  
    """               
    python ${projectDir}/Hg38/filter_vcfs.py -i $vcf -o ${sample_id}_$vcf-filtered.vcf -s ${sample_id}_$vcf-skipped.tsv
                                                              
    """
    stub:
    """
    touch ${sample_id}_$vcf-filtered.vcf 
    """
}

// Create INV vcf file from Manta
process INV {
    publishDir "results/$sample_id/", mode: "copy"
    input:
    tuple val(sample_id), path(vcf)
    output:
    tuple val(sample_id), path("${sample_id}_$vcf-inv.vcf")
    script:
    """
    bcftools view --threads $task.cpus --include 'INFO/SVTYPE="INV"' $vcf > ${sample_id}_$vcf-inv.vcf
    """
    stub:
    """
    touch ${sample_id}_$vcf-inv.vcf
    """
}

// Run SURVIVOR to merge the within-sample variant calls and metrics


process MERGE {
    publishDir "results/$sample_id/", mode: "copy"
    input:
    val(sample_id)
    tuple path(vcf1), path(vcf2), path(vcf3)
    output: 
    tuple val(sample_id), path("${sample_id}_merged_raw.vcf")
    script:
    """
    ls $vcf2    >  vcf_list.txt
    ls $vcf1    >> vcf_list.txt
    ls $vcf3    >> vcf_list.txt
    grep manta vcf_list.txt > vcf_list.sorted.txt
    grep delly vcf_list.txt >> vcf_list.sorted.txt
    grep smoove vcf_list.txt >> vcf_list.sorted.txt
    jasmine file_list=vcf_list.sorted.txt out_file=${sample_id}_merged_raw.vcf max_dist=1000 --ignore_strand 
    """
    stub:
    """
    
    touch ${sample_id}_merged_raw.vcf
    """
}



process NORMALIZATION {
   publishDir "results/$sample_id/", mode: "copy"

   input:
   tuple val(sample_id), path(merged_raw_vcf)
   output:
   tuple val(sample_id), path("${sample_id}_merged.vcf")

   script:
   """
   python ${projectDir}/Hg38/normalize.py -i ${merged_raw_vcf} -o ${sample_id}_merged.vcf  --symbolic-ins   --report-json normalization_report.json   --changes-tsv normalization_changes.tsv 
   """
   stub:
   """
   touch ${sample_id}_merged.vcf
   """

}




process SV_QC {
   publishDir "results/$sample_id/", mode: "copy"
   
   input:
   tuple val(sample_id), path(merged_vcf)
   
   output:
   tuple val(sample_id), path("${sample_id}_SVLEN_summary_mqc.stats"), path("${sample_id}_SVs_summary_mqc.stats")
   
   script: 
   """
    SURVIVOR stats ${sample_id}_merged.vcf -1 -1 -1 ${sample_id}_stats
    mv ${sample_id}_stats  ${sample_id}_SVLEN_summary_mqc.stats  
    mv ${sample_id}_stats_CHR  ${sample_id}_SVs_summary_mqc.stats 
   """
   
   stub: 
   """
   touch ${sample_id}_SVLEN_summary_mqc.stats
   touch ${sample_id}_SVs_summary_mqc.stats
   """ 
}




 process EXPANSION_HUNTER {
   publishDir "results/$sample_id/", mode: "copy"
   input:                                                                   
   tuple val(sample_id), val(alignment_file)
   output: 
   path("${sample_id}.expansion_hunter.vcf")
   script:
    """
    ExpansionHunter --reads $alignment_file --reference $reference_fasta --variant-catalog $expansion_hunter_variant_catalog_json  --output-prefix ${sample_id}.expansion_hunter --analysis-mode streaming --threads $task.cpus
    """
    stub:
    """
    touch ${sample_id}.eh.vcf
    """
 }

process SV_COVERAGE_EXTRACT {
                                                                        
    publishDir "results/$sample_id/", mode:"copy"                               
    input:
    tuple val(sample_id), path(merged_vcf)
    output:
    tuple val(sample_id), path("${sample_id}.sv.bed")                                 
    script:
    """
    bcftools query  -f '%CHROM\t%POS0\t%END\t%SVTYPE\t%SVLEN\n' $merged_vcf > ${sample_id}.sv.bed
    
    """
    
    stub:
    """
    ${sample_id}.sv.bed
    """

 
 }

process SV_COVERAGE {
    publishDir "results/$sample_id/", mode: "copy"                              
    input:                                                                      
    tuple val(sample_id), val(alignment_file), path(sample_sv_bed)
    output:
    tuple val(sample_id), path("${sample_id}.sv.regions.bed.gz") 
    script:
    """
    mosdepth --threads $task.cpus  --by $sample_sv_bed --use-median --mapq 30  --fasta $reference_fasta ${sample_id}.sv $alignment_file
    """
    stub:
    
    """
    ${sample_id}.sv.regions.bed.gz
    """
 
 }


process CG_CONTENT {
    publishDir "results/$sample_id/", mode: "copy"                              
    input:                                                                      
    tuple val(sample_id), path(sample_sv_bed)
    output:
    tuple val(sample_id), path("${sample_id}.sv.CG.txt")
    script:
    """
    bedtools nuc -fi $reference_fasta -bed $sample_sv_bed > ${sample_id}.sv.CG.txt
    """
    stub:
    
    """
    ${sample_id}.sv.CG.txt
    """
 
 }


process EXTRACT_FEATURES {
   
    publishDir "results/$sample_id/", mode: "copy"                              
    input:                                                                      
    tuple val(sample_id),val (alignment_file), path(merged_vcf), path(sample_sv_CG), path(sample_coverage)
    output:
    tuple val(sample_id), path("${sample_id}_features.csv"), path("${sample_id}_features.tsv"), path("${sample_id}_features.tsv_DEL.csv"), path("${sample_id}_features.tsv_DUP.csv"), path("${sample_id}_features.tsv_INS.csv"), path("${sample_id}_features.tsv_INV.csv")
               
    script:
    """
    echo "SAMPLE: $sample_id"
    echo "ALIGNMENT: $alignment_file"
    echo "VCF: $merged_vcf"
    echo "GC: $sample_sv_CG"
    echo "COV: $sample_coverage"
    python ${projectDir}/Hg38/speedxx.py --cram_file $alignment_file --vcf_file $merged_vcf --out_file ${sample_id}_extra_cram.tsv
    
    python ${projectDir}/Hg38/merge_fixed.py --repeat_masker $repeat_masker --segmental_duplications $seg_duplications --vcf_file $merged_vcf  --gc_content  $sample_sv_CG --coverage_mosdepth $sample_coverage --sv_feature_matrix ${sample_id}_extra_cram.tsv --output_path ${sample_id}_features.tsv
    python ${projectDir}/Hg38/convert_sv.py ${sample_id}_features.tsv  -o ${sample_id}_features.tsv  
 
    """
    stub:
    
    """
    ${sample_id}_features.tsv
    ${sample_id}_features.csv
         
    """

   
}



process VARTRUSTML {
    publishDir "results/$sample_id/", mode: "copy"

    input:
    tuple val(sample_id), path(del_file), path(dup_file), path(ins_file), path(inv_file)

    output:
    tuple val(sample_id),
          path("${sample_id}_prediction_del.csv"),
          path("${sample_id}_prediction_dup.csv"),
          path("${sample_id}_prediction_ins.csv"),
          path("${sample_id}_prediction_inv.csv")

    script:
    """
    mkdir prediction
    cp ${del_file} prediction/
    cp ${dup_file} prediction/
    cp ${ins_file} prediction/
    cp ${inv_file} prediction/
    vartrustml predict ${projectDir}/Hg38/VarTrustML_model/fitted_model_DEL/CatBoost_calibrated_model.joblib prediction/${del_file} --output ${sample_id}_prediction_del.csv --proba
    vartrustml predict ${projectDir}/Hg38/VarTrustML_model/fitted_model_DUP/CatBoost_calibrated_model.joblib prediction/${dup_file} --output ${sample_id}_prediction_dup.csv --proba
    vartrustml predict ${projectDir}/Hg38/VarTrustML_model/fitted_model_INS/CatBoost_calibrated_model.joblib prediction/${ins_file} --output ${sample_id}_prediction_ins.csv --proba
    vartrustml predict ${projectDir}/Hg38/VarTrustML_model/fitted_model_INV/CatBoost_calibrated_model.joblib prediction/${inv_file} --output ${sample_id}_prediction_inv.csv --proba
    """
}







process MULTIQC {
    publishDir "results/multiqc", mode: "copy"

    input:
    path(qc_files)
    tuple val(sample_id), path("${sample_id}_SVLEN_summary_mqc.stats"), path("${sample_id}_SVs_summary_mqc.stats")
 
    
    output:
    path "SVs-report_multiqc_report.html" 
    path "SVs-report_multiqc_report_data"

    script:
    """
    multiqc --force --config ${projectDir}/Hg38/multiqc_config_sv.yml .
    """
}


workflow {
    // Get coverage using Mosdepth 
    coverage_out = COVERAGE(data)
    coverage_reports = coverage_out[0]   // summary.txt
                           .mix(coverage_out[1])
                           .collect()

    // Run variant calling (using Delly, Manta, and Smoove)
    delly_output = RECODE(DELLY(data))
    manta_output = MANTA(data)
    smoove_output = SMOOVE(data)
    INV_manta = INV(manta_output) 

    // Filter all the VCFs, and then group tuples (by sample ID), then split into 2 channels using multimap (the first channel is the sample ID, and the second channel is the tuple of filtered VCFs for each sample)
    filtered_vcf_output = FILTER(delly_output.mix(manta_output, smoove_output))
    filtered_vcf_output.groupTuple().multiMap{ it -> sample_id: it[0]; vcfs: it[1] }.set{ filtered_vcfs_to_merge }
    merged_vcf_output_raw = MERGE(filtered_vcfs_to_merge.sample_id, filtered_vcfs_to_merge.vcfs)
    merged_vcf_output = NORMALIZATION(merged_vcf_output_raw)
    sample_alignment_svvcf = data.join(merged_vcf_output)
    sample_sv_bed = SV_COVERAGE_EXTRACT(merged_vcf_output)
    sample_alignment_svbed = data.join(sample_sv_bed)
    cov_mosdepth = SV_COVERAGE(data.join(sample_sv_bed))
    GC = CG_CONTENT(sample_sv_bed)
    SV_stats = SV_QC(merged_vcf_output)
    MULTIQC(coverage_reports,SV_stats)
    input_for_features = merged_vcf_output
    .join(GC)
    .join(cov_mosdepth)
    .join(data)
    .map { sample_id, vcf, cg_txt, bed_gz, cram ->
        tuple(
            sample_id,
            cram,     // alignment_file
            vcf,      // merged_vcf
            cg_txt,   // sample_sv_CG
            bed_gz    // sample_coverage
        )
    }

    features_output = EXTRACT_FEATURES(input_for_features)
    features_sv_only = features_output.map { sample_id, features_csv, features_tsv, del, dup, ins, inv ->
    tuple(sample_id, del, dup, ins, inv)}
    vartrustML_o = VARTRUSTML(features_sv_only)

    //expansion_hunter_output = EXPANSION_HUNTER(data)
}
