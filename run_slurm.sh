#!/bin/bash
#SBATCH --job-name=NextflowSV
#SBATCH --partition=prod
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=nextflowSV-%j.out
#SBATCH --error=nextflowSV-%j.err
#SBATCH --nodes=1

source /shared/conda/miniconda3/bin/activate
conda activate nextflow


export TMPDIR=/shared/work/PI-tommaso.pippucci/paradigm/tmp
export NXF_TEMP=/shared/work/PI-tommaso.pippucci/paradigm/tmp
export SINGULARITY_TMPDIR=/shared/work/PI-tommaso.pippucci/paradigm/tmp
export SINGULARITY_CACHEDIR=/shared/work/PI-tommaso.pippucci/paradigm/tmp
export NXF_SINGULARITY_CACHEDIR=/shared/work/PI-tommaso.pippucci/paradigm/tmp

date
nextflow run /shared/work/PI-tommaso.pippucci/RF-WGS/nextflow_sv_pipeline/sv_main.nf  --reference_fasta /shared/archive/ngsbo/migrated-from-ngsra/db/trioCEU_1KGP_resources/GRCh38_full_analysis_set_plus_decoy_hla.fa --sample_alignments_tsv /shared/work/PI-tommaso.pippucci/RF-WGS/nextflow_sv_pipeline/test_1Kg_4.2.7.csv --outdir /shared/work/PI-tommaso.pippucci/RF-WGS/nextflow_sv_pipeline/results/ --delly_exclude_regions_bed /shared/work/PI-tommaso.pippucci/RF-WGS/nextflow_sv_pipeline/Hg38/exlude.regions.delly.human.hg38.excl.tsv --smoove_exclude_regions_bed /shared/work/PI-tommaso.pippucci/RF-WGS/nextflow_sv_pipeline/Hg38/exclude.smoove.cnvnator_100bp.GRCh38.20170403.bed --manta_inv /shared/work/PI-tommaso.pippucci/RF-WGS/nextflow_sv_pipeline/Hg38/convertInversion.py --expansion_hunter_variant_catalog_json /shared/work/PI-tommaso.pippucci/RF-WGS/nextflow_sv_pipeline/Hg38/variant_catalog_genes.json --recode_delly_python_script /shared/work/PI-tommaso.pippucci/RF-WGS/nextflow_sv_pipeline/Hg38/recode_delly.py --mosdepth_segmental_duplications_bed /shared/work/PI-tommaso.pippucci/RF-WGS/nextflow_sv_pipeline/Hg38/Segmental_dups_hg38_frt_srt.bed --samtools_path /shared/home/emanuela.iovino/.conda/envs/mysamtools/bin/samtools --seg_duplications Hg38/segmental_duplications.bed --repeat_masker Hg38/repeat_masker.tsv --singularity_cache /shared/work/PI-tommaso.pippucci/RF-WGS/my_singularity_container/ --bind_path  /shared/archive/ngsbo/migrated-from-ngsra/cram/genome_iit/,/shared/archive01/PI-tommaso.pippucci/1KG_samples/1000genomes-dragen-v4-2-7/,/shared/work/PI-tommaso.pippucci/RF-WGS/nextflow_sv_pipeline/paper_files,/shared/archive/ngsbo/migrated-from-ngsra/db/,/shared/work/PI-tommaso.pippucci/paradigm/tmp,/shared/home/emanuela.iovino/.conda/envs/mysamtools/bin/samtools,/shared/work/PI-tommaso.pippucci/paradigm/,/shared/work/PI-tommaso.pippucci/RF-WGS/my_singularity_container/,/shared/work/PI-tommaso.pippucci/RF-WGS/nextflow_sv_pipeline/Hg38/,/shared/work/PI-tommaso.pippucci/RF-WGS/brusco/results/,/shared/work/PI-tommaso.pippucci/RF-WGS/nextflow_sv_pipeline/tmp_nx/,/shared/archive/ngsbo/migrated-from-ngsra/cram/genome_iit/,/shared/work/PI-tommaso.pippucci/RF-WGS/benchmark_sv/files_scr_paper -profile slurm -resume  --tmp_dir /shared/work/PI-tommaso.pippucci/RF-WGS/nextflow_sv_pipeline/tmp_nx

date

