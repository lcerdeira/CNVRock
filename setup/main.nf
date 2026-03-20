#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

workflow {
    Channel.fromPath([
        "/lustre/scratch127/gsu/malariagen/pipeline_resources/parasite/pf9/resource-bundle/PlasmoDB-54_Pfalciparum3D7_Genome.fasta",
        "/lustre/scratch127/gsu/malariagen/pipeline_resources/parasite/pf9/resource-bundle/PlasmoDB-54_Pfalciparum3D7_Genome.fasta.fai",
        "/lustre/scratch127/gsu/malariagen/pipeline_resources/parasite/pf9/resource-bundle/PlasmoDB-54_Pfalciparum3D7_Genome.dict"
    ])
        .collect()
        .map { files -> [
            reference_fasta: files[0],
            reference_fai:   files[1],
            reference_dict:  files[2]
        ]}
        .set { reference_genome_ch }

    manifest_ch = Channel.fromPath("../assets/paths_to_bams_crams.tsv")

    interval_list_ch = GenerateWholeGenomeIntervalList(reference_genome_ch)

    sample_interval_pairs_ch = manifest_ch
        .splitCsv(sep: "\t")
        .map { row -> [
            sample_id:        row[0],
            path_to_bam_cram: row[1]
        ]}
        .combine(interval_list_ch)
        .combine(reference_genome_ch)

    GenerateReadCountFile(sample_interval_pairs_ch)

    GenerateReadCountFile.out
        .map { it ->
            it.getName().tokenize(".")[0] + "\treadcounts/${it.getName()}"
        }
        .collectFile(
            name:    "paths_to_readcounts.tsv",
            newLine: true,
            sort:    true
        )
}

process GenerateWholeGenomeIntervalList {
    memory "1GB"
    queue  "normal"

    input:
        val(meta)

    output:
        path "whole_genome.interval_list"

    script:
        def java_Xmx = task.memory.mega - 200
        """
        gatk PreprocessIntervals \
            --java-options          "-Xmx${java_Xmx}m" \
            --reference             ${meta.reference_fasta} \
            --bin-length            1000 \
            --padding               0 \
            --interval-merging-rule OVERLAPPING_ONLY \
            --output                whole_genome.interval_list
        """
}

process GenerateReadCountFile {
    memory        { task.attempt == 1 ? "2GB" : "4GB" }
    queue         "normal"
    errorStrategy { task.attempt <= 2 ? "retry" : "terminate" }
    maxRetries    2
    array         100

    publishDir "readcounts/",
        mode: "copy", overwrite: true, failOnError: true

    input:
        tuple val(sample_meta), path(interval_list), val(reference_meta)

    output:
        path "${sample_meta.sample_id}.counts.tsv"

    script:
        def java_Xmx = task.memory.mega - 200
        """
        gatk CollectReadCounts \
            --java-options            "-Xmx${java_Xmx}m" \
            --reference               ${reference_meta.reference_fasta} \
            --intervals               ${interval_list} \
            --input                   ${sample_meta.path_to_bam_cram} \
            --format                  TSV \
            --read-filter             MappingQualityReadFilter \
            --minimum-mapping-quality 40 \
            --interval-merging-rule   OVERLAPPING_ONLY \
            --output                  ${sample_meta.sample_id}.counts.tsv.TEMP

        mv ${sample_meta.sample_id}.counts.tsv.TEMP \
            ${sample_meta.sample_id}.counts.tsv
        """
}
