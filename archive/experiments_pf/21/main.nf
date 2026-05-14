#!/usr/bin/env nextflow
// Experiment 21 — KpSC GATK CollectReadCounts workflow.
// Adapted from data/setup/main.nf; reference paths come from nextflow.config params.

nextflow.enable.dsl = 2

workflow {
    Channel.of([
        params.reference_fasta,
        params.reference_fai,
        params.reference_dict
    ])
        .map { fasta, fai, dict -> [
            reference_fasta: file(fasta),
            reference_fai:   file(fai),
            reference_dict:  file(dict)
        ]}
        .set { reference_genome_ch }

    // TODO: create assets/kpsc-paths-to-bams.tsv with columns: sample_id, path_to_bam
    // TODO: create assets/kpsc-core-genome.bed from Panaroo core genome analysis
    manifest_ch              = Channel.fromPath("../../../assets/kpsc-paths-to-bams.tsv")
    intervals_of_interest_ch = Channel.fromPath("../../../assets/kpsc-core-genome.bed")

    interval_list_ch = GenerateWholeGenomeIntervalList(reference_genome_ch)

    converted_intervals_ch = intervals_of_interest_ch
        .combine(reference_genome_ch)
        .map { bed, ref_meta -> [
            bed:      bed,
            ref_meta: ref_meta
        ]}

    proper_intervals_ch = ConvertBedToIntervalList(converted_intervals_ch)

    combined_intervals_ch = proper_intervals_ch
        .combine(interval_list_ch)
        .map { intervals_file, genome_intervals -> [
            intervals_of_interest:   intervals_file,
            interval_list_of_genome: genome_intervals
        ]}

    restricted_interval_list_ch = RestrictIntervalList(combined_intervals_ch)

    sample_interval_pairs_ch = manifest_ch
        .splitCsv(sep: "\t")
        .map { row -> [
            sample_id:        row[0],
            path_to_bam_cram: row[1]
        ]}
        .combine(restricted_interval_list_ch)
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
            --bin-length            ${params.bin_length} \
            --padding               0 \
            --interval-merging-rule OVERLAPPING_ONLY \
            --output                whole_genome.interval_list
        """
}

process ConvertBedToIntervalList {
    memory "1GB"

    input:
        val(meta)

    output:
        path "intervals_of_interest.interval_list"

    script:
        def java_Xmx = task.memory.mega - 200
        """
        gatk BedToIntervalList \
            --java-options "-Xmx${java_Xmx}m" \
            --INPUT                 ${meta.bed} \
            --OUTPUT                intervals_of_interest.interval_list \
            --SEQUENCE_DICTIONARY   ${meta.ref_meta.reference_dict}
        """
}

process RestrictIntervalList {
    memory "1GB"

    publishDir "../../../assets/",
        mode: "copy", overwrite: true, failOnError: true

    input:
        val(meta)

    output:
        path "kpsc-core-genome.interval_list"

    script:
        def java_Xmx = task.memory.mega - 200
        """
        gatk IntervalListTools \
            --java-options "-Xmx${java_Xmx}m" \
            --INPUT        ${meta.interval_list_of_genome} \
            --SECOND_INPUT ${meta.intervals_of_interest} \
            --ACTION       OVERLAPS \
            --OUTPUT       kpsc-core-genome.interval_list.TEMP

        mv kpsc-core-genome.interval_list.TEMP kpsc-core-genome.interval_list
        """
}

process GenerateReadCountFile {
    memory        { task.attempt == 1 ? "2GB" : "4GB" }
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
            --minimum-mapping-quality ${params.minimum_mapping_quality} \
            --interval-merging-rule   OVERLAPPING_ONLY \
            --output                  ${sample_meta.sample_id}.counts.tsv.TEMP

        mv ${sample_meta.sample_id}.counts.tsv.TEMP \
            ${sample_meta.sample_id}.counts.tsv
        """
}
