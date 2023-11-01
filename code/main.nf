nextflow.enable.dsl=1
params.publish_dir_mode = 'copy'
params.outdir = 'results'
params.bwa_meth_index = false
params.save_reference = false
params.title = 'Mapping_run'
Channel
    .fromPath(params.samples)
    .splitCsv(header:true)
    .map{ row -> tuple(row.sampleId, file(row.read1), file(row.read2)) }
    .into { ch_fq_for_fastqc; ch_fq_for_trim_galore }

Channel
    .fromPath(params.fasta)
    .into { ch_fasta_for_faidx; ch_fasta_for_bwa_meth_idx }

if( !params.bwa_meth_index ){
    process MakeBwaMethIndex {
        tag "$fasta"
        publishDir path: "${params.outdir}/reference_genome", saveAs: { params.save_reference ? it : null }, mode: params.publish_dir_mode

        input:
        path fasta from ch_fasta_for_bwa_meth_idx

        output:
        path "${fasta}*" into ch_bwa_meth_indices_for_bwa_meth_align
        path fasta

        script:
        """
        bwameth.py index $fasta
        """
    }
}

process fastqc {
    tag "$sampleId"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: { filename -> filename.indexOf('.zip') > 0 ? "zips/$filename" : "$filename"}

    input:
    tuple val(sampleId), path(read1), path(read2) from ch_fq_for_fastqc

    output:
    path '*.{zip,html}' into ch_fqc_for_mqc

    script:
    """
    fastqc -q -t ${task.cpus} ${read1}
    fastqc -q -t ${task.cpus} ${read2}
    """
}

process trim_galore {
    tag "$sampleId"
    publishDir "${params.outdir}/trim_galore", mode: params.publish_dir_mode,
        saveAs: {filename -> 
            if( filename.indexOf("_fastqc") > 0 ) "FastQC/$filename"
            else if( filename.indexOf("trimming_report.txt" )> 0) "logs/$filename"
            else if( params.save_trimmed) filename
            else null    
        }

    input:
    tuple val(sampleId), path(read1), path(read2) from ch_fq_for_trim_galore

    output:
    tuple val(sampleId), path("*fq.gz") into  ch_trimmed_reads_for_alignment
    path "*trimming_report.txt" into ch_trim_galore_results_for_mqc
    path "*_fastqc.{zip,html}"

    script:
    """
    trim_galore --fastqc --gzip --paired $read1 $read2 --cores ${task.cpus}
    """
}

process AlignWithBwaMeth {
    tag "$sampleId"
    publishDir "${params.outdir}/bwa-mem_alignments", mode: params.publish_dir_mode,
        saveAs: {filename -> 
            if( !params.save_align_intermeds ) filename 
            else if( params.save_align_intermeds ) filename
            else null
        }

    input:
    tuple val(sampleId), path(reads) from ch_trimmed_reads_for_alignment
    path bwa_meth_indices from ch_bwa_meth_indices_for_bwa_meth_align.collect()

    output:
    tuple val(sampleId), path('*.bam') into ch_bam_for_samtools_sort_index_flagstat, ch_bam_for_preseq
    //stdout into results

    script:
    fasta = bwa_meth_indices[0].toString() - '.bwameth' - '.c2t' - '.amb' - '.ann' - '.bwt' - '.pac' - '.sa'
    
    """
    bwameth.py --threads ${task.cpus} --reference $fasta ${reads[0]} ${reads[1]} | samtools view -b - > ${sampleId}.bam
    """
}

process samtools_sort_index_flagstat {
    tag "$sampleId"
    publishDir "${params.outdir}/bwa-mem_alignments", mode: params.publish_dir_mode,
        saveAs: {filename -> 
            if(filename.indexOf("report.txt") >0) "logs/$filename"
            else null
        }
    input:
    tuple val(sampleId), path(bam) from ch_bam_for_samtools_sort_index_flagstat

    output:
    tuple val(sampleId), path("${bam.baseName}.sorted.bam") into ch_bam_for_dedup
    path "${bam.baseName}_flagstat_report.txt" into ch_flagstat_results_for_mqc
    path "${bam.baseName}_stats_report.txt" into ch_samtools_stats_for_mqc

    script:
    """
    samtools sort $bam -@ ${task.cpus} -o ${bam.baseName}.sorted.bam
    samtools flagstat ${bam.baseName}.sorted.bam > ${bam.baseName}_flagstat_report.txt
    samtools stats ${bam.baseName}.sorted.bam > ${bam.baseName}_stats_report.txt
    """
}

process markDuplicates {
    tag "$sampleId"
    publishDir "${params.outdir}/bwa-mem_MarkDuplicates", mode: params.publish_dir_mode,
        saveAs: {filename ->
            filename.indexOf(".bam") == -1 ? "logs/$filename" : "$filename" 
        }

    input:
    tuple val(sampleId), path(bam) from ch_bam_for_dedup

    output:
    tuple val(sampleId), path("${sampleId}.markDups.bam") into ch_bam_dedup_for_qualimap
    tuple val(sampleId), path("${sampleId}.markDups.bam.bai") into ch_bam_index
    path "${sampleId}.markDups_metrics.txt" into ch_markDups_results_for_multiqc

    script:
    if( !task.memory ){
        log.info "[Picard MarkDuplicates] Available memory not known - defaulting to 3GB"
        avail_mem = 3
    } else {
        avail_mem = task.memory.toGiga()
    }
    """
    picard -Xmx${avail_mem}g MarkDuplicates \\
        INPUT=${bam} \\
        OUTPUT=${sampleId}.markDups.bam \\
        METRICS_FILE=${sampleId}.markDups_metrics.txt \\
        REMOVE_DUPLICATES=false \\
        ASSUME_SORTED=true \\
        PROGRAM_RECORD_ID=null \\
        VALIDATION_STRINGENCY=LENIENT
    samtools index -@ ${task.cpus} ${sampleId}.markDups.bam
    """
}

process QualiMap {
    tag "$sampleId"
    publishDir "${params.outdir}/qualimap", mode: params.publish_dir_mode

    input:
    tuple val(sampleId), path(bam) from ch_bam_dedup_for_qualimap

    output:
    path "${sampleId}_qualimap" into ch_qualimap_results_for_multiqc

    script:
    """
    samtools sort ${bam} -@ 4 -o ${sampleId}.sorted4qualimap.bam
    qualimap bamqc -gd HUMAN -bam ${sampleId}.sorted4qualimap.bam -outdir ${sampleId}_qualimap --collect-overlap-pairs -nt ${task.cpus}
    """
}


process MultiQC {
    publishDir "${params.outdir}/MultiQC", mode: params.publish_dir_mode

    input:
    path ('fastqc/*') from ch_fqc_for_mqc.collect().ifEmpty([])
    path ('trim_galore/*') from ch_trim_galore_results_for_mqc.collect().ifEmpty([])
    path ('samtools/*') from ch_flagstat_results_for_mqc.flatten().collect().ifEmpty([])
    path ('samtools/*') from ch_samtools_stats_for_mqc.flatten().collect().ifEmpty([])
    path ('picard/*') from ch_markDups_results_for_multiqc.flatten().collect().ifEmpty([])
    path ('qualimap/*') from ch_qualimap_results_for_multiqc.collect().ifEmpty([])

    output:
    path "*multiqc_report.html" into ch_multiqc_report
    path "*_data"
    //path "multiqc_plots"

    script:
    """
    multiqc -f --title ${params.title} --filename ${params.title}_multiqc_report . -m custom_content -m fastqc -m cutadapt -m picard -m samtools -m qualimap
    """
}
workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
