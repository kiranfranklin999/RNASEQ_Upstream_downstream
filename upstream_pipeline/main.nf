// main.nf

params.configName = 'nextflow.config'

// Process: FastQC
process fastqc {
  input:
  set pairId, file(reads) from fastq

  output:
  file "${pairId}_fastqc.zip" into qc


  script:
  """
  fastqc -o qc ${reads}
  """
}

// Process: Trimomatic
process trim {
  input:
  set pairId, file(reads) from fastqc.out

  output:
  file "${pairId}_trimmed.fastq" into trim_qc


  script:
  """
  trimomatic PE -phred33 ${reads} ${pairId}_trimmed.fastq ${pairId}_unpaired.fastq \
    ${pairId}_trimmed_paired.fastq ${pairId}_unpaired_paired.fastq \
    ILLUMINACLIP:${params.trimAdapter} \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
  """
}

// Process: STAR Alignment
process star {
  input:
  set pairId, file(reads) from trim_qc.out

  output:
  file "${pairId}_Aligned.sortedByCoord.out.bam" into align_qc


  script:
  """
  STAR --genomeDir ${params.alignGenomeDir} \
    --readFilesIn ${reads} \
    --runThreadN 8 \
    --sjdbGTFfile ${params.alignGTF} \
    --outFileNamePrefix ${pairId}_ \
    --outSAMtype BAM SortedByCoordinate
  """
}

// Process: FeatureCounts
process featureCounts {
  input:
  file bam from align_qc.out

  output:
  file "featureCounts.txt" into results


  script:
  """
  featureCounts -T 8 -p -t exon -g gene_id -a ${params.alignGTF} -o featureCounts.txt ${bam}
  """
}

// Process: MultiQC
process multiqc {
  input:
  file bam from align_qc.out

  output:
  file "multiqc_report.html" into results


  script:
  """
  multiqc --outdir ${params.multiqcDir} .
  """
}

workflow {
  Channel.fromFilePairs('data/input/*_R{1,2}.fastq')
    .set{ fastq }

  fastqc(fastq)

  trim(fastqc.out)

  star(trim_qc.out)

  featureCounts(align_qc.out)

  multiqc(align_qc.out)
}
