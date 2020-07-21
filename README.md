# bulk_RNASeq
methods used for bulk RNA-seq processing and analysis

### library preparation methods

[SMARTer Stranded Total RNA Sample Prep Kit](https://www.takarabio.com/assets/documents/User%20Manual/SMARTer%20Stranded%20Total%20RNA%20Sample%20Prep%20Kit%20-%20HI%20Mammalian%20User%20Manual_092617.pdf) - HI Mammalian (Takara Bio USA, Inc., Mountain View, CA, USA) incorporates both RiboGone and SMART (Switching Mechanism At 5’ end of RNA Template) technologies to deplete nuclear rRNA and synthesize first-strand cDNA. This along with PCR amplification and AMPure Bead Purification generates Illumina-compatible libraries.

Sequencing Analysis Guidelines
IMPORTANT: The first three nucleotides of the first sequencing read (Read 1) are derived from the template-switching oligo. These three nucleotides must be trimmed prior to mapping. Read 1 is derived from the sense strand of the input RNA. If you are performing paired-end sequencing, Read 2 will correspond to the antisense strand.

Data were sequenced at the Translational Genomics Research Institute on an Illumina NovaSeq6000.

### transcriptome annotation

Custom annotation combining GENCODE 29 and lncipedia 5.2. More details are available in the [corresponding GitHub repo](https://github.com/FOUNDINPD/annotation-RNA).

### packages and versions used
```
  perl: "perl/5.26.2"
  samtools: "samtools/1.9"
  bwa: "bwa/0.7.17"
  ciri: "ciri/2.0.6"
  R: "R/3.4.1"
  rseqc: "RSeQC/2.6.4"
  python2: "python/2.7.13"
  star: "STAR/2.6.1d"
  fastqc: "fastqc/0.11.8"
  featureCounts: "subread/1.6.4"
  multiqc: "multiqc/1.8"
```

### cutadapt (v2.7)
[cutadapt documentation](https://cutadapt.readthedocs.io/en/stable/)

```
cutadapt -u 3 -o ${sampleID}.R1.fastq.gz -p ${sampleID}.R2.fastq.gz ${input_dir}/*${sampleID}*R1*.fastq.gz ${input_dir}/*${sampleID}*R2*.fastq.gz
```

### STAR (v2.6.1d)
[STAR documentation](https://github.com/alexdobin/STAR)
```
    STAR \
    --genomeDir {{ config.file_locations.starIndex }} \
    --sjdbGTFfile {{ config.file_locations.annotation }} \
    --readFilesCommand zcat \
    --readFilesIn {{ fastq_sets | map(attribute='fq1') | join(',') }} {{ fastq_sets | map(attribute='fq2') | join(',') }} \
    --runThreadN ${SLURM_CPUS_PER_TASK} \
    --limitBAMsortRAM 9600000000 \
    --outFileNamePrefix {{ sample }}. \
    --runMode alignReads \
    --outSAMtype BAM Unsorted \
    --outSAMmode Full \
    --outSAMstrandField intronMotif \
    --outFilterType BySJout \
    --outSAMunmapped Within \
    --outSAMmapqUnique 255 \
    --outFilterMultimapNmax 20 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverLmax 0.1 \
    --alignMatesGapMax 1000000 \
    --seedSearchStartLmax 50 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignSJoverhangMin 18 \
    --alignSJDBoverhangMin 18 \
    --chimSegmentMin 18 \
    --chimJunctionOverhangMin 18 \
    --outSJfilterOverhangMin 18 18 18 18 \
    --alignTranscriptsPerReadNmax 50000 \
    --genomeLoad NoSharedMemory
```

### salmon (v1.2.1)
[salmon docs](https://combine-lab.github.io/salmon/)
```
    salmon quant \
      --index {{ config.file_locations.salmonIndex }} \
      --libType A \
      --output {{ sample }} \
      --threads ${SLURM_CPUS_PER_TASK} \
      --numBootstraps 100 \
      --validateMappings \
      --seqBias \
      --gcBias \
      --dumpEq \
      --mates1 {{ fastq_sets | map(attribute='fq1') | join(' ') }} \
      --mates2 {{ fastq_sets | map(attribute='fq2') | join(' ') }} \
      --geneMap {{ config.file_locations.annotation }}
```

### featureCounts (v1.6.4)
[subread/featureCounts docs](http://subread.sourceforge.net/)

```
    featureCounts -p -T ${SLURM_CPUS_PER_TASK} -t exon -g gene_id -a {{ config.file_locations.annotation }} -s 1 -o {{ sample }}.txt {{ config.workdir }}/aligned/{{ sample }}/{{ sample }}.Aligned.out.bam
```

### samtools (v1.9) 
[samtools docs](http://www.htslib.org/)

```
    samtools sort -@ ${SLURM_CPUS_PER_TASK} -m 4G -O bam -o {{ config.workdir }}/aligned/{{ sample }}/{{ sample }}.Aligned.out.sorted.bam {{ config.workdir }}/aligned/{{ sample }}/{{ sample }}.Aligned.out.bam
    
    samtools index {{ config.workdir }}/aligned/{{ sample }}/{{ sample }}.Aligned.out.sorted.bam
```

### CIRI2 (v2.0.6)
[CIRI docs](https://sourceforge.net/projects/ciri/files/CIRI2/)

```
    bwa mem -t ${SLURM_CPUS_PER_TASK} -T 19 {{ config.file_locations.bwaIndex }} <(zcat {{ fastq_sets | map(attribute='fq1') | join(' ') }}) <(zcat {{ fastq_sets | map(attribute='fq2') | join(' ') }}) > {{ sample }}.sam

    perl ciri2.pl -I {{ sample }}.sam -O {{ sample }}.ciri -F {{ config.file_locations.genome }} -A {{ config.file_locations.annotation }} -T ${SLURM_CPUS_PER_TASK}
```
