<p align="center"><h1 align="center">HiLOW (HiC and HiChIP Looping Organization Workflow)</h1></p>

## SYNTAX to execute HiLOW

```java 
java -jar <cromwell.jar> run hilow.wdl --input <inputs.json> --options <options.json>
```


### Required files (included in the <inputs.json>)

1. Genome Fragments file: Restriction fragment file in BED format
    ```bash
    # Arima has multiple restriction enzymes. The one provided in `example/` was generated for hg19.
    # Here is the code to create genome fragments for Arima using our docker image.  

    singularity exec docker://ghcr.io/stjude/abralab/hilow:v1.0 /HiC-Pro_3.1.0/bin/utils/digest_genome.py -r ^GATC G^ANTC -o <output_fragments_name> <genome_fasta>
    ```
2. Fastqfiles: One or More Sample FASTQ files (specified as `fastqfiles_R1` and `fastqfiles_R2`).
3. Genome: Reference FASTA file. Or bowtie2 index and chromosome sizes need to be provided in place of Reference FASTA.
4. Genome Name: Reference Genome Name (eg. hg19, hg38, mm10).
5. Ligation site: Restriction enzyme ligation site used by HiC-Pro.
6. Specify if HiChIP or HiC (i.e. True if HiChIP or False if HiC. Default is True).

#### Example of input json options
```json
{
    "hilow.fastqfiles_R1" : 
    [
        "./fastq/CTCF_HiChIP_Rep1_R1.fastq.gz"
    ],
    "hilow.fastqfiles_R2" : 
    [
        "./fastq/CTCF_HiChIP_Rep1_R2.fastq.gz"
    ],
    "hilow.sampleid" : "TESTING345",
    "hilow.ligationsite" : "GATCGATC",
    "hilow.blacklist" : "./hg19-blacklist.v2.bed",
    "hilow.genomefragment" : "./hg19_MboI_resfrag.bed.gz",
    "hilow.genomename" : "hg19",
    "hilow.chromsizes" : "./hg19.v2.size",
    "hilow.bowtie2index" :
    [
        "./bowtie2index/hg19.1.bt2",
        "./bowtie2index/hg19.2.bt2",
        "./bowtie2index/hg19.3.bt2",
        "./bowtie2index/hg19.4.bt2",
        "./bowtie2index/hg19.rev.1.bt2",
        "./bowtie2index/hg19.rev.2.bt2"
    ]
}
```

Please consult `TestData/bash-hilow.sh` for example usage and execution
