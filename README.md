# gatk4-data-processing

### Purpose :
Workflows for processing high-throughput sequencing data for variant discovery with GATK4 and related tools.

### processing-for-variant-discovery-gatk4 :
The processing-for-variant-discovery-gatk4 WDL pipeline implements data pre-processing according to the [GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035535912). The workflow takes as input an unmapped BAM list file (text file containing paths to unmapped bam files) to perform preprocessing tasks such as mapping, marking duplicates, and base recalibration. It produces a single BAM file and its index suitable for [variant discovery analysis](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-) using tools such as Haplotypecaller.

- If you are starting with FASTQ files visit the [seq-format-conversion](https://github.com/gatk-workflows/seq-format-conversion) repository for workflows to convert FASTQs to unmapped BAMS.
- The processing-for-variant-discovery-gatk4 provides quick and general processing for sequence data using the latest releases of GATK. If users are interested in a more elaborate version of this workflow with quality control tasks and routinely tested for validity (useful in production environments) then visit the [gatk4-genome-processing-pipeline](https://github.com/gatk-workflows/gatk4-genome-processing-pipeline) repository.
- The BAM output from processing-for-variant-discovery-gatk4 can be used to perform a variety of other analysis like somatic short variant discovery, germline short variant discovery, or germline copy number variant discovery. Visit the GATK Best Practices documentation to determine what to do next with the BAM files.

#### Requirements/expectations:
- Pair-end sequencing data in unmapped BAM (uBAM) format
- One or more read groups, one per uBAM file, all belonging to a single sample (SM)
- Input uBAM files must additionally comply with the following requirements:
  - filenames all have the same suffix (we use ".unmapped.bam")
  - files must pass validation by ValidateSamFile 
  - reads are provided in query-sorted order
  - all reads must have an RG tag
- Reference index files must be in the same directory as source (e.g. reference.fasta.fai in the same directory as reference.fasta)

#### Outputs: 
- A clean BAM file and its index, suitable for variant discovery analyses.

### Software version requirements :
- GATK 4 or later
- BWA 0.7.15-r1140
- Picard 2.16.0-SNAPSHOT
- Samtools 1.3.1 (using htslib 1.3.1)
- Python 2.7
- Cromwell version support 
  - Successfully tested on v59 
  
### Important Notes :
- Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
- For help running workflows on the Google Cloud Platform or locally please
view the following tutorial [(How to) Execute Workflows from the gatk-workflows Git Organization](https://gatk.broadinstitute.org/hc/en-us/articles/360035530952).
- Please visit the [User Guide](https://gatk.broadinstitute.org/hc/en-us/categories/360002310591) site for further documentation on our workflows and tools.
- Relevant reference and resources bundles can be accessed in [Resource Bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360036212652).

### Contact Us :
- The following material is provided by the Data Science Platforum group at the Broad Institute. Please direct any questions or concerns to one of our forum sites : [GATK](https://gatk.broadinstitute.org/hc/en-us/community/topics) or [Terra](https://support.terra.bio/hc/en-us/community/topics/360000500432).


### LICENSING :
Copyright Broad Institute, 2021 | BSD-3
This script is released under the WDL open source code license (BSD-3) (full license text at https://github.com/openwdl/wdl/blob/master/LICENSE). Note however that the programs it calls may be subject to different licenses. Users are responsible for checking that they are authorized to run all programs before running this script.
- [GATK](https://software.broadinstitute.org/gatk/download/licensing.php)
- [BWA](http://bio-bwa.sourceforge.net/bwa.shtml#13)
- [Picard](https://broadinstitute.github.io/picard/)
- [Samtools](http://www.htslib.org/terms/)
