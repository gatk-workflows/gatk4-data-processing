# gatk4-data-processing

### Purpose :
Workflows for processing high-throughput sequencing data for variant discovery with GATK4 and related tools.

### processing-for-variant-discovery-gatk4 :
The processing-for-variant-discovery-gatk4 WDL pipeline implements data pre-processing according to the GATK Best Practices.  

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
  - Successfully tested on v37 
  - Does not work on versions < v23 due to output syntax
  
### Important Note :
- The provided JSON is meant to be a ready to use example JSON template of the workflow. It is the user’s responsibility to correctly set the reference and resource input variables using the [GATK Tool and Tutorial Documentations](https://software.broadinstitute.org/gatk/documentation/).
- Relevant reference and resources bundles can be accessed in [Resource Bundle](https://software.broadinstitute.org/gatk/download/bundle).
- Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
- For help running workflows on the Google Cloud Platform or locally please
view the following tutorial [(How to) Execute Workflows from the gatk-workflows Git Organization](https://software.broadinstitute.org/gatk/documentation/article?id=12521).
- The following material is provided by the GATK Team. Please post any questions or concerns to one of our forum sites : [GATK](https://gatkforums.broadinstitute.org/gatk/categories/ask-the-team/) , [FireCloud](https://gatkforums.broadinstitute.org/firecloud/categories/ask-the-firecloud-team) or [Terra](https://broadinstitute.zendesk.com/hc/en-us/community/topics/360000500432-General-Discussion) , [WDL/Cromwell](https://gatkforums.broadinstitute.org/wdl/categories/ask-the-wdl-team).
- Please visit the [User Guide](https://software.broadinstitute.org/gatk/documentation/) site for further documentation on our workflows and tools.

### LICENSING :
Copyright Broad Institute, 2019 | BSD-3
This script is released under the WDL open source code license (BSD-3) (full license text at https://github.com/openwdl/wdl/blob/master/LICENSE). Note however that the programs it calls may be subject to different licenses. Users are responsible for checking that they are authorized to run all programs before running this script.
- [GATK](https://software.broadinstitute.org/gatk/download/licensing.php)
- [BWA](http://bio-bwa.sourceforge.net/bwa.shtml#13)
- [Picard](https://broadinstitute.github.io/picard/)
- [Samtools](http://www.htslib.org/terms/)
