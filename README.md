# gatk4-data-processing

### Purpose :
Workflows for processing high-throughput sequencing data for variant discovery with GATK4 and related tools.

### processing-for-variant-discovery-gatk4 :
The processing-for-variant-discovery-gatk4 WDL pipeline implements data pre-processing according to the GATK Best Practices 
(June 2016).  

#### Requirements/expectations 
- Pair-end sequencing data in unmapped BAM (uBAM) format
- One or more read groups, one per uBAM file, all belonging to a single sample (SM)
- Input uBAM files must additionally comply with the following requirements:
  - filenames all have the same suffix (we use ".unmapped.bam")
  - files must pass validation by ValidateSamFile 
  - reads are provided in query-sorted order
  - all reads must have an RG tag

#### Outputs 
- A clean BAM file and its index, suitable for variant discovery analyses.

### Software version requirements :
- GATK 4 or later
- Picard 2.x
- Samtools (see gotc docker)
- Python 2.7

Cromwell version support 
 - Successfully tested on v32
 - Does not work on versions < v23 due to output syntax

Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
