	#!/usr/bin/env python
"""
RNA TopHat pipeline
documentation.
"""

from argparse import ArgumentParser
import os
import os.path
import sys
from subprocess import call
import subprocess
import re
import pypiper
from datetime import datetime


# Argument Parsing
# #######################################################################################
parser = ArgumentParser(description='Pypiper arguments.')

parser = pypiper.add_pypiper_args(parser, all_args=True)


parser.add_argument('-f', dest='filter', action='store_false', default=True)
parser.add_argument('-d', dest='markDupl', action='store_true', default=False)
parser.add_argument('-w', '--wigsum', default=500000000, dest='wigsum', type=int, help='Target wigsum for track normalisation')


# Core-seq as optional parameter
parser.add_argument('-cs', '--core-seq', default=False, dest='coreseq', action='store_true', help='CORE-seq Mode')

args = parser.parse_args()

if args.single_or_paired == "paired":
	args.paired_end = True
else:
	args.paired_end = False


if not args.input:
	parser.print_help() #or, print_usage() for less verbosity
	raise SystemExit

# Merging
########################################################################################
# If 2 unmapped bam files are given, then these are to be merged.
# Must be done here to initialize the sample name correctly

# This is now deprecated (there is no default sample name implemented)
#merge = False
#if (len(args.input) > 1):
#	merge = True
#	if (args.sample_name == "default"):
#		args.sample_name = "merged";
#else:
#	if (args.sample_name == "default"):
#		args.sample_name = os.path.splitext(os.path.basename(args.input[0]))[0]

# Set up environment path variables
########################################################################################
# Set up an container class to hold paths
class Container:
	pass

paths = Container()
paths.scripts_dir = os.path.dirname(os.path.realpath(__file__))

# import the yaml config (work in progress)...
import yaml
pipelines_config_file = os.path.join(paths.scripts_dir, "rnaTopHat.yaml")
config = yaml.load(open(pipelines_config_file, 'r'))

# Resources
paths.resources_dir = config["resources"]["resources"]
paths.adapter_file = config["resources"]["adapters"]
paths.ref_genome = config["resources"]["genomes"]
paths.ref_genome_fasta = os.path.join(paths.ref_genome, args.genome_assembly, args.genome_assembly + ".fa")
paths.chrom_sizes = os.path.join(paths.ref_genome, args.genome_assembly, args.genome_assembly + ".chromSizes")
paths.bowtie_indexed_genome = os.path.join(paths.ref_genome, args.genome_assembly, "indexed_bowtie2", args.genome_assembly)

# tophat specific resources
paths.gtf = os.path.join(paths.ref_genome, args.genome_assembly , "ucsc_" + args.genome_assembly + "_ensembl_genes.gtf")
paths.gene_model_bed = os.path.join(paths.ref_genome, args.genome_assembly , "ucsc_" + args.genome_assembly + "_ensembl_genes.bed")
paths.gene_model_sub_bed = os.path.join(paths.ref_genome, args.genome_assembly , "ucsc_" + args.genome_assembly + "_ensembl_genes_500rand.bed")


# Tools
paths.trimmomatic_jar = config["tools"]["trimmomatic_epignome"]
paths.bowtie2 = config["tools"]["bowtie2"]
paths.picard_jar = config["tools"]["picard"]
paths.bed2bigBed = config["tools"]["bed2bigBed"]
paths.bed2bigWig = config["tools"]["bed2bigWig"]

paths.tophat = config["tools"]["tophat2"]
paths.bam2wig = config["tools"]["bam2wig"]
paths.wigToBigWig = os.path.join(paths.resources_dir, "tools", "wigToBigWig")
paths.read_distribution = config["tools"]["read_distribution"]
paths.gene_coverage = config["tools"]["gene_coverage"]

# Output
paths.pipeline_outfolder = os.path.join(args.output_parent, args.sample_name + "/")

# Initialize
pm = pypiper.PipelineManager(name="rnaTopHat", outfolder=paths.pipeline_outfolder, args=args)

ngstk = pypiper.NGSTk(pm=pm)

raw_folder = os.path.join(paths.pipeline_outfolder, "raw/")
fastq_folder = os.path.join(paths.pipeline_outfolder, "fastq/")

# Merge/Link sample input and Fastq conversion
# These commands merge (if multiple) or link (if single) input files,
# then convert (if necessary, for bam, fastq, or gz format) files to fastq.
################################################################################
pm.timestamp("### Merge/link and fastq conversion: ")

local_input_files = ngstk.merge_or_link([args.input, args.input2], raw_folder, args.sample_name)
cmd, out_fastq_pre, unaligned_fastq = ngstk.input_to_fastq(local_input_files, args.sample_name, args.paired_end, fastq_folder)
pm.run(cmd, unaligned_fastq, 
	follow=ngstk.check_fastq(local_input_files, unaligned_fastq, args.paired_end))
pm.clean_add(out_fastq_pre + "*.fastq", conditional=True)

pm.report_result("File_mb", ngstk.get_file_size(local_input_files))
pm.report_result("Read_type", args.single_or_paired)
pm.report_result("Genome", args.genome_assembly)

# Adapter trimming
################################################################################
pm.timestamp("### Adapter trimming: ")

cmd = "java -Xmx4g -jar "+ paths.trimmomatic_jar

if not args.paired_end:
	cmd += " SE -phred33 -threads 30"
	cmd += " -trimlog " + os.path.join(paths.pipeline_outfolder, "fastq/") + "trimlog.log "
	cmd += out_fastq_pre + "_R1.fastq "
	cmd += out_fastq_pre + "_R1_trimmed.fastq "

else:
	cmd += " PE -phred33 -threads 30"
	cmd += " -trimlog " + os.path.join(paths.pipeline_outfolder, "fastq/") + "trimlog.log "
	cmd += out_fastq_pre + "_R1.fastq "
	cmd += out_fastq_pre + "_R2.fastq "
	cmd += out_fastq_pre + "_R1_trimmed.fastq "
	cmd += out_fastq_pre + "_R1_unpaired.fastq "
	cmd += out_fastq_pre + "_R2_trimmed.fastq "
	cmd += out_fastq_pre + "_R2_unpaired.fastq "

# for Core-seq, trim off the first 6bp and the bit adjacent to identified adapter sequences:
if args.coreseq:
	cmd += " HEADCROP:6 ILLUMINACLIP:" + paths.adapter_file + ":2:10:4:1:true:epignome:5 SLIDINGWINDOW:4:1 MAXINFO:16:0.40 MINLEN:25"
# otherwise just look for normal adapters:
else:
	cmd += " ILLUMINACLIP:" + paths.adapter_file + ":2:10:4:1:true SLIDINGWINDOW:4:1 MAXINFO:16:0.40 MINLEN:21"

trimmed_fastq = out_fastq_pre + "_R1_trimmed.fastq"
trimmed_fastq_R2 = out_fastq_pre + "_R2_trimmed.fastq"
#pm.run(cmd, out_fastq_pre + "_R1_trimmed.fastq")
#pm.report_result("Trimmed_reads", ngstk.count_reads(trimmed_fastq,args.paired_end))

pm.run(cmd, trimmed_fastq, 
	follow = ngstk.check_trim(trimmed_fastq, trimmed_fastq_R2, args.paired_end,
		fastqc_folder = os.path.join(paths.pipeline_outfolder, "fastqc/")))




# RNA Tophat pipeline.
########################################################################################
pm.timestamp("### TopHat alignment: ")
tophat_folder = paths.pipeline_outfolder + "/tophat_" + args.genome_assembly + "/"
pm.make_sure_path_exists(tophat_folder)
out_tophat = tophat_folder + args.sample_name + ".aln.bam"

align_paired_as_single = True # FH: this appears to be the default behavior of the pipeline at the moment. Should that be configurable by args?

if not args.paired_end:
	cmd = paths.tophat
	cmd += " --GTF " + paths.gtf
	cmd += " --b2-L 15 --library-type fr-unstranded --mate-inner-dist 150 --max-multihits 100 --no-coverage-search --num-threads 2"
	cmd += " --output-dir " + tophat_folder
	cmd += " " + paths.bowtie_indexed_genome
	cmd += " " + out_fastq_pre + "_R1_trimmed.fastq"

else:
	cmd = paths.tophat
	cmd += " --GTF " + paths.gtf
	cmd += " --b2-L 15 --library-type fr-unstranded --mate-inner-dist 150 --max-multihits 100 --no-coverage-search --num-threads 2"
	cmd += " --output-dir " + tophat_folder
	cmd += " " + paths.bowtie_indexed_genome
	# FH: if you use this code, you align both mates separately. As a result, the count_unique_mapped_reads method in paired-end mode will return 0, because the mate flags are not set
	if align_paired_as_single:
		cmd += " " + out_fastq_pre + "_R1_trimmed.fastq" + "," + out_fastq_pre + "_R2_trimmed.fastq"
	else:
		cmd += " " + out_fastq_pre + "_R1_trimmed.fastq"
		cmd += " " + out_fastq_pre + "_R2_trimmed.fastq"

pm.run(cmd, tophat_folder + "/align_summary.txt", shell=False)

pm.timestamp("### renaming tophat aligned bam file ")
cmd = "mv " + tophat_folder + "/accepted_hits.bam " + out_tophat
pm.run(cmd, out_tophat, shell=False, follow=lambda:
	pm.report_result("Aligned_reads", ngstk.count_unique_mapped_reads(out_tophat,args.paired_end and not align_paired_as_single)))





pm.timestamp("### BAM to SAM sorting and indexing: ")
if args.filter:
	cmd = ngstk.bam_conversions(out_tophat,False)
	pm.run(cmd,  re.sub(".bam$", "_sorted.bam", out_tophat) ,shell=True)
else:
	cmd = ngstk.bam_conversions(out_tophat,True)
	pm.run(cmd, re.sub(".bam$", "_sorted.depth", out_tophat),shell=True)

pm.clean_add(out_tophat, conditional=False)
pm.clean_add(re.sub(".bam$" , ".sam", out_tophat), conditional=False)



if not args.filter and args.markDupl:
	pm.timestamp("### MarkDuplicates: ")

	aligned_file = re.sub(".sam$", "_sorted.bam",  out_tophat)
	out_file = re.sub(".sam$", "_dedup.bam", out_tophat)
	metrics_file = re.sub(".sam$", "_dedup.metrics", out_tophat)
	cmd = ngstk.markDuplicates(aligned_file, out_file, metrics_file)
	pm.run(cmd, out_file, follow= lambda:
		pm.report_result("Deduplicated_reads", ngstk.count_unique_mapped_reads(out_file, args.paired_end and not align_paired_as_single)))




#read filtering
########################################################################################

if args.filter:
	pm.timestamp("### Aligned read filtering: ")

	out_sam_filter = tophat_folder + args.sample_name + ".aln.filt.sam"
	headerLines = subprocess.check_output("samtools view -SH " + re.sub(".bam$", ".sam", out_tophat) + "|wc -l", shell=True).strip()
	cmd = "python " + paths.scripts_dir + "/bisulfiteReadFiltering_forRNA.py"
	cmd += " --infile=" + re.sub(".bam$",".sam", out_tophat)
	cmd += " --outfile=" + out_sam_filter
	cmd += " --skipped=" + out_sam_filter.replace(".filt." , ".skipped.")
	cmd += " --skipHeaderLines=" + headerLines
	cmd += " --genome=" + args.genome_assembly
	cmd += " --genomeDir=" + paths.ref_genome
	cmd += " --minNonCpgSites=3"
	cmd += " --minConversionRate=0.9"
	cmd += " --maxConversionRate=0.1"
	cmd += " -r"


	if args.paired_end and not align_paired_as_single:
		cmd = cmd + " --pairedEnd"

	pm.run(cmd, out_sam_filter, follow=lambda:
		pm.report_result("Filtered_reads", ngstk.count_unique_mapped_reads(out_sam_filter, args.paired_end and not align_paired_as_single)))



#	pm.timestamp("### Set flag in skipped reads to 4 (unmapped): ")
#	joined_sam = out_sam_filter.replace(".filt." , ".filtFlag.")
#	skipped_sam = out_sam_filter.replace(".filt." , ".skipped.")
#	cmd = "samtools view -hS " + out_sam_filter + " > " + joined_sam + "\n"
#	cmd += "awk -v skip=" + headerLines + " -v OFS=\"\\t\" '{if (NR>skip){$2=4;$3=\"*\";$4=0;$5=0;$6=\"*\";$7=\"*\";$8=0;$9=0; print}}' " + skipped_sam
#	cmd += " >>" + joined_sam
#
#	pm.call_lock(cmd, joined_sam , shell=True)
#
#	if not args.no_check:
#		x = ngstk.count_reads(joined_sam, args.paired_end and not align_paired_as_single)
#		pm.report_result("Joined_reads", x)
#
#	pm.timestamp("### Filtered: SAM to BAM conversion, sorting and depth calculation: ")
#	cmd = ngstk.sam_conversions(joined_sam)
#	pm.call_lock(cmd, joined_sam.replace(".sam" , "_sorted.depth"),shell=True)
#
#	pm.timestamp("### Skipped: SAM to BAM conversion and sorting: ")
#	cmd = ngstk.sam_conversions(skipped_sam,False)
#	pm.call_lock(cmd, skipped_sam.replace(".sam" , "_sorted.bam"),shell=True)
#
#	if args.markDupl:
#		pm.timestamp("### MarkDuplicates: ")
#
#		aligned_file = joined_sam.replace(".sam" , "_sorted.bam")
#		out_file = joined_sam.replace(".sam" , "_dedup.bam")
#		metrics_file = joined_sam.replace(".sam" , "_dedup.metrics")
#		cmd = ngstk.markDuplicates(paths, aligned_file, out_file, metrics_file)
#		pm.call_lock(cmd, out_file)
#
#		if not args.no_check:
#			x = ngstk.count_mapped_reads(out_file)
#			pm.report_result("Deduplicated_reads", x)

	pm.timestamp("### Filtered: SAM to BAM conversion, sorting and depth calculation: ")
	cmd = ngstk.sam_conversions(out_sam_filter)
	pm.run(cmd, re.sub(".sam$", "_sorted.depth", out_sam_filter),shell=True)

	skipped_sam = out_sam_filter.replace(".filt." , ".skipped.")
	pm.timestamp("### Skipped: SAM to BAM conversion and sorting: ")
	cmd = ngstk.sam_conversions(skipped_sam,False)
	pm.run(cmd, re.sub(".sam$" , "_sorted.bam", skipped_sam),shell=True)

	pm.clean_add(skipped_sam, conditional=False)
	pm.clean_add(re.sub(".sam$" , ".bam", skipped_sam), conditional=False)
	pm.clean_add(out_sam_filter, conditional=False)
	pm.clean_add(re.sub(".sam$" , ".bam", out_sam_filter), conditional=False)



#create tracks
########################################################################################
pm.timestamp("### bam2wig: ")
if args.filter:
	trackFile = re.sub(".sam$", "_sorted.bam", out_sam_filter)
	cmd = paths.bam2wig + " -i" + trackFile
	cmd += " -s " + paths.chrom_sizes
	cmd += " -o " + re.sub(".sam$" , "_sorted", out_sam_filter)
	cmd += " -t " + str(args.wigsum)
	pm.run(cmd, re.sub(".sam$" , "_sorted.wig",out_sam_filter),shell=False)

	pm.timestamp("### wigToBigWig: ")
	cmd = paths.wigToBigWig + " " + re.sub(".sam$" , "_sorted.wig",out_sam_filter)
	cmd += " " + paths.chrom_sizes
	cmd += " " + re.sub(".sam$" , "_sorted.bw",out_sam_filter)
	pm.run(cmd, re.sub(".sam$" , "_sorted.bw",out_sam_filter),shell=False)

else:
	trackFile = re.sub(".bam$", "_sorted.bam",out_tophat)
	cmd = paths.bam2wig + " -i" + trackFile
	cmd += " -s " + paths.chrom_sizes
	cmd += " -o " + re.sub(".bam$" , "_sorted",out_tophat)
	cmd += " -t " + str(args.wigsum)
	pm.run(cmd, re.sub(".bam$" , "_sorted.wig",out_tophat),shell=False)

	pm.timestamp("### wigToBigWig: ")
	cmd = paths.wigToBigWig + " " + re.sub(".bam$" , "_sorted.wig",out_tophat)
	cmd += " " + paths.chrom_sizes
	cmd += " " + re.sub(".bam$" , "_sorted.bw",out_tophat)
	pm.run(cmd, re.sub(".bam$" , "_sorted.bw", out_tophat),shell=False)


pm.timestamp("### read_distribution: ")
cmd = paths.read_distribution + " -i " + trackFile
cmd	+= " -r " + paths.gene_model_bed
cmd += " > " + re.sub("_sorted.bam$", "_read_distribution.txt",trackFile)
pm.run(cmd, re.sub("_sorted.bam$", "_read_distribution.txt",trackFile), shell = True,
	nofail = True)


pm.timestamp("### gene_coverage: ")
cmd = paths.gene_coverage + " -i " + re.sub(".bam$" , ".bw",trackFile)
cmd	+= " -r " + paths.gene_model_sub_bed
cmd += " -o " + re.sub("_sorted.bam$", "",trackFile)
pm.run(cmd, re.sub("_sorted.bam$", ".geneBodyCoverage.png",trackFile), shell = False)


# Cleanup
########################################################################################


# remove temporary marker file:
pm.stop_pipeline()
