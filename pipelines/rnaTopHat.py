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
merge = False
if (len(args.input) > 1):
	merge = True
	if (args.sample_name == "default"):
		args.sample_name = "merged";
else:
	if (args.sample_name == "default"):
		args.sample_name = os.path.splitext(os.path.basename(args.input[0]))[0]

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

# Create a Pypiper object, and start the pipeline (runs initial setting and logging code to begin)
mypiper = pypiper.PipelineManager(name="rnaTopHat", outfolder=paths.pipeline_outfolder, args=args)
ngstk = pypiper.NGSTk(pm=mypiper)

print("N input bams:\t\t" + str(len(args.input)))
print("Sample name:\t\t" + args.sample_name)

sample_merged_bam = args.sample_name + ".merged.bam"
mypiper.make_sure_path_exists(paths.pipeline_outfolder + "unmapped_bam/")

if merge and not os.path.isfile(sample_merged_bam):
	print("Multiple unmapped bams found; merge requested")
	input_bams = args.input;
	print("input bams: " + str(input_bams))
	merge_folder = os.path.join(paths.pipeline_outfolder, "unmapped_bam/")
	input_string = " INPUT=" + " INPUT=".join(input_bams)
	output_merge = os.path.join(merge_folder, sample_merged_bam)
	cmd = "java -Xmx2g -jar " + paths.picard_jar + " MergeSamFiles"
	cmd += input_string
	cmd += " OUTPUT=" + output_merge
	cmd += " ASSUME_SORTED=TRUE"
	cmd += " CREATE_INDEX=TRUE"

	mypiper.run(cmd, output_merge)
	args.input = paths.pipeline_outfolder+"unmapped_bam/" + sample_merged_bam  #update unmapped bam reference
	local_unmapped_bam = paths.pipeline_outfolder+"unmapped_bam/" + sample_merged_bam
else:
	# Link the file into the unmapped_bam directory
	print("Single unmapped bam found; no merge required")
	print("Unmapped bam:\t\t" + str(args.input[0]))
	args.input = args.input[0]
	local_unmapped_bam = paths.pipeline_outfolder+"unmapped_bam/"+args.sample_name+".bam"
	mypiper.callprint("ln -sf " + args.input + " " + local_unmapped_bam, shell=True)


print("Input Unmapped Bam: " + args.input)
#check for file exists:
if not os.path.isfile(local_unmapped_bam):
	print local_unmapped_bam + "is not a file"

# Fastq conversion
########################################################################################
# Uses ngstk module.
mypiper.timestamp("### Fastq conversion: ")

# Sanity checks:
def check_fastq():
	raw_reads = ngstk.count_reads(local_unmapped_bam,args.paired_end)
	mypiper.report_result("Raw_reads", str(raw_reads))
	fastq_reads = ngstk.count_reads(out_fastq_pre + "_R1.fastq", paired_end=args.paired_end)
	mypiper.report_result("Fastq_reads", fastq_reads)
	if (fastq_reads != int(raw_reads)):
		raise Exception("Fastq conversion error? Size doesn't match unaligned bam")

out_fastq_pre = os.path.join(paths.pipeline_outfolder, "fastq/", args.sample_name)
cmd = ngstk.bam_to_fastq(local_unmapped_bam, out_fastq_pre, args.paired_end)
mypiper.run(cmd, out_fastq_pre + "_R1.fastq")

check_fastq()

mypiper.clean_add(out_fastq_pre + "*.fastq", conditional=True)


# Adapter trimming
########################################################################################
mypiper.timestamp("### Adapter trimming: ")

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
mypiper.run(cmd, out_fastq_pre + "_R1_trimmed.fastq")

mypiper.report_result("Trimmed_reads", ngstk.count_reads(trimmed_fastq,args.paired_end))


# RNA Tophat pipeline.
########################################################################################
mypiper.timestamp("### TopHat alignment: ")
tophat_folder = paths.pipeline_outfolder + "/tophat_" + args.genome_assembly + "/"
mypiper.make_sure_path_exists(tophat_folder)
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

mypiper.run(cmd, tophat_folder + "/align_summary.txt", shell=False)

mypiper.timestamp("### renaming tophat aligned bam file ")
cmd = "mv " + tophat_folder + "/accepted_hits.bam " + out_tophat
mypiper.run(cmd, out_tophat, shell=False, follow=lambda:
	mypiper.report_result("Aligned_reads", ngstk.count_unique_mapped_reads(out_tophat,args.paired_end and not align_paired_as_single)))





mypiper.timestamp("### BAM to SAM sorting and indexing: ")
if args.filter:
	cmd = ngstk.bam_conversions(out_tophat,False)
	mypiper.run(cmd,  re.sub(".bam$", "_sorted.bam", out_tophat) ,shell=True)
else:
	cmd = ngstk.bam_conversions(out_tophat,True)
	mypiper.run(cmd, re.sub(".bam$", "_sorted.depth", out_tophat),shell=True)

mypiper.clean_add(out_tophat, conditional=False)
mypiper.clean_add(re.sub(".bam$" , ".sam", out_tophat), conditional=False)



if not args.filter and args.markDupl:
	mypiper.timestamp("### MarkDuplicates: ")

	aligned_file = re.sub(".sam$", "_sorted.bam",  out_tophat)
	out_file = re.sub(".sam$", "_dedup.bam", out_tophat)
	metrics_file = re.sub(".sam$", "_dedup.metrics", out_tophat)
	cmd = ngstk.markDuplicates(aligned_file, out_file, metrics_file)
	mypiper.run(cmd, out_file, follow= lambda:
		mypiper.report_result("Deduplicated_reads", ngstk.count_unique_mapped_reads(out_file, args.paired_end and not align_paired_as_single)))




#read filtering
########################################################################################

if args.filter:
	mypiper.timestamp("### Aligned read filtering: ")

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

	mypiper.run(cmd, out_sam_filter, follow=lambda:
		mypiper.report_result("Filtered_reads", ngstk.count_unique_mapped_reads(out_sam_filter, args.paired_end and not align_paired_as_single)))



#	mypiper.timestamp("### Set flag in skipped reads to 4 (unmapped): ")
#	joined_sam = out_sam_filter.replace(".filt." , ".filtFlag.")
#	skipped_sam = out_sam_filter.replace(".filt." , ".skipped.")
#	cmd = "samtools view -hS " + out_sam_filter + " > " + joined_sam + "\n"
#	cmd += "awk -v skip=" + headerLines + " -v OFS=\"\\t\" '{if (NR>skip){$2=4;$3=\"*\";$4=0;$5=0;$6=\"*\";$7=\"*\";$8=0;$9=0; print}}' " + skipped_sam
#	cmd += " >>" + joined_sam
#
#	mypiper.call_lock(cmd, joined_sam , shell=True)
#
#	if not args.no_check:
#		x = ngstk.count_reads(joined_sam, args.paired_end and not align_paired_as_single)
#		mypiper.report_result("Joined_reads", x)
#
#	mypiper.timestamp("### Filtered: SAM to BAM conversion, sorting and depth calculation: ")
#	cmd = ngstk.sam_conversions(joined_sam)
#	mypiper.call_lock(cmd, joined_sam.replace(".sam" , "_sorted.depth"),shell=True)
#
#	mypiper.timestamp("### Skipped: SAM to BAM conversion and sorting: ")
#	cmd = ngstk.sam_conversions(skipped_sam,False)
#	mypiper.call_lock(cmd, skipped_sam.replace(".sam" , "_sorted.bam"),shell=True)
#
#	if args.markDupl:
#		mypiper.timestamp("### MarkDuplicates: ")
#
#		aligned_file = joined_sam.replace(".sam" , "_sorted.bam")
#		out_file = joined_sam.replace(".sam" , "_dedup.bam")
#		metrics_file = joined_sam.replace(".sam" , "_dedup.metrics")
#		cmd = ngstk.markDuplicates(paths, aligned_file, out_file, metrics_file)
#		mypiper.call_lock(cmd, out_file)
#
#		if not args.no_check:
#			x = ngstk.count_mapped_reads(out_file)
#			mypiper.report_result("Deduplicated_reads", x)

	mypiper.timestamp("### Filtered: SAM to BAM conversion, sorting and depth calculation: ")
	cmd = ngstk.sam_conversions(out_sam_filter)
	mypiper.run(cmd, re.sub(".sam$", "_sorted.depth", out_sam_filter),shell=True)

	skipped_sam = out_sam_filter.replace(".filt." , ".skipped.")
	mypiper.timestamp("### Skipped: SAM to BAM conversion and sorting: ")
	cmd = ngstk.sam_conversions(skipped_sam,False)
	mypiper.run(cmd, re.sub(".sam$" , "_sorted.bam", skipped_sam),shell=True)

	mypiper.clean_add(skipped_sam, conditional=False)
	mypiper.clean_add(re.sub(".sam$" , ".bam", skipped_sam), conditional=False)
	mypiper.clean_add(out_sam_filter, conditional=False)
	mypiper.clean_add(re.sub(".sam$" , ".bam", out_sam_filter), conditional=False)



#create tracks
########################################################################################
mypiper.timestamp("### bam2wig: ")
if args.filter:
	trackFile = re.sub(".sam$", "_sorted.bam", out_sam_filter)
	cmd = paths.bam2wig + " -i" + trackFile
	cmd += " -s " + paths.chrom_sizes
	cmd += " -o " + re.sub(".sam$" , "_sorted", out_sam_filter)
	cmd += " -t " + str(args.wigsum)
	mypiper.run(cmd, re.sub(".sam$" , "_sorted.wig",out_sam_filter),shell=False)

	mypiper.timestamp("### wigToBigWig: ")
	cmd = paths.wigToBigWig + " " + re.sub(".sam$" , "_sorted.wig",out_sam_filter)
	cmd += " " + paths.chrom_sizes
	cmd += " " + re.sub(".sam$" , "_sorted.bw",out_sam_filter)
	mypiper.run(cmd, re.sub(".sam$" , "_sorted.bw",out_sam_filter),shell=False)

else:
	trackFile = re.sub(".bam$", "_sorted.bam",out_tophat)
	cmd = paths.bam2wig + " -i" + trackFile
	cmd += " -s " + paths.chrom_sizes
	cmd += " -o " + re.sub(".bam$" , "_sorted",out_tophat)
	cmd += " -t " + str(args.wigsum)
	mypiper.run(cmd, re.sub(".bam$" , "_sorted.wig",out_tophat),shell=False)

	mypiper.timestamp("### wigToBigWig: ")
	cmd = paths.wigToBigWig + " " + re.sub(".bam$" , "_sorted.wig",out_tophat)
	cmd += " " + paths.chrom_sizes
	cmd += " " + re.sub(".bam$" , "_sorted.bw",out_tophat)
	mypiper.run(cmd, re.sub(".bam$" , "_sorted.bw", out_tophat),shell=False)


mypiper.timestamp("### read_distribution: ")
cmd = paths.read_distribution + " -i " + trackFile
cmd	+= " -r " + paths.gene_model_bed
cmd += " > " + re.sub("_sorted.bam$", "_read_distribution.txt",trackFile)
mypiper.run(cmd, re.sub("_sorted.bam$", "_read_distribution.txt",trackFile),shell=True)


mypiper.timestamp("### gene_coverage: ")
cmd = paths.gene_coverage + " -i " + re.sub(".bam$" , ".bw",trackFile)
cmd	+= " -r " + paths.gene_model_sub_bed
cmd += " -o " + re.sub("_sorted.bam$", "",trackFile)
mypiper.run(cmd, re.sub("_sorted.bam$", ".geneBodyCoverage.png",trackFile),shell=False)


# Cleanup
########################################################################################


# remove temporary marker file:
mypiper.stop_pipeline()
