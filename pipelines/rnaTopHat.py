#!/usr/bin/env python

"""
RNA TopHat pipeline
"""

from argparse import ArgumentParser
import os
import sys
import subprocess
import re
import pypiper


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
	parser.print_help()
	raise SystemExit

# Initialize
pm = pypiper.PipelineManager(name="rnaTopHat", outfolder=os.path.join(args.output_parent, args.sample_name), args=args)

# Tools
pm.config.tools.scripts_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "tools")

# Resources
pm.config.resources.ref_genome = os.path.join(pm.config.resources.genomes, args.genome_assembly)
pm.config.resources.ref_genome_fasta = os.path.join(pm.config.resources.genomes, args.genome_assembly, args.genome_assembly + ".fa")
pm.config.resources.chrom_sizes = os.path.join(pm.config.resources.genomes, args.genome_assembly, args.genome_assembly + ".chromSizes")
pm.config.resources.bowtie_indexed_genome = os.path.join(pm.config.resources.genomes, args.genome_assembly, "indexed_bowtie2", args.genome_assembly)
# tophat specific resources
pm.config.resources.gtf = os.path.join(pm.config.resources.genomes, args.genome_assembly , "ucsc_" + args.genome_assembly + "_ensembl_genes.gtf")
pm.config.resources.gene_model_bed = os.path.join(pm.config.resources.genomes, args.genome_assembly , "ucsc_" + args.genome_assembly + "_ensembl_genes.bed")
pm.config.resources.gene_model_sub_bed = os.path.join(pm.config.resources.genomes, args.genome_assembly , "ucsc_" + args.genome_assembly + "_ensembl_genes_500rand.bed")

# Output
pm.config.parameters.pipeline_outfolder = os.path.join(args.output_parent, args.sample_name)

# Create a ngstk object
ngstk = pypiper.NGSTk(pm=pm)

print "\nPipeline configuration:"
print(pm.config)
tools = pm.config.tools  # Convenience alias
param = pm.config.parameters
resources = pm.config.resources

print("Sample name:\t" + args.sample_name)

# Merge/Link sample input
################################################################################
# This command should now handle all the merging.
local_input_file = ngstk.create_local_input(param.pipeline_outfolder, args.input, args.sample_name)

print("Local input file: " + local_input_file) 

# Make sure file exists:
if not os.path.isfile(local_input_file):
	print local_input_file + " is not a file"

# Fastq conversion
################################################################################
pm.timestamp("### Fastq conversion: ")
# New fastq conversion (can handle .bam or .fastq.gz files)

cmd, fastq_folder, out_fastq_pre, unaligned_fastq = ngstk.input_to_fastq(local_input_file, param.pipeline_outfolder, args.sample_name, args.paired_end)
ngstk.make_sure_path_exists(fastq_folder)
pm.run(cmd, unaligned_fastq, follow=ngstk.check_fastq(local_input_file, unaligned_fastq, args.paired_end))
pm.clean_add(out_fastq_pre + "*.fastq", conditional=True)

# Adapter trimming
################################################################################
pm.timestamp("### Adapter trimming: ")

cmd = tools.java + " -Xmx" + str(pm.mem) + " -jar " + tools.trimmomatic_epignome

if not args.paired_end:
	cmd += " SE -phred33 -threads " + str(pm.cores) + " "
	cmd += out_fastq_pre + "_R1.fastq "
	cmd += out_fastq_pre + "_R1_trimmed.fastq "

else:
	cmd += " PE -phred33 -threads " + str(pm.cores) + " "
	cmd += out_fastq_pre + "_R1.fastq "
	cmd += out_fastq_pre + "_R2.fastq "
	cmd += out_fastq_pre + "_R1_trimmed.fastq "
	cmd += out_fastq_pre + "_R1_unpaired.fastq "
	cmd += out_fastq_pre + "_R2_trimmed.fastq "
	cmd += out_fastq_pre + "_R2_unpaired.fastq "

# for Core-seq, trim off the first 6bp and the bit adjacent to identified adapter sequences:
if args.coreseq:
	cmd += " HEADCROP:6"
	cmd += " ILLUMINACLIP:" + resources.adapters + ":2:10:4:1:true:epignome:5"
	cmd += " SLIDINGWINDOW:4:1"
	cmd += " MAXINFO:16:0.40"
	cmd += " MINLEN:25"
# otherwise just look for normal adapters:
else:
	cmd += " HEADCROP:6"
	cmd += " ILLUMINACLIP:" + resources.adapters + ":2:10:4:1:true"
	cmd += " ILLUMINACLIP:" + "/data/groups/lab_bsf/resources/trimmomatic_adapters/PolyA-SE.fa" + ":2:30:5:1:true"
	cmd += " SLIDINGWINDOW:4:1"
	cmd += " MAXINFO:16:0.40"
	cmd += " MINLEN:21"

trimmed_fastq = out_fastq_pre + "_R1_trimmed.fastq"
pm.run(cmd, out_fastq_pre + "_R1_trimmed.fastq")

pm.report_result("Trimmed_reads", ngstk.count_reads(trimmed_fastq,args.paired_end))


# RNA Tophat pipeline.
########################################################################################
pm.timestamp("### TopHat alignment: ")
tophat_folder = os.path.join(param.pipeline_outfolder,"tophat_" + args.genome_assembly)
pm.make_sure_path_exists(tophat_folder)
out_tophat = os.path.join(tophat_folder,args.sample_name + ".aln.bam")

align_paired_as_single = True # FH: this appears to be the default behavior of the pipeline at the moment. Should that be configurable by args?

if not args.paired_end:
	cmd = tools.tophat2
	cmd += " --GTF " + resources.gtf
	cmd += " --b2-L 15 --library-type fr-unstranded --mate-inner-dist 150 --max-multihits 100 --no-coverage-search --num-threads " + str(pm.cores)
	cmd += " --output-dir " + tophat_folder
	cmd += " " + resources.bowtie_indexed_genome
	cmd += " " + out_fastq_pre + "_R1_trimmed.fastq"

else:
	cmd = tools.tophat2
	cmd += " --GTF " + resources.gtf
	cmd += " --b2-L 15 --library-type fr-unstranded --mate-inner-dist 150 --max-multihits 100 --no-coverage-search --num-threads " + str(pm.cores)
	cmd += " --output-dir " + tophat_folder
	cmd += " " + resources.bowtie_indexed_genome
	# FH: if you use this code, you align both mates separately. As a result, the count_unique_mapped_reads method in paired-end mode will return 0, because the mate flags are not set
	if align_paired_as_single:
		cmd += " " + out_fastq_pre + "_R1_trimmed.fastq" + "," + out_fastq_pre + "_R2_trimmed.fastq"
	else:
		cmd += " " + out_fastq_pre + "_R1_trimmed.fastq"
		cmd += " " + out_fastq_pre + "_R2_trimmed.fastq"

pm.run(cmd, os.path.join(tophat_folder,"align_summary.txt"), shell=False)

pm.timestamp("### renaming tophat aligned bam file ")
cmd = "mv " + os.path.join(tophat_folder,"accepted_hits.bam") + " " + out_tophat
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
	headerLines = subprocess.check_output(tools.samtools + " view -SH " + re.sub(".bam$", ".sam", out_tophat) + "| wc -l", shell=True).strip()
	cmd = tools.python + " " + os.path.join(tools.scripts_dir,"bisulfiteReadFiltering_forRNA.py")
	cmd += " --infile=" + re.sub(".bam$",".sam", out_tophat)
	cmd += " --outfile=" + out_sam_filter
	cmd += " --skipped=" + out_sam_filter.replace(".filt." , ".skipped.")
	cmd += " --skipHeaderLines=" + headerLines
	cmd += " --genome=" + args.genome_assembly
	cmd += " --genomeDir=" + resources.ref_genome
	cmd += " --minNonCpgSites=3"
	cmd += " --minConversionRate=0.9"
	cmd += " --maxConversionRate=0.1"
	cmd += " -r"

	if args.paired_end and not align_paired_as_single:
		cmd = cmd + " --pairedEnd"

	pm.run(cmd, out_sam_filter, follow=lambda:
		pm.report_result("Filtered_reads", ngstk.count_unique_mapped_reads(out_sam_filter, args.paired_end and not align_paired_as_single)))

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
	cmd = tools.bam2wig + " -i " + trackFile
	cmd += " -s " + resources.chrom_sizes
	cmd += " -o " + re.sub(".sam$" , "_sorted", out_sam_filter)
	cmd += " -t " + str(args.wigsum)
	pm.run(cmd, re.sub(".sam$" , "_sorted.wig",out_sam_filter),shell=False)

	pm.timestamp("### wigToBigWig: ")
	cmd = tools.wigToBigWig + " " + re.sub(".sam$" , "_sorted.wig",out_sam_filter)
	cmd += " " + resources.chrom_sizes
	cmd += " " + re.sub(".sam$" , "_sorted.bw",out_sam_filter)
	pm.run(cmd, re.sub(".sam$" , "_sorted.bw",out_sam_filter),shell=False)

else:
	trackFile = re.sub(".bam$", "_sorted.bam",out_tophat)
	cmd = tools.bam2wig + " -i " + trackFile
	cmd += " -s " + resources.chrom_sizes
	cmd += " -o " + re.sub(".bam$" , "_sorted",out_tophat)
	cmd += " -t " + str(args.wigsum)
	pm.run(cmd, re.sub(".bam$" , "_sorted.wig",out_tophat),shell=False)

	pm.timestamp("### wigToBigWig: ")
	cmd = tools.wigToBigWig + " " + re.sub(".bam$" , "_sorted.wig",out_tophat)
	cmd += " " + resources.chrom_sizes
	cmd += " " + re.sub(".bam$" , "_sorted.bw",out_tophat)
	pm.run(cmd, re.sub(".bam$" , "_sorted.bw", out_tophat),shell=False)

pm.timestamp("### read_distribution: ")
cmd = tools.read_distribution + " -i " + trackFile
cmd += " -r " + resources.gene_model_bed
cmd += " > " + re.sub("_sorted.bam$", "_read_distribution.txt",trackFile)
pm.run(cmd, re.sub("_sorted.bam$", "_read_distribution.txt",trackFile),shell=True)

pm.timestamp("### gene_coverage: ")
cmd = tools.gene_coverage + " -i " + re.sub(".bam$" , ".bw",trackFile)
cmd += " -r " + resources.gene_model_sub_bed
cmd += " -o " + re.sub("_sorted.bam$", "",trackFile)
pm.run(cmd, re.sub("_sorted.bam$", ".geneBodyCoverage.png",trackFile),shell=False)


# Cleanup
########################################################################################


# remove temporary marker file:
pm.stop_pipeline()
