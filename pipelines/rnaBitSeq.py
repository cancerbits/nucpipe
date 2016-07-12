#!/usr/bin/env python

"""
RNA BitSeq pipeline
"""

from argparse import ArgumentParser
import os
import sys
import yaml
import subprocess
import re
import pypiper


# Argument Parsing
# #######################################################################################
parser = ArgumentParser(description='Pypiper arguments.')

parser = pypiper.add_pypiper_args(parser, all_args=True)

# Add any pipeline-specific arguments

parser.add_argument('-e', '--ercc',
				default = "ERCC92",
				dest = 'ERCC_assembly',
				type = str,
				help = 'ERCC Assembly')
parser.add_argument('-em', '--ercc-mix',
				default = "False",
				dest = 'ERCC_mix',
				help = 'ERCC mix. If False no ERCC analysis will be performed.')
parser.add_argument('-f', dest='filter', action='store_false', default=True)

# Core-seq as optional parameter
parser.add_argument('-cs', '--core-seq', default=False, dest='coreseq', action='store_true', help='CORE-seq Mode')

# Quant-Seq as optional parameter
parser.add_argument('-qs', '--quantseq', default=False, dest='quantseq', action='store_true', help='Quant-Seq Mode')

args = parser.parse_args()

if args.single_or_paired == "paired":
	args.paired_end = True
else:
	args.paired_end = False

if not args.input:
	parser.print_help()
	raise SystemExit

# Initialize
outfolder = os.path.abspath(os.path.join(args.output_parent, args.sample_name))
pm = pypiper.PipelineManager(name = "rnaBitSeq", outfolder = outfolder, args = args)

# Tools
pm.config.tools.scripts_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "tools")

# Resources
pm.config.resources.ref_genome = os.path.join(pm.config.resources.genomes, args.genome_assembly)
pm.config.resources.ref_genome_fasta = os.path.join(pm.config.resources.genomes, args.genome_assembly, args.genome_assembly + ".fa")
pm.config.resources.ref_ERCC_fasta = os.path.join(pm.config.resources.genomes, args.ERCC_assembly, args.ERCC_assembly + ".fa")
pm.config.resources.chrom_sizes = os.path.join(pm.config.resources.genomes, args.genome_assembly, args.genome_assembly + ".chromSizes")
pm.config.resources.bowtie_indexed_genome = os.path.join(pm.config.resources.genomes, args.genome_assembly, "indexed_bowtie1", args.genome_assembly)
pm.config.resources.bowtie_indexed_ERCC = os.path.join(pm.config.resources.genomes, args.ERCC_assembly, "indexed_bowtie1", args.ERCC_assembly)

# Output
pm.config.parameters.pipeline_outfolder = outfolder

ngstk = pypiper.NGSTk(pm=pm)

tools = pm.config.tools  # Convenience alias
param = pm.config.parameters
resources = pm.config.resources

raw_folder = os.path.join(param.pipeline_outfolder, "raw/")
fastq_folder = os.path.join(param.pipeline_outfolder, "fastq/")

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
	if args.quantseq: cmd += " HEADCROP:6"
	cmd += " ILLUMINACLIP:" + resources.adapters + ":2:10:4:1:true"
	if args.quantseq: cmd += " ILLUMINACLIP:" + "/data/groups/lab_bsf/resources/trimmomatic_adapters/PolyA-SE.fa" + ":2:30:5:1:true"
	cmd += " SLIDINGWINDOW:4:1"
	cmd += " MAXINFO:16:0.40"
	cmd += " MINLEN:21"

trimmed_fastq = out_fastq_pre + "_R1_trimmed.fastq"
trimmed_fastq_R2 = out_fastq_pre + "_R2_trimmed.fastq"

#pm.run(cmd, out_fastq_pre + "_R1_trimmed.fastq",
#	follow = lambda: pm.report_result("Trimmed_reads", ngstk.count_reads(trimmed_fastq,args.paired_end)))

pm.run(cmd, trimmed_fastq, 
	follow = ngstk.check_trim(trimmed_fastq, trimmed_fastq_R2, args.paired_end,
		fastqc_folder = os.path.join(param.pipeline_outfolder, "fastqc/")))


# RNA BitSeq pipeline.
########################################################################################
pm.timestamp("### Bowtie1 alignment: ")

bowtie1_folder = os.path.join(param.pipeline_outfolder,"bowtie1_" + args.genome_assembly)
pm.make_sure_path_exists(bowtie1_folder)
out_bowtie1 = os.path.join(bowtie1_folder, args.sample_name + ".aln.sam")

if not args.paired_end:
	cmd = tools.bowtie1
	cmd += " -q -p " + str(pm.cores) + " -a -m 100 --sam "
	cmd += resources.bowtie_indexed_genome + " "
	cmd += out_fastq_pre + "_R1_trimmed.fastq"
	cmd += " " + out_bowtie1
else:
	cmd = tools.bowtie1
	cmd += " -q -p " + str(pm.cores) + " -a -m 100 --minins 0 --maxins 5000 --fr --sam --chunkmbs 200 "    # also checked --rf (1% aln) and --ff (0% aln) --fr(8% aln)
	cmd += resources.bowtie_indexed_genome
	cmd += " -1 " + out_fastq_pre + "_R1_trimmed.fastq"
	cmd += " -2 " + out_fastq_pre + "_R2_trimmed.fastq"
	cmd += " " + out_bowtie1

pm.run(cmd, out_bowtie1,
	follow=lambda: pm.report_result("Aligned_reads", ngstk.count_unique_mapped_reads(out_bowtie1, args.paired_end)))

pm.timestamp("### Raw: SAM to BAM conversion and sorting: ")

if args.filter:
	cmd = ngstk.sam_conversions(out_bowtie1,False)
	pm.run(cmd,  re.sub(".sam$" , "_sorted.bam",out_bowtie1),shell=True)
else:
	cmd = ngstk.sam_conversions(out_bowtie1,True)
	pm.run(cmd,  re.sub(".sam$" , "_sorted.depth",out_bowtie1),shell=True)

pm.clean_add(out_bowtie1, conditional=False)
pm.clean_add(re.sub(".sam$" , ".bam", out_bowtie1), conditional=False)


if not args.filter:
	pm.timestamp("### MarkDuplicates: ")
	aligned_file = re.sub(".sam$" , "_sorted.bam",out_bowtie1)
	out_file = re.sub(".sam$" , "_dedup.bam",out_bowtie1)
	metrics_file = re.sub(".sam$" , "_dedup.metrics",out_bowtie1)
	cmd = ngstk.markDuplicates(aligned_file, out_file, metrics_file)
	pm.run(cmd, out_file, follow= lambda: pm.report_result("Deduplicated_reads", ngstk.count_unique_mapped_reads(out_file, args.paired_end)))

if args.filter:
	pm.timestamp("### Aligned read filtering: ")
	out_sam_filter = bowtie1_folder + args.sample_name + ".aln.filt.sam"
	skipped_sam = out_sam_filter.replace(".filt." , ".skipped.")
	headerLines = subprocess.check_output("samtools view -SH " + out_bowtie1 + "|wc -l", shell=True).strip()
	cmd = tools.python + " " + os.path.join(tools.scripts_dir,"bisulfiteReadFiltering_forRNA.py")
	cmd += " --infile=" + out_bowtie1
	cmd += " --outfile=" + out_sam_filter
	cmd += " --skipped=" + skipped_sam
	cmd += " --skipHeaderLines=" + headerLines
	cmd += " --genome=" + args.genome_assembly
	cmd += " --genomeDir=" + resources.ref_genome
	cmd += " --minNonCpgSites=3"
	cmd += " --minConversionRate=0.9"
	cmd += " --maxConversionRate=0.1"
	cmd += " -r"

	if args.paired_end:
		cmd = cmd + " --pairedEnd"

	pm.run(cmd, out_sam_filter,follow=pm.report_result("Filtered_reads", ngstk.count_unique_mapped_reads(out_sam_filter, args.paired_end)))

	pm.timestamp("### Filtered: SAM to BAM conversion, sorting and depth calculation: ")
	cmd = ngstk.sam_conversions(out_sam_filter)
	pm.run(cmd, re.sub(".sam$" , "_sorted.depth",out_sam_filter),shell=True)


	pm.timestamp("### Skipped: SAM to BAM conversion and sorting: ")
	cmd = ngstk.sam_conversions(skipped_sam,False)
	pm.run(cmd, re.sub(".sam$", "_sorted.bam", skipped_sam),shell=True)
	
	pm.clean_add(skipped_sam, conditional=False)
	pm.clean_add(re.sub(".sam$" , ".bam", skipped_sam), conditional=False)
	pm.clean_add(out_sam_filter, conditional=False)
	pm.clean_add(re.sub(".sam$" , ".bam", out_sam_filter), conditional=False)

	pm.timestamp("### MarkDuplicates: ")
	
	aligned_file = re.sub(".sam$" , "_sorted.bam",out_sam_filter)
	out_file = re.sub(".sam$" , "_dedup.bam",out_sam_filter)
	metrics_file = re.sub(".sam$" , "_dedup.metrics",out_sam_filter)
	cmd = ngstk.markDuplicates(aligned_file, out_file, metrics_file)
	pm.run(cmd, out_file,follow=lambda: pm.report_result("Deduplicated_reads", ngstk.count_unique_mapped_reads(out_file, args.paired_end)))

# BitSeq
########################################################################################
pm.timestamp("### Expression analysis (BitSeq): ")

bitSeq_dir = os.path.join(bowtie1_folder,"bitSeq")
pm.make_sure_path_exists(bitSeq_dir)
out_bitSeq = os.path.join(bitSeq_dir, args.sample_name + ".counts")

if args.filter:
	cmd = tools.Rscript + " " + os.path.join(tools.scripts_dir,"bitSeq_parallel.R") + " " + out_sam_filter + " " + bitSeq_dir + " " + resources.ref_genome_fasta
else:
	cmd = tools.Rscript + " " + os.path.join(tools.scripts_dir,"bitSeq_parallel.R") + " " + out_bowtie1 + " " + bitSeq_dir + " " + resources.ref_genome_fasta

pm.run(cmd, out_bitSeq)


# ERCC Spike-in alignment
########################################################################################
if not (args.ERCC_mix == "False" ):
	pm.timestamp("### ERCC: Convert unmapped reads into fastq files: ")

	# Sanity checks:
	def check_fastq_ERCC():
		raw_reads = ngstk.count_reads(unmappable_bam + ".bam",args.paired_end)
		pm.report_result("ERCC_raw_reads", str(raw_reads))
		fastq_reads = ngstk.count_reads(unmappable_bam + "_R1.fastq", paired_end=args.paired_end)
		pm.report_result("ERCC_fastq_reads", fastq_reads)
		if (fastq_reads != int(raw_reads)):
			raise Exception("Fastq conversion error? Size doesn't match unaligned bam")

	unmappable_bam = re.sub(".sam$","_unmappable",out_bowtie1)
	cmd = tools.samtools + " view -hbS -f4 " + out_bowtie1 + " > " + unmappable_bam + ".bam"
	pm.run(cmd, unmappable_bam + ".bam", shell=True)

	cmd = ngstk.bam_to_fastq(unmappable_bam + ".bam", unmappable_bam, args.paired_end)
	pm.run(cmd, unmappable_bam + "_R1.fastq",follow=check_fastq_ERCC)

	pm.timestamp("### ERCC: Bowtie1 alignment: ")
	bowtie1_folder = os.path.join(param.pipeline_outfolder,"bowtie1_" + args.ERCC_assembly)
	pm.make_sure_path_exists(bowtie1_folder)
	out_bowtie1 = os.path.join(bowtie1_folder, args.sample_name + "_ERCC.aln.sam")

	if not args.paired_end:
		cmd = tools.bowtie1
		cmd += " -q -p " + str(pm.cores) + " -a -m 100 --sam "
		cmd += resources.bowtie_indexed_ERCC + " "
		cmd += unmappable_bam + "_R1.fastq"
		cmd += " -S " + out_bowtie2
	else:
		cmd = tools.bowtie1
		cmd += " -q -p " + str(pm.cores) + " -a -m 100 --minins 0 --maxins 5000 --fr --sam --chunkmbs 200 "
		cmd += resources.bowtie_indexed_ERCC
		cmd += " -1 " + unmappable_bam + "_R1.fastq"
		cmd += " -2 " + unmappable_bam + "_R2.fastq"
		cmd += " -S " + out_bowtie2


#	if not args.paired_end:
#		cmd = param.bowtie1
#		cmd += " -q -p 6 -a -m 100 --sam "
#		cmd += param.bowtie_indexed_ERCC + " "
#		cmd += unmappable_bam + "_R1.fastq"
#		cmd += " " + out_bowtie1
#	else:
#		cmd = param.bowtie1
#		cmd += " -q -p 6 -a -m 100 --minins 0 --maxins 5000 --fr --sam --chunkmbs 200 "
#		cmd += param.bowtie_indexed_ERCC
#		cmd += " -1 " + unmappable_bam + "_R1.fastq"
#		cmd += " -2 " + unmappable_bam + "_R2.fastq"
#		cmd += " " + out_bowtie1

	pm.run(cmd, out_bowtie1,follow=lambda: pm.report_result("ERCC_aligned_reads", ngstk.count_unique_mapped_reads(out_bowtie1, args.paired_end)))

	pm.timestamp("### ERCC: SAM to BAM conversion, sorting and depth calculation: ")
	cmd = ngstk.sam_conversions(out_bowtie1)
	pm.run(cmd, re.sub(".sam$" , "_sorted.depth", out_bowtie1), shell=True)

	pm.clean_add(out_bowtie1, conditional=False)
	pm.clean_add(re.sub(".sam$" , ".bam", out_bowtie1), conditional=False)
	pm.clean_add(unmappable_bam + "*.fastq", conditional=False)

# BitSeq
########################################################################################
	pm.timestamp("### ERCC: Expression analysis (BitSeq): ")

	bitSeq_dir = os.path.join(bowtie1_folder,"bitSeq")
	pm.make_sure_path_exists(bitSeq_dir)
	out_bitSeq = os.path.join(bitSeq_dir,re.sub(".aln.sam$" , ".counts",out_bowtie1))

	cmd = tools.Rscript + " " + os.path.join(tools.scripts_dir,"bitSeq_parallel.R") + " " + out_bowtie1 + " " + bitSeq_dir + " " + resources.ref_ERCC_fasta
	pm.run(cmd, out_bitSeq)


# Cleanup
########################################################################################
# remove temporary marker file:
pm.stop_pipeline()





