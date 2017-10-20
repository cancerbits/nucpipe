#!/usr/bin/env python
"""
RNA pipeline combining Tophat and the End Sequence Analysis Toolkit (ESAT)
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

parser.add_argument('-d', dest='markDupl', action='store_true', default=False)
parser.add_argument('-w', '--wigsum', default=500000000, dest='wigsum', type=int, help='Target wigsum for track normalisation')

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
pm = pypiper.PipelineManager(name = "rnaESAT", outfolder = outfolder, args = args)

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
pm.config.parameters.pipeline_outfolder = outfolder

ngstk = pypiper.NGSTk(pm=pm)

tools = pm.config.tools  # Convenience aliases
param = pm.config.parameters
resources = pm.config.resources

raw_folder = os.path.join(param.pipeline_outfolder, "raw")
fastq_folder = os.path.join(param.pipeline_outfolder, "fastq")

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
pm.timestamp("### Trimming: ")

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

cmd += " HEADCROP:6"
cmd += " ILLUMINACLIP:" + resources.adapters + ":2:10:4:1:true"
cmd += " ILLUMINACLIP:" + resources.polyA + ":2:30:5:1:true"
cmd += " SLIDINGWINDOW:4:1"
cmd += " MAXINFO:16:0.40"
cmd += " MINLEN:21"

trimmed_fastq = out_fastq_pre + "_R1_trimmed.fastq"
trimmed_fastq_R2 = out_fastq_pre + "_R2_trimmed.fastq"

pm.run(cmd, trimmed_fastq, 
	follow = ngstk.check_trim(trimmed_fastq, args.paired_end, trimmed_fastq_R2,
		fastqc_folder = os.path.join(param.pipeline_outfolder, "fastqc/")))


# Tophat alignment
########################################################################################
pm.timestamp("### TopHat alignment: ")

tophat_folder = os.path.join(param.pipeline_outfolder,"tophat_" + args.genome_assembly)
pm.make_sure_path_exists(tophat_folder)
out_tophat = os.path.join(tophat_folder,args.sample_name + ".aln.bam")

align_paired_as_single = True # FH: this appears to be the default behavior of the pipeline at the moment. Should that be configurable by args?

cmd = tools.tophat2
cmd += " --GTF " + resources.gtf
cmd += " --b2-L " + str(param.tophat.b2L)
cmd += " --library-type " + str(param.tophat.librarytype)
cmd += " --mate-inner-dist " + str(param.tophat.mateinnerdist)
cmd += " --max-multihits " + str(param.tophat.maxmultihits)
cmd += " --no-coverage-search"
cmd += " --num-threads " + str(pm.cores)
cmd += " --output-dir " + tophat_folder
cmd += " " + resources.bowtie_indexed_genome
if not args.paired_end:
	cmd += " " + out_fastq_pre + "_R1_trimmed.fastq"
else:
	# FH: if you use this code, you align both mates separately. As a result, the count_unique_mapped_reads method in paired-end mode will return 0, because the mate flags are not set
	if align_paired_as_single:
		cmd += " " + out_fastq_pre + "_R1_trimmed.fastq" + "," + out_fastq_pre + "_R2_trimmed.fastq"
	else:
		cmd += " " + out_fastq_pre + "_R1_trimmed.fastq"
		cmd += " " + out_fastq_pre + "_R2_trimmed.fastq"

pm.run(cmd, os.path.join(tophat_folder,"align_summary.txt"), shell=False)

pm.timestamp("### renaming tophat aligned bam file ")

cmd = "mv " + os.path.join(tophat_folder,"accepted_hits.bam") + " " + out_tophat

def check_tophat():
	ar = ngstk.count_unique_mapped_reads(out_tophat,args.paired_end and not align_paired_as_single)
	pm.report_result("Aligned_reads", ar)
	rr = float(pm.get_stat("Raw_reads"))
	tr = float(pm.get_stat("Trimmed_reads"))
	pm.report_result("Alignment_rate", round(float(ar) * 100 / float(tr), 2))
	pm.report_result("Total_efficiency", round(float(ar) * 100 / float(rr), 2))
	mr = ngstk.count_multimapping_reads(out_tophat, args.paired_end)
	pm.report_result("Multimap_reads", mr)
	pm.report_result("Multimap_rate", round(float(mr) * 100 / float(tr), 2))

pm.run(cmd, re.sub(".bam$", "_sorted.bam", out_tophat), shell=False, follow=check_tophat)

pm.timestamp("### BAM to SAM sorting and indexing: ")

cmd = tools.samtools + " view -h " + out_tophat + " > " + out_tophat.replace(".bam", ".sam") + "\n"
cmd += tools.samtools + " sort --threads " + str(pm.cores) + " " + out_tophat + " -o " + out_tophat.replace(".bam", "_sorted.bam") + "\n"
cmd += tools.samtools + " index " + out_tophat.replace(".bam", "_sorted.bam") + "\n"
cmd += tools.samtools + " depth " + out_tophat.replace(".bam", "_sorted.bam") + " > " + out_tophat.replace(".bam", "_sorted.depth") + "\n"
pm.run(cmd, re.sub(".bam$", "_sorted.depth", out_tophat),shell=True)

pm.clean_add(out_tophat, conditional=False)
pm.clean_add(re.sub(".bam$" , ".sam", out_tophat), conditional=False)

if args.markDupl:
	pm.timestamp("### MarkDuplicates: ")

	aligned_file = re.sub(".sam$", "_sorted.bam",  out_tophat)
	out_file = re.sub(".sam$", "_dedup.bam", out_tophat)
	metrics_file = re.sub(".sam$", "_dedup.metrics", out_tophat)
	cmd = ngstk.markDuplicates(aligned_file, out_file, metrics_file)
	pm.run(cmd, out_file, follow= lambda:
		pm.report_result("Deduplicated_reads", ngstk.count_unique_mapped_reads(out_file, args.paired_end and not align_paired_as_single)))


# Create tracks
########################################################################################

pm.timestamp("### bam2wig: ")
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
cmd += " -r " + param.ESAT.refGen + args.genome_assembly + "_refGene.bed"
cmd += " > " + re.sub("_sorted.bam$", "_read_distribution.txt",trackFile)
pm.run(cmd, re.sub("_sorted.bam$", "_read_distribution.txt",trackFile),shell=True, nofail=True)

#pm.timestamp("### gene_coverage: ")
#cmd = tools.gene_coverage + " -i " + re.sub(".bam$" , ".bw",trackFile)
#cmd += " -r " + param.ESAT.refGen + args.genome_assembly + "_refGene.bed"
#cmd += " -o " + re.sub("_sorted.bam$", "",trackFile)
#pm.run(cmd, re.sub("_sorted.bam$", ".geneBodyCoverage.png",trackFile),shell=False)


# ESAT pipeline
########################################################################################
pm.timestamp("### End Sequence Analysis Toolkit (ESAT): ")

ESAT_folder = os.path.join(param.pipeline_outfolder,"ESAT_" + args.genome_assembly)
pm.make_sure_path_exists(ESAT_folder)
out_ESAT_gene = os.path.join(ESAT_folder,args.sample_name + ".gene.txt")
out_ESAT_window = os.path.join(ESAT_folder,args.sample_name + ".window.txt")

cmd = tools.java + " -Xmx" + str(pm.mem) + " -jar " + tools.ESAT
cmd += " -task " + str(param.ESAT.task)
cmd += " -in " + out_tophat
cmd += " -geneMapping " + param.ESAT.refGen + args.genome_assembly + "_refGene.tsv"
cmd += " -out " + args.sample_name
cmd += " -wLen " + str(param.ESAT.wLen)
cmd += " -wOlap " + str(param.ESAT.wOlap)
cmd += " -wExt " + str(param.ESAT.wExt)
cmd += " -sigTest " + str(param.ESAT.sigTest)
cmd += " -quality " + str(param.ESAT.quality)
cmd += " -multimap " + str(param.ESAT.multimap)

os.chdir(ESAT_folder)
pm.run(cmd, out_ESAT_gene, shell=False)


# Cleanup
########################################################################################

pm.stop_pipeline()
