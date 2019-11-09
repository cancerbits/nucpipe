#!/usr/bin/env python

## https://www.ncbi.nlm.nih.gov/pmc/articles/pmid/26890679/

"""
RNA NucSeq pipeline
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
pm = pypiper.PipelineManager(name = "rnaNucSeq", outfolder = outfolder, args = args)

# Tools
# pm.config.tools.scripts_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "tools")

# Resources
# pm.config.resources.ref_genome = os.path.join(pm.config.resources.genomes, args.genome_assembly)
# pm.config.resources.ref_genome_fasta = os.path.join(pm.config.resources.genomes, args.genome_assembly, args.genome_assembly + ".fa")
# pm.config.resources.chrom_sizes = os.path.join(pm.config.resources.genomes, args.genome_assembly, args.genome_assembly + ".chromSizes")


# Output
pm.config.parameters.pipeline_outfolder = outfolder

ngstk = pypiper.NGSTk(pm=pm)

tools = pm.config.tools
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

cmd = tools.java + " -Xmx" + str(pm.mem) + " -jar " + tools.trimmomatic

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

#cmd += " ILLUMINACLIP:"+resources.adapters+":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:60"
#cmd += " ILLUMINACLIP:"+resources.adapters+":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:10 MINLEN:25"
cmd += " ILLUMINACLIP:"+resources.adapters+":2:10:5 LEADING:1 TRAILING:1 SLIDINGWINDOW:4:1 MINLEN:21"

	# kallisto pipe:
	
	# cmd += " ILLUMINACLIP:" + resources.adapters + ":2:10:4:1:true"
	# # TODO: generalize the path to the adapters.
	# if args.quantseq: cmd += " ILLUMINACLIP:" + "/data/groups/lab_bsf/resources/trimmomatic_adapters/PolyA-SE.fa" + ":2:30:5:1:true"
	# cmd += " SLIDINGWINDOW:4:1"
	# cmd += " MAXINFO:16:0.40"
	# cmd += " MINLEN:21"
	
	# if args.quantseq: cmd += " HEADCROP:6"
	# cmd += " ILLUMINACLIP:" + resources.adapters + ":2:10:4:1:true"
	# if args.quantseq: cmd += " ILLUMINACLIP:" + "/data/groups/lab_bsf/resources/trimmomatic_adapters/PolyA-SE.fa" + ":2:30:5:1:true"
	# cmd += " SLIDINGWINDOW:4:1"
	# cmd += " MAXINFO:16:0.40"
	# cmd += " MINLEN:21"


trimmed_fastq = out_fastq_pre + "_R1_trimmed.fastq"
trimmed_fastq_R2 = out_fastq_pre + "_R2_trimmed.fastq"

pm.run(cmd, trimmed_fastq, 
	follow = ngstk.check_trim(trimmed_fastq, args.paired_end, trimmed_fastq_R2,
		fastqc_folder = os.path.join(param.pipeline_outfolder, "fastqc/")))




# RSEM pipeline.
########################################################################################
pm.timestamp("### RSEM alignment: ")

indexed_genome = os.path.join(pm.config.resources.genomes, args.genome_assembly, "indexed_rsem", args.genome_assembly)
if not (args.ERCC_mix == "False" ):
	indexed_genome = os.path.join(pm.config.resources.genomes, args.genome_assembly, "indexed_rsem_"+args.ERCC_assembly, args.genome_assembly)


cmd = tools.rsem
cmd += " --bowtie2 --calc-pme --calc-ci --estimate-rspd --time  --tag MA:i:2 -fragment-length-min 1 --fragment-length-max 500 --output-genome-bam"
cmd += " -p " + str(pm.cores)

rsem_folder = os.path.join(param.pipeline_outfolder,"rsem_" + args.genome_assembly)
pm.make_sure_path_exists(rsem_folder)
out_rsem = os.path.join(rsem_folder, args.sample_name + ".transcript.bam")

if not args.paired_end:
	cmd += " " + out_fastq_pre + "_R1_trimmed.fastq "
else:
	cmd += " --paired-end " + out_fastq_pre + "_R1_trimmed.fastq " + out_fastq_pre + "_R2_trimmed.fastq "

cmd += indexed_genome + " " + rsem_folder + "/" + args.sample_name

pm.run(cmd, out_rsem,
	follow=lambda: pm.report_result("Aligned_reads", ngstk.count_unique_mapped_reads(out_rsem, args.paired_end)))


# Cleanup
########################################################################################
# remove temporary marker file:
pm.stop_pipeline()





