#!/usr/bin/env python3.6
# Top script for pipeline of GUESSmyLT.
# The script handles user arguments, validates inputs and executes the pipeline by calling Snakemake.
#
# Example runs:
#   With minimum inputs: python GUESSmyLT.py --reads read1 read2 --organism euk
#   With some optional inputs: python GUESSmyLT.py --reads read1 read2 --organism euk --threads 10 --memory 5G --reference reference.fasta
#
# For more options type: python GUESSmyLT.py --help

import argparse, sys, json, os, re, subprocess
from BCBio import GFF

# The path to the directory where the script is executed.
script_dir = os.path.dirname(os.path.abspath(__file__)) + "/"
# Working dir, where files will be written if no --ouput arg given
working_dir = os.getcwd() + "/"

# ---Functions for checking user arguments---

def check_organism():
    """
    Checks that organism argument is correct.
    Can be any string within prokaryote or eukaryote.
    For example, euk, euka, eu are all valid arguments.
    """
    if args.organism.lower() in "prokaryote" and args.organism.lower() not in "eukaryote" :
        args.organism="prokaryote"

    elif args.organism.lower() in "eukaryote" and args.organism.lower() not in "prokaryote":
        args.organism="eukaryote"
    else:
        print("Error. Unrecognized organism --organism "+args.organism+". Only eukaryote/prokaryote are valid organism.")
        sys.exit()

def check_annotation(annotation_file):
    """
    Checks that an optional annotation file is correct by containing genes.
    """
    in_handle=open(annotation_file)
    for record in GFF.parse(in_handle, limit_info=dict(gff_type=["gene"])):
        if record:
            return True
    # If no genes are found, gff file cannot be used in analysis.
    print("Error. No genes could be found in --annotation "+annotation_file+". Please, submit a .gff file containing genes or no .gff file.")
    sys.exit()



def check_subsample(num):
    """
    Checks that number of reads used for subsampling is even.
    Otherwise, it will not work if we work with paired end data.
    """
    if num%2==0:
        # Number is even --> Valid
        pass
    else:
        # Number is odd --> Not valid
        print("Error. Number of reads used for subsampling (--subsample) must be even.")
        sys.exit()



def check_mapped():
    # Not developed yet. Added in to do list.
    print("Checker for mapper has not been developed yet.")
    print("Right now you cannot provide a map-file (.bam).")
    print("Therefore, skip the map-file and let GUESSmyLT do the mapping for you.")
    sys.exit()


def check_reference(ref):
    """
    Checks that provided reference is valid according to:
        1. Can be opened.
        2. Has at least one line that begins with '>'.

    Also checks if we are dealing with genome or transcriptome by looking at
    the number of lines that begin with '>'.
    If one line begins with '>' we are dealing with genome.
    If more than one line begin with '>' we are dealing with transcriptome.
    """
    try:
        open(ref)
    except:
        print("Error. Cannot open --reference '"+ref+"'. Make sure it exists.")
        sys.exit()
    if ".gz" in ref:
        zip_command="zcat <"
    else:
        zip_command="cat"
    num_headers=subprocess.check_output(
    zip_command+" "+ref+" | grep '^>' | wc -l",
    shell=True,
    encoding="utf8")
    if int(num_headers.split("\n")[0]) > 1:
        # Multiple headers --> dealing with transcriptome.
        # So far we return genome anyway because we have not optimized busco for transcriptome.
        # Need to optimize busco script!
        # This is however a later improvement since we can lie to busco that we are dealing with genome.
        #return "transcriptome"
        return "genome"
    elif int(num_headers.split("\n")[0])==1:
        # One header --> dealing with genome.
        return "genome"
    else:
        print("Error. Reference file is not in fasta format. Missing '>' in beginning of fasta header.")
        sys.exit()


def check_memory(args):
    """
    Checks that the maximum memory used are written correctly in gigabytes, e.g. 8G.
    """
    if "G" == args.memory[-1] and "0" != args.memory[0] and args.memory[:-1].isdigit():
        return True
    else:
        print("Invalid --memory argument: '"+args.memory+"'. Memory should be given in GIGABYTES. For example --memory '4G'")
        sys.exit()


# ---Main from here ---

def main():

    # Use argparse for handling user arguments.
    # Mandatory arguments are:
    # --organism and --reads
    # Optional arguments are:
    # --reference, --annotation, --mapped (For providing ref, ann, map files).
    # --threads, --memory (How many threads and maximum memory that GUESSmyLT can use.)

    parser=argparse.ArgumentParser(description="GUESSmyLT, GUESS my Library Type. Can predict the library type used for RNA-Seq. The prediction is based on the orientaion of your read files in .fastq format. Knowing the library type helps you with downstream analyses since it greatly improves the assembly.")
    parser._action_groups.pop()
    required=parser.add_argument_group("REQUIRED ARGUMENTS")
    required.add_argument("--organism",type=str, help="What organism are you dealing with? prokaryote or eukaryote.")
    optional=parser.add_argument_group("OPTIONAL ARGUMENTS")
    required.add_argument("--reads",nargs="+",type=str,help="One or two read files in .fastq format. Files can be compressed or uncrompressed. Handles interleaved read files and any known .fastq header format. ")
    optional.add_argument("--subsample",type=int,default=100000,help="Number of subsampled reads that will be used for analysis. Must be an even number. Default value is 100,000 reads.")
    optional.add_argument("--reference",type=str,help="Reference file in .fasta format. Reference can be either transcriptome or genome.")
    optional.add_argument("--annotation",type=str,help="Annotation file in .gff format. Needs to contain genes. Not implemented yet.")
    optional.add_argument("--mapped",type=str,help="Mapped file in sorted .bam format. Reference that reads have been mapped to has to be provided. Checker for this has not been developed yet. Therefore, do not provide a mapping file.") #Maybe add this?
    optional.add_argument("--threads", type=int, default=10,help="The number of threads that can be used by GUESSmyLT. Needs to be an integer. Defualt value is 10.")
    optional.add_argument("--memory",type=str, default="8G",help="Maximum memory that can be used by GUESSmyLT in GB. E.g. '10G'. Default value is 8G.")
    optional.add_argument("--output",type=str,default=working_dir,help="Full path to output directory. Default is working directory.")
    args=parser.parse_args()

    global output_dir
    output_dir = args.output
    if not os.path.exists(args.output):
        os.system("mkdir "+args.output)
    if not os.path.exists(args.output+"intermediate_data"):
        os.system("mkdir "+args.output+"intermediate_data")

    # Copy the snakemake config file from the script directory (install or git clone dir)
    if not os.path.exists(args.output+"config.json"):
        os.system("cp "+script_dir+"config.json "+args.output+"config.json")
    # ---Check subsampling---
    check_subsample(args.subsample)

    # ---Check reads---
    # At least one readfile must be provided. (single end or interleaved paired end reads)
    # At most two readiles can be provided. (paired end reads)
    readname1=""
    readname2=""
    subsampled_names = []
    if not args.reads:
        print("Error. No read files provided. At least one read file should be provided.")
        sys.exit()
    elif len(args.reads) == 1:
        readname1=re.split(".fq|.fastq",os.path.basename(args.reads[0]))[0]+".sub."+str(int(args.subsample))
        subsampled_names.append(args.output+"intermediate_data/"+readname1+".fastq.gz")
    elif len(args.reads) == 2:
        readname1=re.split(".fq|.fastq",os.path.basename(args.reads[0]))[0]+".sub."+str(int(args.subsample))
        subsampled_names.append(args.output+"intermediate_data/"+readname1+".fastq.gz")
        readname2=re.split(".fq|.fastq",os.path.basename(args.reads[1]))[0]+".sub."+str(int(args.subsample))
        subsampled_names.append(args.output+"intermediate_data/"+readname2+".fastq.gz")
    elif len(args.reads) > 2:
        print("Error. Too many read files. Only one or two read files can be provided.")
        sys.exit()

    # ---Check mapping file---
    # If a mapping file (.bam/.sam) is provided, the reference used for mapping
    # must also be provided.
    if args.mapped and not args.reference:
        print("Error. If a mapping file is provided, the reference used for mapping must also be provided.")
        sys.exit()

    # Tells busco what type of reference: genome or transcriptome.
    # Right now we always tell busco that it's working with a genome,
    # because it gives less result files if it is a transcriptome.
    # This could be optimized.
    busco_reference_mode="genome"
    if args.reference:
        busco_reference_mode=check_reference(args.reference)
    else:
        # Add Trinity assembly to config file.
        args.reference="intermediate_data/"+readname1+".fasta"

    # Get refname from the path, basis for name of BUSCO run
    refname=re.split(".fa|.fasta",os.path.basename(args.reference))[0]

    if args.mapped:
        check_mapped()
    else:
        args.mapped="intermediate_data/"+readname1+"_on_ref_"+re.split('/|\.',args.reference)[-2]+"_sorted.bam"

    if args.annotation:
        print("Inputting annotation files not implemeted yet. Exiting...")
        sys.exit()
        check_annotation(args.annotation)
    else:
        args.annotation="intermediate_data/run_"+refname

    check_memory(args)

    # Update snakemake config file
    config_path = args.output+"config.json"
    with open(config_path,"r+") as configfile:
        data=json.load(configfile)
        data["trinity"]["reference"]=args.reference
        data["input"]["to_subsample"]=args.reads
        data["input"]["reads"]=subsampled_names
        data["input"]["readname"]=readname1
        data["subsample"]=int(args.subsample)
        data["busco"]["annotation"]=args.annotation
        data["input"]["organism"]=args.organism
        data["bowtie2"]["mapped-reads"]=args.mapped
        data["busco"]["lineage"]=data["prokaryote_db"] if args.organism == "prokaryote" else data["eukaryote_db"]
        data["busco"]["mode"]=busco_reference_mode
        data["input"]["threads"]=args.threads
        data["input"]["memory"]=args.memory
        data["output"]=args.output
        data["script_dir"]=script_dir
        configfile.seek(0)
        json.dump(data,configfile,indent=4)
        configfile.truncate()

    # Execute Snakemake
    os.system("snakemake -s "+script_dir+"Snakefile -d "+args.output+" result_"+readname1+"_on_"+refname+".txt --cores "+str(args.threads))

if __name__ == "__main__":
    main()
