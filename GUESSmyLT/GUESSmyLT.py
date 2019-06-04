#!/usr/bin/env python3
# Top script for pipeline of GUESSmyLT.
# The script handles user arguments, validates inputs and executes the pipeline by calling Snakemake.
#
# Example runs:
#   With minimum inputs: python GUESSmyLT.py --reads read1 read2 --organism euk
#   With some optional inputs: python GUESSmyLT.py --reads read1 read2 --organism euk --threads 10 --memory 5G --reference reference.fasta
#
# For more options type: python GUESSmyLT.py --help

import argparse, sys, json, os, re, subprocess, shlex
from BCBio import GFF
import pysam
import tarfile
from Bio import SeqIO
from shutil import copyfile

# The path to the directory where the script is executed.
script_dir = os.path.dirname(os.path.abspath(__file__)) + "/"
# Working dir, where files will be written if no --ouput arg given
working_dir = os.getcwd() + "/"

# --- check BUSCO dataset ---

# Set busco eukaryote 
busco_euk = script_dir+"data/eukaryota_odb9"
if not os.path.isfile(busco_euk):
    tf = tarfile.open(script_dir+"data/eukaryota_odb9.tar.gz")
    tf.extractall(path=script_dir+"data/")
# Set busco prokaryote
busco_prok = script_dir+"data/bacteria_odb9"
if not os.path.isfile(busco_prok):
    tf = tarfile.open(script_dir+"data/bacteria_odb9.tar.gz")
    tf.extractall(path=script_dir+"data/")


# --- Miscellaneous functions ---

def is_interleaved(read,pattern1,pattern2):
    """
    Function for checking if a single read file is interleaved.
    If both pattern1 and pattern2 can be found in the first 10000 lines,
    True is returned. Else False.
    If two read headers are identical, True is returned. Else false.
    """
    if read.endswith(".gz"):
        f=gzip.open(read,"rt")
    else:
        f=open(read,"r")

    num_read1=0
    num_read2=0
    counter=0
    headers=[]

    for line in f:
        # Only look at read headers. They come every fourth line.
        if counter%4==0:
            if re.match(pattern1,line):
                num_read1+=1
            elif re.match(pattern2,line):
                num_read2+=1

            if line.startswith("@"):
                ID=line.split(" ")[0]
                if ID in headers:
                    return True
                else:
                    headers.append(ID)
        counter+=1
        if counter==10000:
            break
    f.close()

    if num_read1>0 and num_read2>0:
        return True
    else:
        return False

def get_full_path(path):
    #create absolute path to it
    if not path.startswith("/"):
        if path.startswith("./"):
            return ref[2:]

        return working_dir+path
    else:
        return path

def get_type(bam):
    cmd = "samtools view -c -f 1 " + str(bam)
    proc = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = proc.communicate()
    nb_reads_paired = int(str(out.decode("utf-8")).split("\n")[0])
    if nb_reads_paired  > 0:
        return "paired"
    else:
        return "single"


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
    #check first if we have a value
    if num:
        if num%2==0:
            # Number is even --> Valid
            return num
        else:
            # Number is odd --> Not valid
            print("Error. Number of reads used for subsampling (--subsample) must be even.")
            sys.exit()
    else:
        return "full"


def check_mapped(bam_file, ref):
    # check if mapped file provided contains same sequence ID as the reference provided
    pysam.index(bam_file)
    samfile = pysam.AlignmentFile(bam_file, "rb")
    record_dict = SeqIO.to_dict(SeqIO.parse(open(ref),'fasta'))
    read_found = None   

    for record_id in record_dict:
        try:
            for read in samfile.fetch(record_id):
                read_found=1
                break          
            samfile.close()
            
            if read_found:
                print("Reference and bam have Sequence ID matching. We can continue.")
                break
        except:
            pass

    if not read_found:
        print("Not match found between sequence id from the bam file and the reference file. Please check provided files.")
        sys.exit()

def check_mode(mode):
    """
    checks if we are dealing with genome or transcriptome by looking at the provided option    
    """
    #check that the mode asked exists
    if not mode:
        print("Mode parameter (--mode) must be filled when no annotation provided (--annotation). Accepted value: genome or transcriptome.")
        sys.exit()
    elif mode.lower() in "genome" and mode.lower() not in "transcriptome":
        return "genome"
    elif mode.lower() in "transcriptome" and mode.lower() not in "genome" :
        print("Transcriptome mode is not yet implemented.")
        sys.exit()
        return "transcriptome"
    else:
        print("Unrecognized mode --mode "+mode+". Only <genome> or <transcriptome> are valid mode.")
        sys.exit()

def check_reference(ref):
    """
    Checks that provided reference is valid according to:
        1. Can be opened.
        2. Has at least one line that begins with '>'.

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
    shell=True)

    if int(str(num_headers.decode("utf-8")).split("\n")[0]) == 0:
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

    parser=argparse.ArgumentParser(description="GUESSmyLT, GUESS my Library Type. Can predict the library type used for RNA-Seq. The prediction is based on the orientaion of your read file(s) in .fastq/.fastq.gz/.bam format. Knowing the library type helps you with downstream analyses since it greatly improves the assembly.")
    parser.add_argument("--organism",type=str, help="Mandatory when no annotation provided. What organism are you dealing with? prokaryote or eukaryote.")
    parser.add_argument("--reads",nargs="+",type=str,help="One or two read files in .fastq or .fastq.gz format. Files can be compressed or uncrompressed. Handles interleaved read files and any known .fastq header format. ")
    parser.add_argument("--subsample", type=int,help="Number of subsampled reads that will be used for analysis. Must be an even number.")
    parser.add_argument("--reference", type=str,help="Mandatory when --mapped used or when no reads provided (--reads). Reference file in .fa/.fasta format. Reference can be either transcriptome or genome.")
    parser.add_argument("--mode", default="genome", type=str,help="Mode can be genome or transcriptome (default genome). It defines how the reference fasta file will be handled by BUSCO. This option is used when no annotation is provided (--annotation).")
    parser.add_argument("--annotation", type=str,help="Annotation file in .gff format. Needs to contain genes.")
    parser.add_argument("--mapped",type=str,help="Mapped file in .bam format (Will be sorted). Reference that reads have been mapped to has to be provided.")
    parser.add_argument("--threads", type=int, default=2,help="The number of threads that can be used by GUESSmyLT. Needs to be an integer. Defualt value is 2.")
    parser.add_argument("--memory",type=str, default="8G",help="Maximum memory that can be used by GUESSmyLT in GB. E.g. '10G'. Default value is 8G.")
    parser.add_argument("--output",type=str,default=working_dir,help="Full path to output directory. Default is working directory.")
    parser.add_argument("-n", action="store_true", help="(Snakemake dryrun option) Allows to see the scheduling plan including the assigned priorities.")
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args=parser.parse_args()

    # ---fill the full path---
    if args.reads:
        args.reads[0] = get_full_path(args.reads[0])
        if len(args.reads) > 1:
            args.reads[1] = get_full_path(args.reads[1])
    if args.reference:
        args.reference = get_full_path(args.reference)
    if args.mapped:
        args.mapped = get_full_path(args.mapped)    
    if args.annotation:
        args.annotation = get_full_path(args.annotation)      

    # ---Deal with output dir---
    global output_dir
    if not args.output.endswith("/"):
       args.output=args.output+"/" 
    output_dir =  get_full_path(args.output)
    if not os.path.exists(output_dir):
        os.system( "mkdir " + output_dir)
    if not os.path.exists( output_dir + "intermediate_data" ):
        os.system( "mkdir "+ output_dir + "intermediate_data" )

    # ---Check subsampling---
    subsample_value=check_subsample(args.subsample)

        # ---Get reads name and set sample name---
    sample_name = None
    reads_names = []
    if args.reads:
        read1_name=re.split(".fq|.fastq",os.path.basename(args.reads[0]))[0]
        reads_names=[read1_name]
        sample_name = read1_name
        if len(args.reads) > 1:
            read2_name=re.split(".fq|.fastq",os.path.basename(args.reads[1]))[0]
            reads_names.append(read2_name)
    else:
        sample_name = re.split(".sam|.bam",os.path.basename(args.mapped))[0]


    # ---Check mapping file---
    # If a mapping file (.bam/.sam) is provided, the reference used for mapping
    # must also be provided. If not we set the expected mapped file that will be produce by the pipeline
    mapped_name = None
    if args.mapped:
        if not args.reference:
            print("Error. If a mapping file is provided, the reference used for mapping must also be provided.")
            sys.exit()
        else:
            check_mapped(args.mapped, args.reference)
            
        if get_type(args.mapped) == "single":
            if not args.reads:
                print ("The provided mapped reads file is Single End type. When providing such input, we need the fastq file to correctly guess the library type. Why ? Here the explanation:")
                print ("""
                    Inferring the library type of single end data turned out to be more challenging than paired end
                    data. Just as the lacking information of single end compared to paired end impairs analysis of
                    RNA-Seq data it makes predicting the library type harder. The information about which strand
                    the reads were generated from is in the case of paired end gained from the relative orientation of
                    paired reads. In the case of single end reads this information has to be sourced from looking at
                    the original reads. The library type is predicted by looking at the mapped direction of the reads and whether the mapped
                    reads have been reverse complemented. This is due to the fact that Bowtie2 saves reads that map
                    in reverse direction as the corresponding sequence in the reference, not the original reads sequence.
                """)
                sys.exit()
            elif len(args.reads) >= 2:
                print ("The provided mapped reads file is Single End type. When providing such input, we expect only one fastq file !") 
                sys.exit()
            elif is_interleaved(args.reads[0],"^@.+/1","^@.+/2") or is_interleaved(args.reads[0],"^@.+ 1:","^@.+ 2:"):
                print ("The provided mapped reads file is Single End type, while the fastq file is paired end (interleaved one). Please provide proper files.") 
                sys.exit()
        mapped_name=re.split(".bam|.sam",os.path.basename(args.mapped))[0]        

    else:
        # ---Check reads---
        # At least one readfile must be provided. (single end or interleaved paired end reads)
        # At most two readiles can be provided. (paired end reads)
        if not args.reads:
            print("Error. No read file(s) (.fastq/.fastq.gz) or mapped file (.bam) provided.")
            sys.exit()
        elif len(args.reads) > 2:
            print("Error. Too many read files. Only one or two read files can be provided.")
            sys.exit()
        
    # ---check memory---
    check_memory(args)
        
    # check reference
    if args.reference:    
        check_reference(args.reference)

    # check BUSCO information i.e: type of reference: genome or transcriptome.
    busco_reference_mode = None
    refname = None
    lineage = None
    if not args.annotation:
        busco_reference_mode=check_mode(args.mode)

        # Get refname from the path, basis for name of BUSCO run
        if args.reference:
            refname=re.split(".fa|.fasta",os.path.basename(args.reference))[0]
        else:
            refname = sample_name

        # set busco dataset
        if args.organism == "prokaryote":
            lineage = busco_prok
        else:
            lineage = busco_euk

    # check annotation file
    if args.annotation:
        print("Inputting annotation files not implemeted yet. Exiting...")
        check_annotation(args.annotation)

    # Update snakemake config file
    config_path = args.output+"config.json"
    with open(config_path,"w") as configfile:
        data = {
                "reference_path" : args.reference,
                "reference_name" : refname,
                "reads" : args.reads,
                "reads_names" : reads_names,
                "name" : sample_name,
                "mapped" : args.mapped,
                "mapped_name" : mapped_name,
                "annotation" : args.annotation,
                "organism" : args.organism,
                "lineage" : lineage,
                "mode" : busco_reference_mode,
                "subsample" : subsample_value,
                "threads" : args.threads,
                "memory" : args.memory,
                "output" : output_dir,
                "script_dir" : script_dir
               }
        configfile.seek(0)
        json.dump(data,configfile,indent=4)
        configfile.truncate()

    # Execute Snakemake
    if args.n:
        os.system("snakemake -nps "+script_dir+"Snakefile -d "+args.output+" result_"+sample_name+"_on_"+refname+".txt --cores "+str(args.threads))
    else:
        os.system("snakemake -s "+script_dir+"Snakefile -d "+args.output+" result_"+sample_name+"_on_"+refname+".txt --cores "+str(args.threads))


if __name__ == "__main__":
    main()
