#!/usr/bin/env python3.6

# Test run for GUESSmyLT

import subprocess
import os
import sys

def fill_path(file):
    path = os.path.realpath(__file__)
    tail = path.rsplit('/',1)
    return tail[0]+"/"+file

def main():
    # PATH to the reference FASTA file
    GENOME = "4.fa"

    # PATH to the read files
    READ_1 = "4_r1.fastq.gz"
    READ_2 = "4_r2.fastq.gz"

    # Name of ouput folder
    OUTPUT = os.getcwd() + "/GUESSmyLT_example_out/"

    # Type of organism
    ORGANISM = "eukaroyte"

    # Create the command
    command = "GUESSmyLT --reads "+fill_path(READ_1)+" "+fill_path(READ_2)+" --organism "+ORGANISM+" --reference "+fill_path(GENOME)+" --output "+OUTPUT
    print("Running the following command: "+command)

    # Execute the command
    os.system(command)


if __name__ == '__main__':
    main()
