#!/usr/bin/env python3 
"""
Script for inferring library type of mapped reads.
Uses augustus-predicted genes via BUSCO and bowtie2-mapped reads.
$ python inferr_lib.py {reference} {bam} {readname} {readpath} {PE_SE} {single OR paired}
"""

import sys
import os  
import re
import fnmatch
from difflib import SequenceMatcher
from BCBio import GFF
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import pysam
import gzip

# ---------- FUNCTIONS ----------

def extract_genes(annotation, run_name):
    '''
    Function for extracting genes corresponding to BUSCO hits (genes).
    Returns a SeqRecord object callesd correct_genes with one feature per gene.
    '''
    #where correct genes will be saved
    correct_genes = SeqRecord(seq='', id='correct_genes')
    limit_infos = dict( gff_type = ["gene"]) # Only want genes

    #################################################
    # if annotation is a folder data comes from busco #
    if os.path.isdir(annotation):
        file_tsv = open(annotation+"/full_table_"+run_name+".tsv", 'r')

        # Extract BUSCO IDs, start and end from table of hits into SeqRecord object, each BUSCO as a SeqFeature
        busco_record = SeqRecord(seq='', id='hits')
        hit_found = None
        for line in file_tsv.readlines():
            hit = (re.search(r'(\S*)\s(Complete)\s(\S*)\s(\S*)\s(\S*)\s\S*', line))
            if hit:
                hit_found = 1
                busco_record.features.append(SeqFeature(FeatureLocation(int(hit.group(4)), int(hit.group(5))), id=hit.group(1), type='gene', qualifiers={'contig': hit.group(3)}))
        file_tsv.close()
        # No Complete BUSCO found let's try to find a least one duplicated one /!\take only one of the duplicated
        if not hit_found:
            print("No Complete BUSCO found, let's look for Duplicated")
            file_tsv = open(annotation+"/full_table_"+run_name+".tsv", 'r')
            known_busco_ids={}
            for line in file_tsv.readlines():
                columns = line.split('\t') 

                if len(columns) > 1 and columns[1].lower() == "duplicated":
                    if not columns[0] in known_busco_ids:
                        known_busco_ids[columns[0]]=1
                        #print("add BUSCO: "+columns[0]+"\n")
                        busco_record.features.append(SeqFeature(FeatureLocation(int(columns[3]), int(columns[4])), id=columns[0], type='gene', qualifiers={'contig': columns[2]}))
            file_tsv.close()

        # Match the BUSCOs to augustus predicted genes in gff files
        gff_records = []

        # Extract genes from .gffs
        if hit_found:           
            for busco in busco_record.features:
                filename = busco.id+".gff" # gff filenames are [busco_id].gff
                try:
                    file_gff = open(annotation+"/augustus_output/gffs/"+filename)
                    for record in GFF.parse(file_gff, limit_info=limit_infos):
                        gff_records.append(record)
                    file_gff.close()
                except:
                    pass
        else: # No Complete BUSCO found let's try to find a least one duplicated one /!\take only one of the duplicated
            for busco in busco_record.features:
                #print("busco: "+busco.id)
                for filename in fnmatch.filter(os.listdir(annotation+"/augustus_output/predicted_genes/"), busco.id+"*"):                    
                    #print("filename: "+filename)
                    try:
                        file_gff = open(annotation+"/augustus_output/predicted_genes/"+filename)
                        for record in GFF.parse(file_gff, limit_info=limit_infos):
                            #print("add a record "+ record.id)
                            gff_records.append(record)
                        file_gff.close()
                    except:
                        pass

        # Find augustus predicted genes from the gff that match BUSCOs
        known_busco_ids={}
        for hit in busco_record.features:
            if not busco_record.id in known_busco_ids:
                for rec in gff_records:
                   for feature in rec.features:
                        if hit.location.start-1 == feature.location.start and hit.location.end == feature.location.end: # For some reason start has 1 nt diff
                            feature.id = rec.id
                            correct_genes.features.append(feature)
                            known_busco_ids[busco_record.id]=1
                            break
                   if busco_record.id in known_busco_ids:
                        break


    else:
        for record in GFF.parse(annotation, limit_info=limit_infos):
            for feature in record.features:
                feature.id = record.id
                correct_genes.features.append(feature)

    print("Number of genes extracted: %d\n" % (len(correct_genes.features)))
    return correct_genes

def infer_paired_region(genes):
    '''
    Function for inferring paired library-type by looking at a regions corresponding to genes
    '''
    # Counters for the different lib-types
    libs = {
        'fr_first': 0,
        'fr_second': 0,
        'rf_first': 0,
        'rf_second': 0,
        'ff_first': 0,
        'ff_second': 0,
        'undecided': 0
    }

    # For every gene extract reads that map to gene region and find lib-type
    for gene in genes.features:
        contig = gene.id
        start = int(gene.location.start)
        stop = int(gene.location.end)
        strand = gene.strand
        reads = []
       
        # Get reads mapped to a specific contig and in a sequence range
        # TODO: Look into optimizing this step, only take a subset (1000ish) reads? samfile.mate is not made for high throughput
        try:
            for read in samfile.fetch(contig, start, stop):
                if not read.mate_is_unmapped and read.is_read1:
                    reads.append([read, samfile.mate(read)])
        except ValueError as message:
            print(message)

        # Check lib-type of reads
        for read in reads:
            first = read[0]
            second = read[1]
            try:
                lib = ''
                if not first.is_reverse:
                    lib += 'f'
                else:
                    lib += 'r'
                if not second.is_reverse:
                    lib += 'f'
                else:
                    lib += 'r'
                # Gene on sense strand
                if strand == 1:
                    if first.reference_start > second.reference_start:
                        # Flip order of reads
                        lib = lib[::-1]
                        lib += '_first'
                    elif first.reference_start < second.reference_start:
                        lib += '_second'
                    else:
                        lib = 'undecided'
                # Gene on antisense
                elif strand == -1:
                    if first.reference_start > second.reference_start:
                        # Flip order of reads
                        lib = lib[::-1]
                        lib += '_second'
                    elif first.reference_start < second.reference_start:
                        lib += '_first'
                    else:
                        lib = 'undecided'
                libs[lib] += 1
            except: libs['undecided'] +=1 #Some reads missing start or end-values
    return libs

def infer_single_region(genes):
    """
    Function for inferring library type of single-ended library types
    r-library prediction not validated
    """

    libs = {
        'f_first': 0,
        'f_second': 0,
        'r_first': 0,
        'r_second': 0,
        'undecided': 0
    }

    # Read original read IDs and sequences from the fastq files
    og_reads = {}
    with gzip.open(read_path, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            og_reads[record.id]= str(record.seq)

    # For every gene extract reads that map to gene region and find lib-type
    for gene in genes.features:
        contig = gene.id
        start = int(gene.location.start)
        stop = int(gene.location.end)
        strand = gene.strand
        reads = []
        # Get reads mapped to a specific contig and in a sequence range
        for read in samfile.fetch(contig, start, stop):
            if not read.is_unmapped:
                reads.append(read)

        # Check lib-type of reads
        for read in reads:
            try:
                flag = SequenceMatcher(None, og_reads[read.query_name], read.query_sequence).ratio() >= 0.8
                lib = ''
                if strand == 1:
                    if flag and not read.is_reverse:
                        lib += 'f_second'
                    elif not flag and read.is_reverse:
                        lib += 'f_first'
                    elif flag and read.is_reverse:
                        lib += 'r_second'
                    elif not flag and not read.is_reverse:
                        lib += 'r_first'
                    else:
                        lib = 'undecided'
                elif strand == -1:
                    if not flag and read.is_reverse:
                        lib += 'f_second'
                    elif flag and  not read.is_reverse:
                        lib += 'f_first'
                    elif flag and read.is_reverse:
                        lib += 'r_first'
                    elif not flag and  not read.is_reverse:
                        lib += 'r_second'
                    else:
                        lib = 'undecided'
                else:
                    lib = 'undecided'
                libs[lib] += 1
            except: libs['undecided'] +=1 # Some reads missing start or end-values
    return libs

def write_result(lib_dict):
    """
    Function for calculating precentages and writing results of inferring to resultfile.
    """
    lib_viz = {
        'fr_first': ["3' ----------<==1== 5'","5' ==2==>---------- 3'" ],
        'fr_second': ["3' ----------<==2== 5'", "5' ==1==>---------- 3'"],
        'rf_first': ["3' ----------==1==> 5'","5' <==2==---------- 3'" ],
        'rf_second': ["3' ----------==2==> 5'","5' <==1==---------- 3'" ],
        'ff_first': ["3' ----------<==1== 5'","5' <==2==---------- 3'" ],
        'ff_second': ["3' ----------==2==> 5'","5' ==1==>---------- 3'" ],
        'f_first': ["3' ----------<==1== 5'","5' ---------------- 3'" ],
        'f_second': ["3' ---------------- 5'","5' ==1==>---------- 3'" ],
        'r_first': ["3' ----------==1==> 5'","5' ---------------- 3'" ],
        'r_second': ["3' ---------------- 5'","5' <==1==---------- 3'" ],
        'undecided': ["3' -------??------- 5'","5' -------??------- 3'" ]
    }
    output = open(output_file, 'w+')
    output.write("Results of %s library inferring of reads %s on ref %s: \n\nLibrary type    Reads     Percent     Vizualization according to firststrand\n" % (state, run_name, reference))
    print("Results of %s library inferring of reads %s on ref %s: \n\nLibrary type    Reads     Percent     Vizualization according to firststrand\n" % (state, run_name, reference))

    total_reads = 0
    for i in lib_dict:
        total_reads += lib_dict[i]
    max_len = max(8, len(str(total_reads)))
    max_len2 = 52+max_len
    if total_reads > 0:
        for i in sorted(lib_dict):
            percent = '{:.1%}'.format(lib_dict[i]/total_reads)
            output.write("%12s %*d %11s %26s\n%*s\n\n" % (i, max_len, lib_dict[i], percent, lib_viz[i][0], max_len2, lib_viz[i][1]))
            print("%12s %*d %11s %26s\n%*s\n\n" % (i, max_len, lib_dict[i], percent, lib_viz[i][0], max_len2, lib_viz[i][1]))
    else:
        for i in sorted(lib_dict):
            output.write("%12s %*d %11s %26s\n%*s\n\n" % (i, max_len, 0, 0, lib_viz[i][0], max_len2, lib_viz[i][1]))
            print("%12s %*d %11s %26s\n%*s\n\n" % (i, max_len, 0, 0, lib_viz[i][0], max_len2, lib_viz[i][1]))
    output.write("Roughly 50/50 split between the strands of the same library orientation should be interpreted as unstranded.")
    print("Roughly 50/50 split between the strands of the same library orientation should be interpreted as unstranded.")

# ---------- RUNNING ----------
print (sys.argv)
annotation = sys.argv[1]
bam = sys.argv[2]
read_path = sys.argv[3]
state = sys.argv[4]
run_name = sys.argv[5]
reference = sys.argv[6]
output_file = sys.argv[7]

samfile = pysam.AlignmentFile(bam, "rb")

print("Extracting genes...")
genes = extract_genes(annotation, reference)

if state == 'single':
    print("Running single end inferring")
    result = infer_single_region(genes)
elif state == 'paired':
    print("Running paired end inferring")
    result = infer_paired_region(genes)

print("Prediction finished:\n")
write_result(result)
print("\nResults also written to file " + output_file )
