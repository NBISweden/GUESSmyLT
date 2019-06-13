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

def infer_paired_region_with(genes):
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

        # Get reads mapped to a specific contig and in a sequence range
        # TODO: Look into optimizing this step, only take a subset (1000ish) reads? samfile.mate is not made for high throughput
        try:
            for read in samfile.fetch(contig, start, stop):
                if not read.mate_is_unmapped and read.is_read1:
                    first = read
                    second = samfile.mate(read)
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
        except ValueError as message:
            print(message)         

    return libs

def infer_paired_region_without():
    '''
    Function for inferring paired library-type by looking at a regions corresponding to genes
    '''
    # Counters for the different lib-types
    libs = {
        'fr': 0,
        'rf': 0,
        'ff': 0,
        'undecided0': 0
    }

    # Get reads 
    # TODO: Look into optimizing this step, only take a subset (10000ish) reads? samfile.mate is not made for high throughput
    try:
        max_read=20000
        count=0
        print ("As no gene was found we focus only on "+str(max_read)+" first reads (both reads mapped)")
        for read in samfile.fetch():       
            if not read.mate_is_unmapped and read.is_read1:
                count+=1
                if count > max_read:
                    break
                first = read
                second = samfile.mate(read)
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
                    if first.reference_start > second.reference_start:
                        # Flip order of reads
                        lib = lib[::-1]

                    libs[lib] += 1
                except: libs['undecided0'] +=1 #Some reads missing start or end-values               
    except ValueError as message:
        print(message)

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

        # Get reads mapped to a specific contig and in a sequence range
        for read in samfile.fetch(contig, start, stop):
            if not read.is_unmapped:
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

def write_result(lib_dict, state, genes):
    """
    Function for calculating precentages and writing results of inferring to resultfile.
    """
    lib_viz = {
        'fr_first': ["fr_firststrand", "3' ----------<==1== 5'","5' ==2==>---------- 3'", "inward"],
        'fr_second': ["fr_secondstrand", "3' ----------<==2== 5'", "5' ==1==>---------- 3'", "inward"],
        'fr': ["fr", "----------<=====","=====>----------", "inward"],
        'rf_first': ["rf_firststrand", "3' <==2==---------- 5'","5' ----------==1==> 3'", "outward"],
        'rf_second': ["rf_secondstrand", "3' <==1==---------- 5'","5' ----------==2==> 3'", "outward"],
        'rf': ["rf", "<=====----------","----------=====>" , "outward"],
        'ff_first': ["ff_firststrand", "3' <==2==----<==1== 5'","5' ---------------- 3'", "matching"],
        'ff_second': ["ff_secondstrand", "3' ---------------- 5'","5' ==1==>----==2==> 3'", "matching"],
        'ff': ["ff", "<=====----<=====","", "matching"],
        'f_first': ["f_firststrand", "3' ----------<==1== 5'","5' ---------------- 3'" ],
        'f_second': ["f_secondstrand", "3' ---------------- 5'","5' ==1==>---------- 3'" ],
        'r_first': ["r_firststrand", "3' ----------==1==> 5'","5' ---------------- 3'" ],
        'r_second': ["r_secondstrand", "3' ---------------- 5'","5' <==1==---------- 3'" ],
        'undecided0': ["undecided", "-------??-------","-------??-------", "NA"],
        'undecided': ["undecided", "3' -------??------- 5'","5' -------??------- 3'", "NA" ]
    }
    output = open(output_file, 'w+')
     # --------- Print GENERAL message ----------
    if state=='single' and not len(genes.features) > 0: #stop here nothing to do
        output.write("Sorry but without any gene and with single end reads, we cannot infer anything.")
        print ("Sorry but without any gene and with single end reads, we cannot infer anything.")
    else:
        output.write("Results of %s library inferred from reads: \n\n" % (state))
        print("Results of %s library inferred from reads: \n\n" % (state))

        # --------- Print HEADER ----------
        #paired vs single => relative orientation column
        # paired without gene vs paired with gene => 
        if 'fr' in lib_dict: # paired without gene Vizualization according to firststrand or not
            output.write("{:>15}    {:>20}    {:>10}    {:>8}    {:>}\n".format("Library type","Relative orientation","Reads","Percent","Vizualization"))
            print ('{:>15}    {:>20}    {:>10}    {:>8}    {:>}\n'.format("Library type","Relative orientation","Reads","Percent","Vizualization"))
        elif 'fr_first' in lib_dict: # paired with gene
            output.write("{:>15}    {:>20}    {:>10}    {:>8}    {:>20}\n".format("Library type","Relative orientation","Reads","Percent","Vizualization according to firststrand"))
            print ('{:>15}    {:>20}    {:>10}    {:>8}    {:>20}\n'.format("Library type","Relative orientation","Reads","Percent","Vizualization according to firststrand"))
        else:#single end
            output.write("{:>15}    {:>10}    {:>8}    {:>20}\n".format("Library type","Reads","Percent","Vizualization according to firststrand"))
            print("{:>15}    {:>10}    {:>8}    {:>20}\n".format("Library type","Reads","Percent","Vizualization according to firststrand"))


        # --------- Print VALUES ----------
        total_reads = 0
        for i in lib_dict:
            total_reads += lib_dict[i]

        for i in sorted(lib_dict):
            percent = "0%"
            if total_reads > 0:
                percent = '{:.1%}'.format(lib_dict[i]/total_reads)
            if 'fr' in lib_dict : #paired end
                output.write("{:>15}    {:>20}    {:>10}    {:>8}    {:>}\n{:>85}\n\n".format(lib_viz[i][0], lib_viz[i][3], lib_dict[i], percent, lib_viz[i][1], lib_viz[i][2]))
                print ('{:>15}    {:>20}    {:>10}    {:>8}    {:>}\n{:>85}\n\n'.format(lib_viz[i][0], lib_viz[i][3], lib_dict[i], percent, lib_viz[i][1], lib_viz[i][2]))
            elif 'fr_first' in lib_dict: #paired end
                output.write("{:>15}    {:>20}    {:>10}    {:>8}    {:>}\n{:>91}\n\n".format(lib_viz[i][0], lib_viz[i][3], lib_dict[i], percent, lib_viz[i][1], lib_viz[i][2]))
                print ('{:>15}    {:>20}    {:>10}    {:>8}    {:>}\n{:>91}\n\n'.format(lib_viz[i][0], lib_viz[i][3], lib_dict[i], percent, lib_viz[i][1], lib_viz[i][2]))
            else: #single end
                output.write("{:>15}    {:>10}    {:>8}    {:>20}\n{:>67}\n\n".format(lib_viz[i][0], lib_dict[i], percent, lib_viz[i][1], lib_viz[i][2]))
                print ('{:>15}    {:>10}    {:>8}    {:>20}\n{:>67}\n\n'.format(lib_viz[i][0], lib_dict[i], percent, lib_viz[i][1], lib_viz[i][2]))
                

        if 'fr' in lib_dict:
            output.write("Sorry but without annotated gene we cannot stand which strand is the referential (*_firststrand / *_secondstrand).\nYou cannot guess neither if the library is stranded or unstranded, we can only observe the relative orientation of reads between pairs.\n")
            print ("Sorry but without annotated gene we cannot stand which strand is the referential (*_firststrand / *_secondstrand).\nYou cannot guess neither if the library is stranded or unstranded, we can only observe the relative orientation of reads between pairs.\n")
        else:
            output.write("Roughly 50/50 split between the strands of the same library orientation should be interpreted as unstranded.")
            print("Roughly 50/50 split between the strands of the same library orientation should be interpreted as unstranded.")




# -----------------------------
# ---------- RUNNING ----------
print (sys.argv)
annotation = sys.argv[1]
bam = sys.argv[2]
read_path = sys.argv[3]
state = sys.argv[4]
run_name = sys.argv[5]
output_file = sys.argv[6]

samfile = pysam.AlignmentFile(bam, "rb")

print("Extracting genes...")
genes = extract_genes(annotation, run_name)
result=None

if state == 'single':
    print("Running single end reads inference")
    if len(genes.features) > 0:       
        result = infer_single_region(genes)
elif state == 'paired':
    print("Running paired end reads inference")
    if len(genes.features) > 0:
        result = infer_paired_region_with(genes)
    else:
        result = infer_paired_region_without()

print("Inference finished\n")
write_result(result, state, genes)
print("\nResults also written to file " + output_file )