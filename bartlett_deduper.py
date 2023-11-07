#!/usr/bin/env python

#-----------------------------------------------------------------------------------------------------------------------------
#                                            Implement Argparse
#-----------------------------------------------------------------------------------------------------------------------------

import re
import argparse

def get_args():
    parser = argparse.ArgumentParser(description="This python code takes a sorted SAM file, and outputs a sorted sam file without PCR duplicates")
    parser.add_argument("-f", "--file", help="Location of the SAM file to deduplicate, MSUT BE SORTED", required=True)
    parser.add_argument("-o", "--out", help="The output file path/name")
    parser.add_argument("-u", "--umi", help="The umi file path/name")

    return parser.parse_args()

#Global Variables
args=get_args()
f=args.file
o=args.out
u=args.umi

#to run
#bartlett_deduper.py -u /Users/emilybartlett/bioinfo/Bi624/Deduper-ebart/STL96.txt -f /Users/emilybartlett/bioinfo/Bi624/Deduper-ebart/test.sam -o /Users/emilybartlett/bioinfo/Bi624/Deduper-ebart/try

#-----------------------------------------------------------------------------------------------------------------------------
#                                            Functions
#-----------------------------------------------------------------------------------------------------------------------------

def get_UMI(file):
    #```This function takes in a file that includes the known UMIs. It will then return all known UMIs as a list ```
    with open(u, "rt") as fh1:
        known_UMI = [] #create empty list to hold UMIs
        for line in fh1:
            known_UMI.append(line.rstrip('\n'))#append UMIs to the known_UMI list and remove the new line character
    return(known_UMI)

def get_flag(flag):
    '''Looks at bitwise flag to determine is strand is "-" or "+"'''
    if ((int(flag) & 16) != 16):
        strand = '+'
    else:
        strand = '-'
    return strand

def QNAME_parse(QNAME):
    #```Looks at the QNAME (column 1 on SAM) and split the qname by the ":" and return the last element```
    qname_list = QNAME.split(":") #split the first column in SAM file by ":" and store
    return qname_list[-1] #return barcode, last indexed item in list

#NOTE: re.findall regex command will iterate over a string (in our case the CIGAR string) and find characters that match a specified pattern.
#syntax: re.findall(<pattern>, string)

def  adjust_position(line):
    '''Adjust the start position if there is soft clipping for both forward and reverse strands'''
    #The input will be the whole read line, so assign variables
    CIGAR=line.split()[5]#(col6/index5)
    pos = int(line.split()[3])#(col4/index3)
    flag = line.split()[1]#(col2/index1)
    
    if ((int(flag) & 16) != 16): #Check the flag, proceed with this section if it is forward
        #find the following pattern in the CIGAR string:one or more digits(/d+) followed by M N D or S [MNS]
        cigar_1 = re.findall(r'^(\d+)([MNDS])', CIGAR)[0]
        if cigar_1[1] == 'M': #if there is a "M" in the first tuple dont adjust pos
            pos_adj = pos
        if cigar_1[1] == 'S':
            pos_adj = pos - int(cigar_1[0]) #if there is a "S" in the first tuple, adjusts its corresponding digit to an in and subract it from the original position to get adjusted position

    else: #If its not forward proceed as if the strand is reverse
        cigar_2 = re.findall(r'(\d+)([MNDS])$', CIGAR)[0]
        if cigar_2[1] == 'M':#if there is a "M" 
            pos_adj = pos + sum([int(i) for i in re.findall(r'(\d+)[MND]', CIGAR)]) - 1#sum all matching bases minus 1 to adjust position
        if cigar_2[1] == 'S':#soft clipping at the end of CIGAR string, position is adjusted by adding the sum of all bases, AND clipping, and subtracting 1
            pos_adj = pos + sum([int(i) for i in re.findall(r'(\d+)[MND]', CIGAR)]) + int(cigar_2[0]) -1
        
    return str(pos_adj)

#-----------------------------------------------------------------------------------------------------------------------------
#                                            De Duplicating Code
#-----------------------------------------------------------------------------------------------------------------------------

#call get_UMI function to return the list of known UMIs from arparsed "--umi" file
known_umis = get_UMI(u)
print(known_umis)
#Set Gloabl Variables to Tally
total_reads = 0
total_duplicates = 0
Unknown_UMIs = 0

with open(f, 'rt') as fh_i, open(o, 'wt') as fh_o:

    READID_set = set()#create empty set to compare reads to
    prev_chromosome = ''#Initiate a way to track chromosome (when you jump from one chromosome to another there is no way it is a duplicate)

    for line in fh_i:
        if line.startswith('@'):#If the line starts with an @ it is a header and can immediatly be written to our out file
            fh_o.write(line)
            continue
        total_reads += 1
        chromosome = str(line.split()[2]) #get chromosome from RNAME (col3/index2)
        if chromosome != prev_chromosome:
            # Reset READID_set
            #If our current read chromosome doesnt equal our previous chromosome, we have moved onto the next chromosome so reset the set
            READID_set = set()
        UMI = QNAME_parse(line.split()[0])#Assign reads UMI to variable
        #print(UMI)
        if UMI not in known_umis:#Make sure it is a known UMI
            Unknown_UMIs += 1
            continue        
        pos_adj = adjust_position(line)#call adjust_position function to compare start locations
        strand = get_flag(line.split()[1])
        ID = chromosome+'_'+pos_adj+'_'+strand+UMI#create unique ID for read
       
        if ID not in READID_set:
            READID_set.add(ID)
            fh_o.write(line)
        else:
            total_duplicates += 1
        
        # store rname_prev
        prev_chromosome = chromosome

#-----------------------------------------------------------------------------------------------------------------------------
#                                            Write A Summary File
#-----------------------------------------------------------------------------------------------------------------------------

#Maths
percent_dup = ((total_duplicates)/total_reads)*100

with open(f'Summary.txt', 'w') as fout:
    fout.write(f'Summary:\n\nTotal Number of Reads Parsed:\n{total_reads}\n\n')
    fout.write(f'Total Number of Reads Considered Duplicates:\n{total_duplicates}\n\n')
    fout.write(f'Total Number of Reads With Unknown UMIs:\n{Unknown_UMIs}\n\n')
    fout.write(f'Percent of Total Reads Encountered That Were Discarded Due to Unknown UMI or Suspected Duplicate:\n{percent_dup}%\n\n')
