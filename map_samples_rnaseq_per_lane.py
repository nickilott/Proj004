
################################################
# map wellcome trust to Sumeets identifiers in
# the 3' mRNA experiment
################################################

import re
import sys
import os

datadirs = ["/gfs/archive/powrie/proj004/rnaseq/bsg-ftp.well.ox.ac.uk/180511_K00198_0309_BHTTWFBBXX",
            "/gfs/archive/powrie/proj004/rnaseq/bsg-ftp.well.ox.ac.uk/180517_K00150_0325_AHV5N3BBXX"]

# just keep Sumeets samples (there are other data there)
sumeets = ["D" + str(i) for i in range(1,5)]

for treatment in ["CPZ", "Sal", "CPZ-Sal"]:
    sumeets = sumeets + [x + treatment for x in sumeets]

###############################
###############################
###############################

def oid2newid(wtid, oid):
    '''
    create new id per lane
    '''
    rep = oid[1]
    rep = "R" + rep
    if len(oid) == 2:
        treatment = "Control"
    elif len(oid) == 5:
        if "CPZ" in oid:
            treatment = "cpz"
        else:
            treatment = "Salmonella"
    else:
        treatment = "cpzSalmonella"

    lane = wtid.split("_")[1]
    newid = "-".join([treatment, lane, rep])
    return(newid)

###############################
###############################
###############################

def makeLinkRead1(datadirs, wtid, newid):
    '''
    create soft link from original file
    to new file name in current directory
    '''
    original_fastq = wtid + "_1.fastq.gz"
    
    # find original file
    for d in datadirs:
        tofind = os.path.join(d, original_fastq)
        if os.path.exists(tofind):
            os.symlink(tofind, newid + ".fastq.1.gz")

###############################
###############################
###############################

def makeLinkRead2(datadirs, wtid, newid):
    '''
    create soft link from original file
    to new file name in current directory
    '''
    original_fastq = wtid + "_2.fastq.gz"
    
    # find original file
    for d in datadirs:
        tofind = os.path.join(d, original_fastq)
        if os.path.exists(tofind):
            os.symlink(tofind, newid + ".fastq.2.gz")
            

###############################
###############################
###############################

# do it!

for line in open("/gfs/archive/powrie/proj004/rnaseq/identifiers/identifiers.tsv"):
    data = line[:-1].split("\t")
    wtid, oid = data
    if oid not in sumeets:
        continue
    newid = oid2newid(wtid, oid)
    makeLinkRead1(datadirs, wtid, newid)
    makeLinkRead2(datadirs, wtid, newid)



            
    
