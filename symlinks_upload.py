##################################################################
##################################################################
##################################################################
# This is a script to simply take the symlinked files with
# meaningful names and get the original file names to create
# a set of data to be uploaded to arrayExpress
##################################################################
##################################################################
##################################################################

import glob
import os
import sys

DIR=os.path.abspath(sys.argv[1])

linkedfiles = glob.glob(os.path.join(DIR, "*.fastq*.gz"))

# get total number of linked files
l = len(linkedfiles)

# check its the correct number -> (nlanes*nsamples)*2 (paired end)
assert len(linkedfiles) == (3*16)*2, "Not expected number of fastq files actual count = %i" % l

# output a text file with the full path to the original data files
# as well as a mapping to the newname for reference

outf = open("original_files.txt", "w")
outf_map = open("sample_mapping.txt", "w")

for inf in linkedfiles:
    orig = os.readlink(inf)
    outf.write(orig + " ")
    outf_map.write("\t".join([os.path.basename(inf), orig]) + "\n")
outf.close()
