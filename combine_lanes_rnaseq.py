##################################################
##################################################
# A script to combine samples across three lanes
# of data
##################################################
##################################################

import itertools
import glob
import os

cgatdir="/gfs/devel/nilott/py36-v1/cgat-scripts/CGAT/scripts/"
datadir="/gfs/work/nilott/proj004/data/rnaseq2"

def getSample(filename):
    '''
    return the condition-replicate
    filename prefix from per-lane
    filename
    '''
    sample = os.path.basename(filename).split("-")
    sample = "-".join(["macrophage", sample[0], sample[2]])
    return(sample)

def combineAcrossLanes(infiles):
    '''
    combine files from the same sample
    across the three lanes
    '''
    for file1, file2, file3 in itertools.combinations(infiles, 3):
        s1 = getSample(file1)
        s2 = getSample(file2)
        s3 = getSample(file3)
        if s1 == s2 == s3:
            statement = """zcat %s %s %s | python %s/fastq2fastq.py --method=sort --log=%s.log | gzip > %s""" % (file1, file2, file3, cgatdir, s1, s1)
            os.system(statement)

combineAcrossLanes(glob.glob(os.path.join(datadir, "*.fastq.1.gz")))
combineAcrossLanes(glob.glob(os.path.join(datadir, "*.fastq.2.gz")))
