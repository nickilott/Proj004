
import sys
import collections


sample_map = collections.OrderedDict([("Donor11", "UntreatedUnInfected"),
                                      ("Donor12", "UntreatedInfected"),
                                      ("Donor13", "CPZUninfected"),
                                      ("Donor14", "CPZInfected"),
                                      ("Donor21", "UntreatedUnInfected"),
                                      ("Donor22", "UntreatedInfected"),
                                      ("Donor23", "CPZUninfected"),
                                      ("Donor24", "CPZInfected"),
                                      ("Donor31", "UntreatedUnInfected"),
                                      ("Donor32", "UntreatedInfected"),
                                      ("Donor33", "CPZUninfected"),
                                      ("Donor34", "CPZInfected")])

for id, condition in sample_map.iteritems():
    sys.stdout.write("%s\t%s\n" % (id, condition))
