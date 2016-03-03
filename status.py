#!/usr/bin/env python
from __future__ import print_function
import argparse
import glob


iloc2prot = dict()
for infile in glob.glob('species/*/*.protein2ilocus.txt'):
    with open(infile, 'r') as instream:
        for line in instream:
            protid, locusid = line.strip().split('\t')
            iloc2prot[locusid] = protid

chlorophytes = ['Apro', 'Crei', 'Cvar', 'Csub', 'Mpus', 'Msrc', 'Oluc', 'Otau', 'Vcar']

class HiLocus(object):

    def __init__(self, record):
        self._rawdata = record
        values = record.strip().split('\t')
        self.iloci = values[2].split(',')
        self.species = values[3].split(',')

    def status(self):
        assert len(self.species) >= 1 and len(self.species) <= 13
        chlorophyte_count = 0
        nonchlorophyte_count = 0
        for species in self.species:
            if species in chlorophytes:
                chlorophyte_count += 1
            else:
                nonchlorophyte_count += 1
        assert chlorophyte_count <= 9 and nonchlorophyte_count <= 4 and chlorophyte_count + nonchlorophyte_count > 0

        if chlorophyte_count == 9:
            return 'HighlyConserved'
        elif chlorophyte_count > 3:
            return 'Conserved'
        elif chlorophyte_count + nonchlorophyte_count > 1:
            return 'Matched'
        else:
            return 'Orphan'


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('hiloci', type=argparse.FileType('r'),
                        help='hilocus data table')
    return parser


def main(args):
    hilocus_count = 0
    print('hiLocus', 'iLocus', 'Status', 'Species', 'Protein', sep='\t')
    for line in args.hiloci:
        hilocus = HiLocus(line)
        hilocus_count += 1
        hid = 'HILC-{}'.format(hilocus_count)
        for locusid in hilocus.iloci:
            print(hid, locusid, hilocus.status(), locusid[0:4],
                  iloc2prot[locusid], sep='\t')

if __name__ == '__main__':
    main(get_parser().parse_args())
