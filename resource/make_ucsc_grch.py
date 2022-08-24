#! /usr/bin/env python

from __future__ import print_function
import sys

input_file = sys.argv[1]

with open(input_file, 'r') as hin:
    for line in hin:
        if line.startswith('#'): continue
        F = line.rstrip('\n\r').split('\t')
        if F[4].startswith('CM'):
            print(F[2] + '\t' + F[9])
        else:
            print(F[4] + '\t' + F[9])

