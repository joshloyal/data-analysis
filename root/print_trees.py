#!/usr/bin/env python

import argparse
import rootpy
from root_numpy import list_trees, list_branches, list_structures
from rootpy.io import root_open

parser = argparse.ArgumentParser(description=
            'Print the trees and their branches which are contained in a rootfile')
parser.add_argument('filename', help='name of the root file')
args = parser.parse_args()

print args.filename


f = root_open(args.filename)
trees = list_trees(args.filename)
ntrees = len(trees)

for i, tree in enumerate(trees):
    T = f[tree]
    print ' |'
    print ' +-- %s, entries = %i' %(tree, T.GetEntries()) 
    if i == (ntrees - 1):
        print '    |'
    else:
        print ' |  |'
    branches = list_structures(args.filename, treename=tree)
    for key, item in branches.items():
        branchname, btype = item[0][0], item[0][1]
        if i == (ntrees - 1):
            print '    \-- %s : %s'%(branchname, btype)
        else:
            print ' |  \-- %s : %s'%(branchname, btype)
