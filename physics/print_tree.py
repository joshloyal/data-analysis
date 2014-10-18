#!/usr/bin/env python

import getopt, os, sys
from ROOT import TFile
from ROOT import gDirectory
from prettyprint import make_table

def usage():
    print ''
    print '-i | --input: input root file'
    print '-t | --tree: tree name'

try:
    shortops = 'i:t:'
    longopts = ['input=','tree=']
    opts, args = getopt.getopt(sys.argv[1:], shortops, longopts)

except getopt.GetoptError:
    print >> sys.stderr, 'ERROR: options unknown in %s' %sys.argv[1:]
    usage()
    sys.exit(1)

for o, a in opts:
    if o in('--input', '-i'):
        file_name = a
    if o in ('--tree', '-t'):
        tree_name = a

if not file_name or not tree_name:
    usage()
    sys.exit(1)

# open the file
myfile = TFile( file_name )

# retrieve the tree of interest
mytree = gDirectory.Get( tree_name )

# loop over all branches in the tree
branches = []
for branch in mytree.GetListOfBranches():

    # get the variable name and class type of each branch
    branch_name = branch.GetName()
    type_name = branch.GetClassName()
    branch_size =  branch.GetEntries()

    # if the branch is not a class (i.e. vector<object>) look at the
    # data types of the individual leaves
    if not type_name:
        for l in branch.GetListOfLeaves():
            type_name = branch.GetLeaf(l.GetName()).GetTypeName()
            break
    
    branches.append( (branch_name, type_name, branch_size) )

table = make_table(branches, header=('Variable','Type', '# Entries'), justify='L')
print table
