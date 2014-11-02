#!/usr/bin/env python
import sys
import numpy as np

if __name__ == '__main__':
	if len(sys.argv) < 3:
		print 'Usage: {} INFILE TREENAME OUTFILE'.format(sys.argv[0])
		sys.exit(0)
	else:
		import rootpy
		from rootpy.tree import Tree, TreeModel, TreeChain, Cut
		from rootpy.io import root_open
		rootpy.log.basic_config_colorized()
		from root_numpy import tree2rec, array2tree

		f = root_open(sys.argv[1])
		treename = (sys.argv[2])
		f_out = root_open(sys.argv[3], 'recreate')

	T = f[treename]

	print 'Generating...' 

	evts = T.to_array()

	mid = evts.shape[0] / 2

	idx = range(0, evts.shape[0])
	np.random.shuffle(idx)
	evts = evts[idx]

	train_evts = evts[:mid]
	test_evts = evts[mid:]

	print 'Writing...'
	T_train = array2tree(train_evts, 'train_ntup')
	T_test = array2tree(test_evts, 'test_ntup')

	T_train.Write()
	T_test.Write()

	f_out.Close()