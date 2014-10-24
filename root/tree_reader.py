import rootpy
from rootpy.io import root_open
import root_numpy
from utils.functions import identity
from utils.matrix import rec_to_mat

class tree_reader(object):
    """
    tree_reader: a class provding some basic wrapper functionality around
    record arrays, binning, and TTrees in general.

    Adopted from tree_reader class written by Luke de Oliveria
    """
    def __init__(self, infile = None, tree = None):
        super(tree_reader, self).__init__()
        if infile is not None:
            if infile.__class__ == str:
                self.files = infile
                self.f = root_open(infile)
            else:
                raise TypeError('need to specify ROOT file as a string')
        else:
            self.files = None
        self.tree_name = tree
        self.binning = {}
        self.array = None
        self.T = None

    def set_file(self, filename):
        if filename.__class__ == str:
            self.files = filename
            self.f = root_open(filename)
        else:
            raise TypeError('need to specify ROOT file as a string')
    
    def get_tree(self, treename):
        if treename is None and self.tree_name is None:
            raise NameError('No tree name found.')
        else:
            self.tree_name = treename
            self.T = self.f[treename]

    def get_array(self, branches=None, selection=None, start=None, stop=None, step=None):
        self.array = root_numpy.tree2array(self.T, branches=branches, selection=selection, start=start, stop=stop, step=step)

    def to_array(self, branches=None, selection=None, start=None, stop=None, step=None):
        if self.array is not None:
            return self.array
        return root_numpy.tree2array(self.T, branches=branches, selection=selection, start=start, stop=stop, step=step)

    def add_binning(self, rule, binning):
        if rule.__class__ is not str:
            raise TypeError('Binning rule needs to be a string.')
        if binning.__class__ is not list:
            raise TypeError('Bins need to be a numeric list')

        self.binning.update( {rule : binning} )

    def generate(self):
        self.array = (recfunctions.append_fields(
            self.array, ['categ_' + _abs_strip(name) for name in self.binning.keys()],
            [_generate_bins(data, rule, binning) for rule, binning in self.binning.iteritems()],
            usemask=False, asrecarray=True))
        for v in ('categ_' + _abs_strip(name) for name in self.binning.keys()):
            self.array = self.array[self.array[v] > -1]

    def pull_matrix(self, varlist = None):
        if varlist is not None:
            return rec_to_mat(self.array[varlist])
        return rec_to_mat(self.array)

    def _find_bin(point, strategy, abs_val = identity):
        if strategy.__class__ == list:
            nbins = len(strategy)
        else:
            nbins = strategy.shape[0]
        for i in range(1, nbins):
            if abs_val(point) >= strategy[i - 1] and abs_val(point) < strategy[i]:
                return i - 1;
        return -1

    def _abs_strip(string):
        return string.replace('abs(', '').replace(')', '')

    def _generate_bins(data, rule, binning):
        if 'abs(' in rule:
            name = _abs_strip(rule)
            return np.array([float(_find_bin(x, binning, np.abs)) for x in data[name]])
        return np.array([float(_find_bin(x, binning)) for x in data[rule]])
