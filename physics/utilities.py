from os import listdir
from os.path import isfile, join
from ROOT import TLegend

def get_histogram_names(tfile):
    names = []
    keys = tfile.GetListOfKeys()
    for key in keys:
        object = key.ReadObj()
        if object.InheritsFrom('TH1'):
            names.append( object.GetName() )

    return names

# list the contents of a directory
def ls_directory(path):
    files = [ f for f in listdir(path) if isfile(join(path,f)) ]
    return files

def make_legend(x1, y1, x2, y2):
    leg = TLegend(x1, y1, x2, y2)
    leg.SetBorderSize(0)
    leg.SetTextFont(62)
    leg.SetFillColor(10)
    leg.SetTextSize(0.03)

    return leg

def normalize_histogram(histo):
    integral = histo.Integral()
    if integral != 0:
        histo.Scale(1./integral)
