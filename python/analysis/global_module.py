# global constants and functions for Compressed SUSY Analysis
# J.Loyal (Yale)

import sys
sys.path.append("modules/")
# import atlas styles
from AtlasStyle import *
from AtlasUtils import *
from ROOT import *
from array import *
from math import *
import ROOT
from ROOT import TFile
import math


################################################################
print  "Import global functions from global_module.py. S.Chekanov"
size_symbol=1.2
size_label=0.035
size_title=0.045
size_line=3

dataLumi=5835.21
fdataLumi=dataLumi/1000
fdataLumi=" %.1f fb^{-1}" % fdataLumi 

figdir="figs/"

TotLumi=str(dataLumi)
DataLab="Data 2012 #sqrt{s}=8 TeV "+fdataLumi
KinemCuts="p_{T}>100 GeV |#eta|<1.37";
KinemCutsF="p_{T}>100 GeV 1.52<|#eta|<2.37";
ATLASlab="ATLAS Work In Progress"
ElectronKinemCuts="p^{e}_{T}>25 GeV |#eta|<2.47"
MuonKinemCuts="p^{#mu}_{T}>25 |#eta|<2.4"
PhotonKinemCuts="p^{#gamma}_{T}>15 GeV 1.37<|#eta| and 1.52<|#eta|<2.37"

'''
# files with systematics 
data=[]
dp70=[]
dp140=[]
dp280=[]
dp500=[]

jf70=[]
jf140=[]
jf240=[]
jf500=[]

# generate file list 
for i in range(0,Nsys):
  sdir='../out/sys'+str(i)+"/"
  if (i != 5 and i != 6 ):
    data.append( TFile(sdir+"data12.root") )
    dp70.append( TFile(sdir+"dp70.root"))
    dp140.append( TFile(sdir+"dp140.root"))
    dp280.append( TFile(sdir+"dp280.root"))
    dp500.append( TFile(sdir+"dp500.root"))
    jf70.append( TFile(sdir+"jf70.root"))
    jf140.append( TFile(sdir+"jf140.root"))
    jf240.append( TFile(sdir+"jf240.root"))
    jf500.append( TFile(sdir+"jf500.root"))

  print "Read=",i, sdir
'''

def getHistoStyle(h1,style=1,col=21,xlab="none",ylab="none"):
           global size_symbol,size_label,size_title,size_line
           h1.Sumw2()
           h1.SetTitle("");
           h1.SetFillColor(col)
           h1.SetLineStyle(style)
           h1.SetLineWidth(size_line)
           h1.SetMarkerStyle(20)
           h1.SetMarkerSize(size_symbol)
           h1.SetStats(0)
           # if (h1.Integral()>0):  h1.Scale(1/h1.Integral());
           h1.SetAxisRange(0,0.15,"y");
           ax=h1.GetXaxis(); ax.SetTitleOffset(0.8)
           ax.SetTitle( xlab );
           ay=h1.GetYaxis(); ay.SetTitleOffset(1.2)
           ay.SetTitle( ylab );
           ax.SetTitleColor(1); ay.SetTitleColor(1);
           ax.SetTitleOffset(1.1); ay.SetTitleOffset(1.6)
           ax.SetTitleSize(size_title); ay.SetTitleSize(size_title);
           ax.SetLabelSize(size_label); ay.SetLabelSize(size_label);
           ax.SetLabelOffset(.015); ay.SetLabelOffset(.015);
           return h1,ax,ay


def drawXAxis(sf,gPad,XMIN,YMIN,XMAX,YMAX,nameX,nameY):
 h=gPad.DrawFrame(XMIN,YMIN,XMAX,YMAX);
 ay=h.GetYaxis();
 if (sf==1): ay.SetLabelSize(0.05)
 if (sf==2): ay.SetLabelSize(0.10)
# ay.SetTitleSize(0.1)
 ay.SetNdivisions(505);
 if (sf==1): 
             ay.SetTitle( nameY )
             ay.SetTitleSize(0.05)
             ay.SetTitleOffset(1.5)

 # ay.Draw("same")
 ax=h.GetXaxis(); ax.SetTitle( nameX );
 ax.SetTitleOffset(1.0)
 # if (sf==1): ax.SetTitleOffset(0.3)

 ay.SetTitleOffset(1.2)
 ax.SetLabelFont(42)
 ax.SetTitleFont(42)
 ay.SetLabelFont(42)
 ay.SetTitleFont(42)
 ax.SetLabelSize(0.10)
 ax.SetTitleSize(0.13)

 ax.Draw("same");
