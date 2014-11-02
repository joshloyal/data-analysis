#!/usr/bin/env python

'''abcd data driven utils and tasks for Wg and Zg'''

__author__ = 'Mia Liu(ml149@duke.edu) '

import sys
import re
import math
import ROOT as rt
rt.gROOT.SetBatch(rt.kTRUE)
rt.gROOT.ProcessLine( '.L ./loadVector.C+')

__temp = sys.argv
sys.argv = []
# Replace sys.argv for user parsing
sys.argv = __temp
###set photon ID tool location to read the function !
#RootCoreDir="/Users/miaoyuan/Analysis/packages/RootCore/lib/"
RootCoreDir="/afs/cern.ch/work/m/mliu/private/packages/RootCore/lib/"
rt.gROOT.LoadMacro(RootCoreDir+'libPhotonIDTool.so')
rt.gROOT.LoadMacro(RootCoreDir+'libFudgeMCTool.so')

#photon tightness definition
def idLevels(phIDObj):
    ph_isEM = phIDObj.isEM(4,2012)
    isEM = not ((ph_isEM & 130050))
#    isLoose =  not ((ph_isEM & 0x04fc01))
#    isLoose = phIDObj.PhotonCutsLoose(4) 
    isLoose = False
    isTight = False
    wantBits = [17,18,19,21] #bits can be changed for systematics
    bitsum = sum(((ph_isEM >> bb) & 1) for bb in wantBits)
    if bitsum > 1:
       isLoose = True
#    isTight = isEM
    isTight = phIDObj.PhotonCutsTight(2012)    
#    isTight = not isLoose
#    isLoose = not isEM
#    isLoose = False
    return isLoose, isTight

def performID(ph_pt,shower,isconv,isMC=False):
    from array import array
    #isEMVals = set()
    eta2=rt.Double(shower.at(0)) 
#    eta2=shower[0]
    rhad1= rt.Double(shower.at(1))
    rhad=rt.Double(shower.at(2))
    e277=rt.Double(shower.at(3))
    reta=rt.Double(shower.at(4))
    rphi= rt.Double(shower.at(5))
    weta2=rt.Double(shower.at(6))
    f1=rt.Double(shower.at(7))
    fside=rt.Double(shower.at(8))
    wtot=rt.Double(shower.at(9))
    w1 = rt.Double(shower.at(10))
    deltae = rt.Double(shower.at(11))
    eratio = rt.Double(shower.at(12))
    isconv = int(shower.at(13))
    GammaSelection = rt.PhotonIDTool(ph_pt,eta2,rhad1,rhad,e277,reta,rphi,weta2,f1,fside,wtot,w1,deltae,eratio,isconv)
    isLoose, isTight = idLevels(GammaSelection)
    return isLoose, isTight

#ABCD loop

def countLoop(chain, maxEntries=None, nFB=-1,IsoVar='ph_TopoIso20',isocut=4000,nonisocut=6000,isMC=False):
    from array import array
    chain.SetBranchStatus("*", 1)
    nEntries = chain.GetEntries()
    print nEntries
    if maxEntries is not None and maxEntries >= 0:
        nEntries = min(maxEntries, nEntries)
    printPoint = max(1, nEntries/10)

    counts = {}
    elists = {}
    numbercounts = {}
    isos = {}

    _lowPtEdges  = array('d', [10,20,30,40,60,100,1000])
    _lowEtaEdges = array('d', [0,1.37,1.8,2.37])

    #name four regions and store numbers in 2D hists for pt/eta bin purpose
    for kk in ['sigReg', 'dblRev', 'nonIso', 'nonTight']:
        counts[kk] = rt.TH2F(kk,'', len(_lowPtEdges ) - 1, _lowPtEdges,len(_lowEtaEdges ) - 1,_lowEtaEdges)
        elists[kk] = rt.TEventList(kk,'',0,0)
        numbercounts[kk] = []
    #check whether Isolation has pile up dependence
    
    isEMVals = set()

    for iEntry in range(nEntries):
        if iEntry%printPoint == 0:
            print 'Working on % 6d of % 6d events' % (iEntry, nEntries)
            sys.stdout.flush()
        chain.GetEntry(iEntry)
        TightPhotonCollection = []
        LoosePhotonCollection = []
        for iPhoton in range(chain.ph_E.size()):
            eta = abs(chain.ph_Eta.at(iPhoton))
            pt  = chain.ph_Pt.at(iPhoton)/1000
            #excluding crack region?.... check
            if  1.37 < eta < 1.52           \
                or eta > _lowEtaEdges[-1]   \
                or pt < _lowPtEdges[0] or pt > _lowPtEdges[-1]:
                continue
            ph_iso = rt.Double(chain.__getattr__(IsoVar).at(iPhoton))
            ph_pt = rt.Double(chain.ph_Pt.at(iPhoton))
            shower=chain.ph_showershape.at(iPhoton)
            isconv = chain.ph_conversion[0][0]
            isLoose,isTight = performID(ph_pt,shower,isconv,isMC=False)

            if isTight: 
               TightPhoton = pt,eta,abs(ph_iso)
               TightPhotonCollection.append(TightPhoton)
            if isLoose and not isTight:
               LoosePhoton = pt,eta,abs(ph_iso)
               LoosePhotonCollection.append(LoosePhoton)
            if not isLoose:
                continue
        if not len(TightPhotonCollection) and not len(LoosePhotonCollection):
            continue

        pileupkey = ''
        if chain.averageIntPerXing <= 15:
           pileupkey = '_mu_0_15'
        else:
            if chain.averageIntPerXing < 25 and chain.averageIntPerXing > 15:
               pileupkey = '_mu_15_25'
            else:
                if chain.averageIntPerXing >= 25:
                   pileupkey = '_mu_25_50'
     
        if len(TightPhotonCollection):
            TightPhotonCollection.sort(reverse = True)
            isos['Tight_iso_total'].Fill(TightPhotonCollection[0][2])
            isos['Tight_iso'+pileupkey].Fill(TightPhotonCollection[0][2])
            isIso    = TightPhotonCollection[0][2] < isocut
            isNonIso = TightPhotonCollection[0][2] > nonisocut
            select_pt =  TightPhotonCollection[0][0]
            select_eta =  TightPhotonCollection[0][1]
            #if isIso:
            #   mllg.Fill(chain.Mllg)
        else:
            LoosePhotonCollection.sort(reverse = True)
            isos['nonTight_iso_total'].Fill(LoosePhotonCollection[0][2])
            isos['nonTight_iso'+pileupkey].Fill(LoosePhotonCollection[0][2])
            isIso    = LoosePhotonCollection[0][2] < isocut
            isNonIso = LoosePhotonCollection[0][2] > nonisocut
            select_pt =  LoosePhotonCollection[0][0]
            select_eta =  LoosePhotonCollection[0][1]
            
        key = None

        if len(TightPhotonCollection):
           if isIso:
              key = 'sigReg'
           if isNonIso:
              key = 'nonIso'
        else:
           if isIso:
              key = 'nonTight'
           if isNonIso:
              key = 'dblRev'
        if key is not None:
           weight = 1 if nFB < 0 else nFB*1000*chain.weights.at(0)
           counts[key].Fill(select_pt, select_eta, weight)
           elists[key].Enter(iEntry)

    for kk, vv in counts.iteritems():
        print kk, vv.GetEntries(), vv.GetSumOfWeights()
        for i in range(6):
            numbercounts[kk].append(vv.GetBinContent(i+1,1)+vv.GetBinContent(i+1,2)+vv.GetBinContent(i+1,3))
#    for i in range(6):    
#        print 'purity :',(numbercounts['sigReg'][i]-numbercounts['nonIso'][i]*numbercounts['nonTight'][i]/numbercounts['dblRev'][i])/numbercounts['sigReg'][i]
#        print 'R: ', (numbercounts['sigReg'][i]*numbercounts['dblRev'][i]/(numbercounts['nonIso'][i]*numbercounts['nonTight'][i]))
    n_bkg = 1/(float(elists['dblRev'].GetN())/(float(elists['nonTight'].GetN())*float(elists['nonIso'].GetN())))
    n_sigReg = float(elists['sigReg'].GetN())
    n_nonIso = float(elists['nonIso'].GetN())
    n_nonTight = float(elists['nonTight'].GetN())
    n_dblRev = float(elists['dblRev'].GetN())

    #c_nonTight = n_nonTight/n_sigReg
    #c_nonIso = n_nonIso/n_sigReg
    #c_dblRev = n_dblRev/n_sigReg

    if isMC:
       print 'calculate signal leakage from SM zgamma events'
    print 'signal leakage in region nonTight: ',n_nonTight/n_sigReg 
    print 'signal leakage in region nonIso  : ',n_nonIso/n_sigReg
    print 'signal leakage in region dblRev  : ',n_dblRev/n_sigReg

#    c_nonTight = 0.113
#    c_nonIso = 0.012
#    c_dblRev = 0.00259
    c_nonTight = 0.12954
    c_nonIso = 0.00345
    c_dblRev = 0.0007472
 
    n_bkg = (n_nonTight-c_nonTight*n_sigReg)*(n_nonIso-c_nonIso*n_sigReg)/(n_dblRev-c_dblRev*n_sigReg)
    n_nonTight_bkg = n_nonTight-c_nonTight*n_sigReg
    n_nonIso_bkg = n_nonIso-c_nonIso*n_sigReg
    n_dblRev_bkg = n_dblRev-c_dblRev*n_sigReg
    
    scale = n_bkg/(n_nonTight+n_nonIso+n_dblRev)

    print 'n_bkg = ',n_bkg
    print 'n_sig = ',n_sigReg-n_bkg
    return counts,isos,elists,scale

def count2dABCD(chain, maxEntries=None,scale=1,IsoVar='ph_TopoIso20',isocut=3000,nonisocut=6000,isMC=False):
    from array import array

    chain.SetBranchStatus("*", 1)
    nEntries = chain.GetEntries()
    print nEntries
    if maxEntries is not None and maxEntries >= 0:
        nEntries = min(maxEntries, nEntries)
    printPoint = max(1, nEntries/10)

    counts = {}
    counts1 = {}
    elists = {}
    elists1 = {}
    numbercounts = {}
    numbercounts1 = {}
    isos = {}

    _lowPtEdges  = array('d', [15,20,30,40,60,100,10000])
    #_lowPtEdges  = array('d', [40,45,55,60,65,100,1000])
    _lowEtaEdges = array('d', [0,1.37,1.8,3.4])

    #name four regions and store numbers in 2D hists for pt/eta bin purpose
    for kk in ['obs','sig','nonTight','nonIso','dblRev']:        
        counts[kk] = rt.TH2F(kk,'', len(_lowPtEdges ) - 1, _lowPtEdges,len(_lowEtaEdges ) - 1,_lowEtaEdges)
        numbercounts[kk] = []
        elists[kk] = rt.TEventList(kk,'',0,0)
    for kk in ['ssig','snonTight','snonIso','sdblRev','subsig','subnonTight','subnonIso','subdblRev']:        
        counts1[kk] = rt.TH2F(kk,'', len(_lowPtEdges ) - 1, _lowPtEdges,len(_lowEtaEdges ) - 1,_lowEtaEdges)
        numbercounts1[kk] = []
        elists1[kk] = rt.TEventList(kk,'',0,0)
    
    isEMVals = set()

    for iEntry in range(nEntries):
        if iEntry%printPoint == 0:
            print 'Working on % 6d of % 6d events' % (iEntry, nEntries)
            sys.stdout.flush()
        chain.GetEntry(iEntry)
        TightPhotonCollection = []
        PhotonCollection = []
        nonTightPhotonCollection = []
        LoosePhotonCollection = []
       # if chain.dRgg < 0.7 or chain.dRlg[0] < 0.7 or chain.dRlg[1]< 0.7:
       #    continue
        for iPhoton in range(chain.ph_E.size()):
            eta = abs(chain.ph_Eta.at(iPhoton))
            pt  = chain.ph_Pt.at(iPhoton)/1000
            #excluding crack region?.... yes, always do this for photons
            ph_iso = rt.Double(chain.__getattr__(IsoVar).at(iPhoton))
            ph_pt = chain.ph_Pt.at(iPhoton)
            ph_phi = chain.ph_Phi.at(iPhoton)
            ph_E = chain.ph_E.at(iPhoton)
            #ph_pt = rt.Double(chain.ph_Pt.at(iPhoton))
            #ph_phi = rt.Double(chain.ph_Phi.at(iPhoton))
            #ph_E = rt.Double(chain.ph_E.at(iPhoton))
            shower=chain.ph_showershape.at(iPhoton)
            isconv = chain.ph_conversion[iPhoton][0]
            isLoose,isTight = performID(ph_pt,shower,isconv,isMC=False)
            #if isLoose and not isTight:
#            isTight = True
            phlv = rt.TLorentzVector()
            phlv.SetPtEtaPhiE(ph_pt,chain.ph_Eta.at(iPhoton),ph_phi,ph_E)
            if isLoose and not isTight:
               LoosePhoton = pt,eta,abs(ph_iso),isTight,phlv
               LoosePhotonCollection.append(LoosePhoton)
            if isTight:
               TightPhoton = pt,eta,abs(ph_iso),isTight,phlv
               TightPhotonCollection.append(TightPhoton)
        #if len(TightPhotonCollection)!=2:
           #print 'tight ph',len(TightPhotonCollection)
        #print 'loose ph',len(LoosePhotonCollection)
        if len(TightPhotonCollection) > 1:
           TightPhotonCollection.sort(reverse = True)
           PhotonCollection.append(TightPhotonCollection[0])
           PhotonCollection.append(TightPhotonCollection[1])
        if len(TightPhotonCollection) == 1 and len(LoosePhotonCollection) >= 1:
           TightPhotonCollection.sort(reverse = True)
           LoosePhotonCollection.sort(reverse = True)
           if TightPhotonCollection[0][0]>LoosePhotonCollection[0][0]:
              PhotonCollection.append(TightPhotonCollection[0])
              PhotonCollection.append(LoosePhotonCollection[0])
           if TightPhotonCollection[0][0]<LoosePhotonCollection[0][0]:
              PhotonCollection.append(LoosePhotonCollection[0])
              PhotonCollection.append(TightPhotonCollection[0])
        if len(LoosePhotonCollection) > 1 and len(TightPhotonCollection)==0:
           LoosePhotonCollection.sort(reverse = True)
           PhotonCollection.append(LoosePhotonCollection[0])
           PhotonCollection.append(LoosePhotonCollection[1])
           
 
        if len(PhotonCollection)!=2:
           continue

        PhotonCollection.sort(reverse = True)
        isIso1    = PhotonCollection[0][2] <= isocut
        isNonIso1 = PhotonCollection[0][2] >= nonisocut
        isTight1  = PhotonCollection[0][3]  
        isIso2    = PhotonCollection[1][2] <= isocut
        isNonIso2 = PhotonCollection[1][2] >= nonisocut
        isTight2  = PhotonCollection[1][3]  
        select_pt1  = PhotonCollection[0][0]  
        select_eta1  = PhotonCollection[0][1]  
        select_pt2  = PhotonCollection[1][0]  
        select_eta2  = PhotonCollection[1][1]  

        elv = rt.TLorentzVector()
        elv.SetPtEtaPhiE(chain.tightlep_Pt[0],chain.tightlep_Eta[0],chain.tightlep_Phi[0],chain.tightLep_E[0])
        mlgg = rt.Double((elv+PhotonCollection[0][4]+PhotonCollection[1][4]).M())
        Ptlgg = (elv+PhotonCollection[0][4]+PhotonCollection[1][4]).Pt()
        mlg1 = (elv+PhotonCollection[0][4]).M()
        mlg2 = (elv+PhotonCollection[1][4]).M()
        dRgg = PhotonCollection[0][4].DeltaR(PhotonCollection[1][4])
        if mlgg > 80000 and mlgg < 100000 or Ptlgg < 20000 or mlg1 > 80000 and mlg1<100000 or mlg2 > 80000 and mlg2<100000 or dRgg < 0.7:
        #if mlgg > 80 and mlgg < 1000:
        #if chain.tightlep_Pt[0] < 25000 or mlgg > 80000 and mlgg < 100000:
           continue

        key = None
        key2 = None
        key3 = None

        if len(PhotonCollection) == 2:
           if mlgg>80000:
              if isIso1 and isTight1:
                 key = 'sig'
              if isIso1 and not isTight1:
                 key = 'nonTight'
              if isNonIso1 and isTight1:
                 key = 'nonIso'
              if isNonIso1 and not isTight1:
                 key = 'dblRev'
           if isIso2 and isTight2:
              key2 = 'subsig'
           if isIso2 and not isTight2:
              key2 = 'subnonTight'
           if isNonIso2 and isTight2:
              key2 = 'subnonIso'
           if isNonIso2 and not isTight2:
              key2 = 'subdblRev'
           if isIso1 and isTight1:
              if isIso2 and isTight2:
                 key3='ssig'
              if isNonIso2 and isTight2:
                 key3='snonIso'
              if isIso2 and not isTight2:
                 key3='snonTight'
              if isNonIso2 and not isTight2:
                 key3='sdblRev'
        
        weight = chain.weights.at(1)*scale
#        weight = scale
        if isIso1 and isTight1 and isIso2 and isTight2:
           counts['obs'].Fill(select_pt1, select_eta1, weight)
        if key is not None:
           counts[key].Fill(select_pt1, select_eta1, weight)
           elists[key].Enter(iEntry)
        if key2 is not None:
           counts1[key2].Fill(select_pt2, select_eta2, weight)
           elists1[key2].Enter(iEntry)
        if key3 is not None:
           counts1[key3].Fill(select_pt2, select_eta2, weight)
           elists1[key3].Enter(iEntry)

    for kk, vv in counts.iteritems():
#        print kk, vv.GetEntries(), vv.GetSumOfWeights()
        number = float(0)
        for i in range(6):
            number= number + vv.GetBinContent(i+1,1)+vv.GetBinContent(i+1,2)+vv.GetBinContent(i+1,3)
            #print kk,'bin1,2,3',vv.GetBinContent(i+1,1),vv.GetBinContent(i+1,2),vv.GetBinContent(i+1,3)
            #number = number + vv.GetSumOfWeights(i+1,1)+vv.GetSumOfWeights(i+1,2)+vv.GetSumOfWeights(i+1,3)
        #numbercounts[kk].append(number)
        numbercounts[kk]=number
    for kk, vv in counts1.iteritems():
#        print kk, vv.GetEntries(), vv.GetSumOfWeights()
        number = float(0)
        for i in range(6):
            number = number+vv.GetBinContent(i+1,1)+vv.GetBinContent(i+1,2)+vv.GetBinContent(i+1,3)
        numbercounts1[kk]=number

    return numbercounts,numbercounts1,elists,elists1

if __name__ == '__main__':
#LOADER probably still needed.  I think its on smuhpc:~/Loader.C
    dataCh = rt.TChain("WggCandidates")
    datadir ='../SFrameOutputs/'
    #dataCh.Add(datadir+'Zeegamma2011EventSelection.data.Zg.B.v0.root')
    #print type(dataCh)
    dataCh.Add(datadir+'WenugammaEventSelection.data.v0.root')
    #dataCh.Add(datadir+'Zeegamma2011EventSelection.MC.Zjets.v0.root')
   # dataCounts,isos,elists,scale = countLoop(dataCh, maxEntries=-1, nFB=-1)
#    dataCounts,isos,elists,scale = countLoop(dataCh, maxEntries=None, nFB=-1,IsoVar = 'ph_TopoIso20',isocut=4000,nonisocut=6000,isMC=False)
    dataCounts,dataCounts1,elists,elists2 = count2dABCD(dataCh, maxEntries=None,scale=1,IsoVar='ph_TopoIso20',isocut=3000,nonisocut=6000,isMC=False)
    print dataCounts,dataCounts1 
    #dataCounts1,isos1,elists1,scale1 = countLoop(dataCh, maxEntries=None, nFB=-1,IsoVar = 'ph_TopoIso30',isocut=4000,nonisocut=5000,isMC=False)
    #dataCounts2,isos2,elists2,scale2 = countLoop(dataCh, maxEntries=None, nFB=1,IsoVar = 'ph_TopoIso20',isocut=3000,nonisocut=4000,isMC=False)
