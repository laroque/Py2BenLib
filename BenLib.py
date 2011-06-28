'''
 This is a library of short general use functions.
 It includes lots of generic python 2-15 liners that I use often.
 Also included are several pieces of frequently used data handling and analysis.
'''

import numpy
import datetime
from math import *
import ROOT
import sys

############################################################
############################################################
def str2table(stringname):
    '''
      takes a string possibly from open(filename).readlines()
      outputs an array of arrays formed by the lines of the string
    '''
    filteredlines=[map(float,line.split()) for line in stringname if not 
                             line.startswith('##') and not line == '\n']
    table=numpy.array(filteredlines,dtype=numpy.float64)
    return table

############################################################
############################################################
def compact(x,y,dummy):
    '''
      Removes entries from arrays x and y at the same index when y[i]==dummy
    '''
    outx=[]
    outy=[]
    for x1,y1 in zip(x,y):
        if y1!=dummy:
          outx.append(x1)
          outy.append(y1)
    return numpy.array(outx,dtype=numpy.float64),numpy.array(outy,dtype=numpy.float64)

############################################################
############################################################
def table2dic(elements,table):
    '''
      Takes an input 2D array, table, and column lables, elements:
      Output is a dictionary with each key being the name of a column 
                 its value is an array of the contents
    '''
    cols={}
    for colname in elements:
        cols[colname]=[]
    for row in table:
        for i,name in enumerate(elements):
            cols[name].append(row[i])
    for key,value in cols.items():
        cols[key]=numpy.array(value,dtype=numpy.float64)
    return cols

############################################################
############################################################
def eff_calc(NetArea,Uncertainty,Activity,Branch,LiveTime):
    '''
      Net Area and Uncertainty are for the peak in question in counts
      Activity is in Becquerels
      Branch is the branching ratio
      LiveTime is in seconds

      returns abs_ef,abs_ef_unc,rel_ef,rel_ef_unc
     '''
    eff_NaI = 1.2E-3 #eff of IEEE standard NaI detector geometry
    abs_ef = NetArea/(Activity*Branch*LiveTime)
    abs_ef_unc = Uncertainty/(Activity*Branch*LiveTime)
    rel_ef = 100.*abs_ef/eff_NaI
    rel_ef_unc = 100.*abs_ef_unc/eff_NaI
    return abs_ef,abs_ef_unc,rel_ef,rel_ef_unc

############################################################
############################################################
def new_activity(RunDate,CalibrationDate,HalfLife,CalibrationActivity):
    '''
      RunDate is the date the sourse was used
      CalibrationDate is the yyyymmdd the quoted activity was measured
      HalfLife is in days
      
      output is the activity on the use date in the same units as quoted activity
    '''
    runDtInt = int(RunDate)
    runyear = runDtInt/10000
    runmonth = (runDtInt/100)%100
    runday = runDtInt%100
    runDt = datetime.date(year=runyear,month=runmonth,day=runday)
    calDtInt = int(CalibrationDate)
    runyear = calDtInt/10000
    runmonth = (calDtInt/100)%100
    runday = calDtInt%100
    calDt = datetime.date(year=runyear,month=runmonth,day=runday)
    time=runDt-calDt
    activity=float(CalibrationActivity)*e**(-log(2.)*(float(time.days)/float(HalfLife)))
    return activity

############################################################
############################################################
#def res_plot(shapetimes,FWHMs,FWHMS_ers,filename):
#    '''
#      Creaets a plot of FWHM^2 vs shaping time and fits to the
#        standard a*x+b+c/x to give components of noise
#    '''
#    elements = ['Shaping','FWHM','FWHM_unc']
#    values = table2dic(elements,zip(shapetimes,FWHMs,FWHMS_ers))
#    xerrs = numpy.zeros(len(values['Shaping']),dtype=numpy.float64)
#    yerrs = values['FWHM_unc']
#    yerrs = [yerr**2 for yerr in yerrs]
#    noisecurve = ROOT.TGraphErrors(len(values['Shaping']),values['Shaping'],
#                                   numpy.array(values['FWHM'],'d')**2,xerrs,
#                                   numpy.array(yerrs,'d'))
#    hist = ROOT.TH2D("hist","FWHM^2 vs Shaping Time",1,0.9*min(values['Shaping']),
#                     1.1*max(values['Shaping']),1,0.9*min(numpy.array(values['FWHM'],
#                     'd')**2),1.1*max(numpy.array(values['FWHM'],'d')**2))
#    hist.SetXTitle("Shaping Time [\mus]")
#    hist.SetYTitle("FWHM^2 [kev^2]")
#
#    # make a fit
#    noisefit = ROOT.TF1('noisefit','[0]*x+[1]+[2]*(1/x)',0.2,11.)
#    noisefit.SetParameters(25,150,400)
#    noisefit.SetParLimits(0,0,1000)
#    noisefit.SetParLimits(1,0,1000)
#    noisefit.SetParLimits(2,0,1000)
#    noisecurve.Fit('noisefit','RMBQ')
#    fitparams = []
#    fitparerr = []
#    for pnum in range(0,3):
#        fitparams.append(noisefit.GetParameter(pnum))
#        fitparerr.append(noisefit.GetParError(pnum))
#    paramstr = 'FWHM(#tau)=(%1.4e)#tau+(%1.3e)+(%1.4e)#frac{1}{#tau}' % (fitparams[0],
#                                                            fitparams[1],fitparams[2])
#    pt1 = ROOT.TLatex(.3,.8,paramstr)
#    pt1.SetTextSize(0.04)
#    pt1.SetNDC()
#
#    # Set styles and draw everything
#    c1 = ROOT.TCanvas('c1','',800,600)
#    c1.SetLogx(1)
#    c1.SetLogy(1)
#    hist.SetTitleSize(5,'')
#    hist.SetTitleSize(0.45,'yx')
#    hist.SetLabelSize(0.04,'yx')
#    hist.SetTitleOffset(1.4,'x')
#    hist.SetTitleOffset(1.4,'y')
#    hist.GetYaxis().SetMoreLogLabels()
#    hist.GetXaxis().SetMoreLogLabels()
#    hist.Draw()
#    noisecurve.SetMarkerStyle(3)
#    noisecurve.Draw('P')
#    noisefit.SetLineColor(2)
#    noisefit.Draw('same')
#    pt1.Draw()
#
#    # Write Output
#    plotype='.pdf'
#    c1.Print(filename.replace('.data','_res'+plotype))
#    outfile = open(filename.replace('.data','.out'),'a')
#    paramsers = '(%1.4e +/- %1.4f)tau+(%1.4e +/- %1.4e)+(%1.4e +/- %1.4e)(1/tau)' % (fitparams[0],
#                                 fitparerr[0],fitparams[1],fitparerr[1],fitparams[2],fitparerr[2])
#    print >>outfile, '\n'
#    print >>outfile, "noise curve fit parameters"
#    print >>outfile, paramsers
#    outfile.close()

#########################################
#########################################
def CreateCalibratedHist(unCalHist,slope,intercept,outputtype='f'):
    '''
      Takes an input TH1 <unCalHist> and linear calibration [Energy = ADCchannel*<slope> + <intercept>]
      Returns a TH1F or TH1D as flagged by <outputtype> of the same data with the x-axis in calibrated energy
      Binning is done to match bin numbers between input and output hists
      (No promises are made about rebinning artifacts)
    '''
    binVals=[]
    [binVals.append(unCalHist.GetBinContent(bin)) for bin in range(0,unCalHist.GetNbinsX())]
    fltset=set(['f','F'])
    dblset=set(['d','F'])
    if fltset.issuperset(outputtype): CalHist=ROOT.TH1F('Calibrated_'+unCalHist.GetName(),
                                  unCalHist.GetTitle(),len(binVals)+1,intercept,len(binVals)*slope+intercept)
    if dblset.issuperset(outputtype): CalHist=ROOT.TH1D('Calibrated_'+unCalHist.GetName(),
                                  unCalHist.GetTitle(),len(binVals)+1,intercept,len(binVals)*slope+intercept)
    if not (fltset.union(dblset)).issuperset(outputtype):
        print('outputtype not a flag for Float or Double, options are F,f,D,d')
        sys.exit(1)
    binCenters=[unCalHist.GetBinCenter(val) for val in range(0,len(binVals))]
    CalHist.FillN(len(binVals),slope*numpy.array(binCenters,'double')+intercept,numpy.array(binVals,'double'))
    return CalHist

#########################################
#########################################
def FitAPeak(fithist,lowbin,highbin):
    '''
      Takes a hist and upper and lower bins between which to do a fit.
      Returns centroid, sigma, peak area, continuum area
    '''
    ROOT.TROOT.gApplication.ExecuteFile("~/bin/MGTFitter/LoadMGTFitterClasses.C")
    ROOT.gROOT.SetBatch(0)
    from ROOT import MGTPeakFit
    pkfit = MGTPeakFit(fithist)
    pkfit.SetPause(0)
    pkfit.SetPrintFit(0)
    pkfit.SetFitVerbosity(1)
    pkfit.SetUseTail(0)
    pkfit.SetUseStep(0)
    pkfit.LinearPlusGauss(numpy.float(lowbin),numpy.float(highbin))
    return pkfit.GetFitMean(),pkfit.GetFitSigma(),pkfit.GetPeakCounts(),pkfit.GetContinuumCounts()

#########################################
#########################################
def IntegrateAPeak(inthist,lowbin,highbin,nsigma=5):
    '''
      Takes a hist and closed bin interval in which to integrate a peak.
      Finds the interval's absolute max and literal FWHM about that max.
      Integrates +/- 5 sigma about the max based on the determined FWHM.
      Background is taken as an additional +/- 5 sigma both up and down.
      Different +/- sigma may be defined by the optional sigma argument.
      If +/- 10 sigma extends past the interval then just returns False.
      
      returns area,uncertainty,sigma[bins]
    '''
    values=[]
    for i in range(0,inthist.GetNbinsX()): values.append(inthist.GetBinContent(i))
    window=values[lowbin:highbin+1]
    fwhmtest=[]
    for (i,xi) in enumerate(window):
        if xi<=(max(window)/2.): fwhmtest.append(1)
        else: fwhmtest.append(0)
    fwhmtop = fwhmtest[window.index(max(window)):].index(1)+window.index(max(window))+lowbin
    fwhmtest.reverse()
    fwhmbot = window.index(max(window))-fwhmtest[-window.index(max(window))-1:].index(1)+lowbin
    center = int(fwhmbot + (fwhmtop - fwhmbot)/2.)
    sigma = int(ceil((fwhmtop-fwhmbot)/2.35))
    if (center-2*nsigma*sigma)<lowbin or (center+2*nsigma*sigma)>highbin: (area,uncty)=(False,False)
    else:
        area = (2*inthist.Integral(center-nsigma*sigma,center+nsigma*sigma)-
               inthist.Integral(center-2*nsigma*sigma,center+2*nsigma*sigma))
        uncty = sqrt(4*inthist.Integral(center-nsigma*sigma,center+nsigma*sigma)+
                    inthist.Integral(center-2*nsigma*sigma,center+2*nsigma*sigma))
    return area,uncty,sigma
        

#########################################
#########################################
def GraphicHistCal(unCalHist):
    '''
      Take input TH1 and allow user to fit an arbitrary number of peaks and provide an energy for each.
      Store Energy,Bin pairs and create a linear calibration E=f(bin).
      Pass everything to CreateCalibratedHist.
      Return the calibrated hist and the calibration used.
    '''
    done_all_peaks = 'N'
    energyvals=[]
    adcvals=[]
    while set(done_all_peaks).issubset(['N','n','O','o']):
        (Mean,Sigma,Area,Background)=GraphicPeakFit(unCalHist)
        Energy = raw_input('What is the energy (in keV) of that peak? \n Enter 0 to discard this fit.\n')
        if Energy:
            energyvals.append(Energy)
            adcvals.append(Mean)
        done_all_peaks = raw_input('Have all desired peaks been fit?')
    energyvals=numpy.array(energyvals,'double')
    adcvals=numpy.array(adcvals,'double')
    slope=(len(energyvals)*sum(energyvals*adcvals)-sum(energyvals)*sum(adcvals))/(len(energyvals)*
                                                            sum(adcvals**2)-sum(adcvals)**2)
    intercept=(sum(adcvals**2)*sum(energyvals)-sum(adcvals)*sum(energyvals*adcvals))/(len(energyvals)*
                                                                   sum(adcvals**2)-sum(adcvals)**2)
    CalHist=CreateCalibratedHist(unCalHist,slope,intercept)
    return CalHist,slope,intercept

#########################################
#########################################
def GraphicPeakFit(HistIn):
    '''
      Take an input hist and draw it.
      Let user zoom in on an window and fit for one peak in that window.
      Allow user to approve the fit or do a refit.
      Return all fit parameters (centroid, sigma, peak_counts, continuum_counts).
    '''
    ROOT.TROOT.gApplication.ExecuteFile("~/bin/MGTFitter/LoadMGTFitterClasses.C")
    ROOT.gROOT.SetBatch(0)
    from ROOT import MGTPeakFit
    fit_is_good = 'N'
    while set(fit_is_good).issubset(['N','n','O','o']):
        HistIn.Draw()
        print('You must enter a region to fit. You may zoom on the histogram to find it.\n')
        xmin = raw_input('What is the ROI lower bound?\n')
        xmax = raw_input('What is the ROI upper bound?\n')
        (mean,sigma,peakcounts,continuumcounts)=FitAPeak(HistIn,xmin,xmax)
        fit_is_good = raw_input('Is this fit good?\n')
    return mean,sigma,peakcounts,continuumcounts

#########################################
#########################################
def IndexFromSet(listname,setname):
    '''
      Checks to see if exactly one element of <setname> exists in <listname>.
      If so, returns the result of listname.index() on that element.
      Is good for matching a pattern with finite permutations.
      Does a sys.exit(1) if there is not one and only one element.
      Perhaps this should be an exception rather than and exit?
    '''
    if len(setname.intersection(listname)) is 1: indexint= listname.index(setname.intersection(
                                                                               listname).pop())
    else:
        print('failed to find and element of: '+str(setname))
        sys.exit(1)
    return indexint

#########################################
#########################################
def GetSpeLive(filename):
    '''
      Takes a Maestro or ORCA produced .Spe and processes meta data

      return LiveTime
    '''
    lines = open(filename).readlines()
    timeset = set(['$MEAS_TIM:\n','$MEAS_TIM:\r\n'])
    timeindex = IndexFromSet(lines,timeset)+1
    LiveTime = float(lines[timeindex].split()[0])
    return LiveTime

#########################################
#########################################
def GetSpeHist(filename,rate=False):
    '''
      Take a Maestro or ORCA produced .Spe file and extract the histogram data.
      Requires $DATA or $Data before the start of the data.
      Requires a line starting with '$' after the data ends.
      Outputs a TH1F of the data.
      The TH1F has title set to the comments for the Hist.
      The TH1F y-axis is total events by default.
      If rate==True, uses the livetime to convert to event rate.

      return is a ROOT TH1
    '''
    lines = open(filename).readlines()
    data_start = False
    data_stop = False
    spec_name = False
    dataset = set(['$DATA:\n','$DATA:\r\n','$Data:\n','$Data\r\n'])
    roiflag = set(['$ROI:\n','$ROI:\r\n'])
    nameset = set(['$SPEC_ID:\n','$SPEC_ID:\r\n'])
    timeset = set(['$MEAS_TIM:\n','$MEAS_TIM:\r\n'])
    data_start = IndexFromSet(lines,dataset)+1
    data_stop = IndexFromSet(lines,roiflag)
    spec_name = lines[IndexFromSet(lines,nameset)+1]
    if not rate: livetime = 1
    else: livetime = float(lines[IndexFromSet(lines,timeset)+1].split()[0])
    spe_vals = [map(int,line.split())[0] for line in lines[data_start:data_stop]]
    spe_hist = ROOT.TH1F('spe_hist',spec_name.split()[0],len(spe_vals),0,len(spe_vals))
    spe_hist.FillN(len(spe_vals),numpy.array(range(0,len(spe_vals)),dtype='double'),
                                      numpy.array(spe_vals,dtype='double')/livetime)
    return spe_hist

#########################################
#########################################
