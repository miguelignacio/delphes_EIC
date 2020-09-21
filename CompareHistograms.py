#!/usr/bin/env python
from AtlasCommonUtils import *


def my_callback(option, opt, value, parser):
  setattr(parser.values, option.dest, value.split(','))

import optparse,os
op = optparse.OptionParser(usage="%prog [opts] ")
op.add_option("","--outputDir",dest="outdir",
              default=".",
              help="Output Directory")
op.add_option("-i","--inputFileList",dest="inlist",
              default='input.txt',
              help="Input File list with xsections")
op.add_option("-n","--Normalized",dest="Normalized",
              default=None,
              help="Normalized")

op.add_option("-l","--RatioLimit",dest="RatioLimit",
              default=4,
              help="The limit for the ratio plot (in %)")
op.add_option("-t","--Tag",dest="Tag",
              default="",
              help="tag to add to plots")

op.add_option("","--hists",type='string', dest = "hists",
              action='callback', callback=my_callback,
              help="comma separated list of histograms you would like to Compare, e.g. hists,hists_i,hists_o,hists_ss" )


(p_ops, p_args) = op.parse_args()


output_dir = p_ops.outdir
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
else:
    print "Using output dir ",output_dir,", overwriting existing contents"


if not p_ops.inlist:
    raise RuntimeError("Must specify infile list")

hists = []
if not p_ops.hists:
    hists = ["nhits"]
else:
    hists = p_ops.hists



RatioLimit = float(p_ops.RatioLimit)
Tag = p_ops.Tag

print 'The ratio plot will be plotted at maximum' , RatioLimit
import ROOT 
from Legend import Legend
from ROOT import TLatex
SetAtlasStyle()


color = [1,2,4,6,7,8]
style = [20,20,20,26] 

def DivideGraphs(den,num): 
    x = ROOT.Double() 
    y1= ROOT.Double()
    y2= ROOT.Double()
    ratio = ROOT.TGraphErrors()
    for i in range(den.GetN()):
        num.GetPoint(i,x,y2)
        den.GetPoint(i,x,y1)
        if(y1!=0):
            #print "SETTING POINT ", y2/y1
            ratio.SetPoint(i,x,(y2/y1-1)*100)
            ratio.SetPointError(i,0,0)
    #num.Print()
    #den.Print()
    #ratio.Print()
    return ratio
     

def plotHists(hists,files, xtitle, ytitle, title, pdfname,ly=False):
    c = ROOT.TCanvas()
    start = 0
    label = Legend("")
    multi = ROOT.TMultiGraph()

    for i in range(len(hists)): #This loops over the files. nhists(f1,f2,f3) is the format of hists.
        temp = files[i].GetName()
        #print 'printing Temp ' , temp
        temp = temp.split('_MinBias')[0]
        if temp in tag :
            temp = temp + " , version %s " %(i+1)		
        hists[i].SetLineColor(color[i])
        hists[i].SetMarkerColor(color[i])
        hists[i].SetMarkerStyle(style[i])
        hists[i].SetMarkerSize(1)
        hists[i].SetLineWidth(2)
       # print 'i ' , i , ' MINIMUM:  ' , hists[i].GetHistogram().GetMinimum() , ' MAXIMUM : ' , hists[i].GetHistogram().GetMaximum(), 
        label.Add(hists[i],tag[i],"L")
        
      
        multi.Add(hists[i])

    

    multi.Draw("AP")
    ROOT.gStyle.SetErrorX(0.0001);
    #multi.SetMaximum(1.02*multi.GetHistogram().GetMaximum())
    multi.SetMaximum(0.25)
    multi.SetMinimum(0.0)
    multi.SetTitle(title)
    multi.GetXaxis().SetTitle(xtitle)
    multi.GetYaxis().SetTitle(ytitle)
    multi.GetXaxis().SetNdivisions(10)

    label.Draw(0.2,.95)

    if ly==True:
        ROOT.gPad.SetLogy(1)
    latex = TLatex()
    latex.SetNDC()
    #latex.DrawLatex(0.2,0.90,"#font[62]{#scale[0.9]{lucidEvtOr}}")
    c.SaveAs(pdfname)
    c.Clear()
    print "END OF PLOTTING FUNCTION " 

def get_hists(files,hists):
    allhists = []
    for h in hists:
        hlist = []
        for f in files:
            f.cd()
            hlist.append(f.Get(h))
        allhists.append(hlist)
    return allhists


files = []
tag = []
xsecs = []
legs = []
weight = []

print "Abriendo el file " 
inlist = open(p_ops.inlist)

for line in inlist:
    if line[0] == "#":
       continue
    print line
    tag.append(line.split(",")[1])
    weight.append(1)
    line = line.split(",")[0]
    print 'Opening Line', line
    files.append(ROOT.TFile(line.strip("\n")))
    line = line.split("hists")[0].split("_")[-1]
    tagOutputPlots = line
    


print "Getting all hists"
c = ROOT.TCanvas()

allhists = get_hists(files, hists)
temp = get_hists(files,hists)
for i in range(len(hists)):
    plotHists(temp[i], files, allhists[i][0].GetXaxis().GetTitle(),allhists[i][0].GetYaxis().GetTitle() ,  allhists[i][0].GetTitle(),  output_dir + "/"+hists[i]+Tag+'_Comparison.pdf',ly=False)

c = ROOT.TCanvas()
c.SetGrid(0,1)
#allhists[0][0].Print()
zero = ROOT.TF1('zero','[0]+[1]*x',-5,4000)
zero.SetLineColor(2)

g_ratio = []
multi = ROOT.TMultiGraph()
label = Legend("")

print len(tag)
if len(tag)>1 :
    for i in range(len(hists)):
        g_ratio.append(DivideGraphs(allhists[i][0],allhists[i][1]))
        #zero.SetLineColor(1)
        #g_ratio[0].Fit(zero)
        g_ratio[0].SetMarkerSize(1)
        g_ratio[0].SetMarkerColor(color[1])
        g_ratio[0].Print()
        print 'Mean of ratios ' , '%.3f' %g_ratio[0].GetMean(2) , ' %,  RMS of ratios ' ,  '%.3f' %g_ratio[0].GetRMS(2) , " % " 
        multi.Add(g_ratio[0])
        #label.Add(g_ratio[0],tag[1]+"/"+tag[0])
        for j in range(2,len(tag)):
            g_ratio.append(DivideGraphs(allhists[i][0],allhists[i][j]))
            g_ratio[j-1].SetMarkerColor(color[j])
            g_ratio[j-1].SetMarkerSize(1)
            #zero.SetLineColor(color[j])
            #g_ratio[j-1].Fit(zero)
            
            multi.Add(g_ratio[j-1])
           # label.Add(g_ratio[j-1],tag[j]+"/"+tag[0])

      
        multi.Draw("AP")
        label.Draw(0.5,.35)
        multi.SetMaximum(RatioLimit)
        multi.SetMinimum(-1*RatioLimit)
        multi.GetYaxis().SetTitle("Ratio -1 [%]")
        multi.GetXaxis().SetTitle(allhists[i][0].GetXaxis().GetTitle())
        #latex = TLatex()
        #latex.SetNDC()
        #latex.DrawLatex(0.5,0.90,tag[1]+"/" + tag[0])              

        #c.SaveAs(output_dir +"/"+hists[i]+"_"+Tag+"Ratios.pdf")






