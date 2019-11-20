
from rootpy.plotting import Hist,Hist2D
from rootpy.plotting.views import ScaleView,StyleView
import matplotlib.pyplot as plt
from rootpy.io import File
try:
    from collections import OrderedDict
except ImportError:
    from ordered import OrderedDict
import logging; logging.basicConfig(level=logging.DEBUG)
from rootpy import log
from rootpy import asrootpy
from rounding import rounding
import re
from configobj import ConfigObj

log.basic_config_colorized()
class NoDictMessagesFilter(logging.Filter):
    def filter(self, record):
        #return "is not an exact divider of nbins=" not in record.msg
        return " is not an exact divider of nbins=" not in record.msg
log["/ROOT.TH1D.Rebin"].addFilter(NoDictMessagesFilter())
class LayoutMessagesFilter(logging.Filter):
    def filter(self, record):
        return "the on-file layout version 4 of class 'TAttAxis' differs from" not in record.msg
log["/ROOT.TStreamerInfo.CompareContent"].addFilter(LayoutMessagesFilter())


log_plotlib = log["plotlib"]

rnd=rounding()

## helper methods
#
# some times one needs float lists
def xfrange(start, stop, step):
    while start < stop:
        yield start
        start += step

# get n random colors (via hls)
def getColorList(ncolors):
    from colorsys import hls_to_rgb,rgb_to_hls
    from math import pi

    from random import uniform
    colors=[]
    for i in xfrange(0,2.*pi,(2.*pi)/ncolors):
        colors.append(hls_to_rgb(i, uniform(0.25,0.75) ,uniform(0.,1.) ))
    return colors

# helper to search if a hist contains a key
def getDictValue(hist,parmDict):
    if isinstance(parmDict,dict):
        for par in parmDict:
            if par in hist:
                return parmDict[par]
        return None
    if isinstance(parmDict,list):
        for par in parmDict:
            if par in hist:
                return True
        return False

def scale_xaxis(hist,factor):
    a=hist.GetXaxis()
    if (a.IsVariableBinSize()):
        import array
        tmpArray=a.GetXbins()
        bins=[tmpArray[i]*factor for i in range(tmpArray.GetSize())]
        a.Set(len(bins)-1, array("d",bins))
    else:
        a.Set(a.GetNbins(),a.GetXmin()*factor,a.GetXmax()*factor)
    return




# stolen from rootpy
# ignore_binns are hist if the bin values=0 the bins will be ignored
def duke_errorbar(hists,
             xerr=True, yerr=True,
             xpadding=0, ypadding=.1,
             xerror_in_padding=True,
             yerror_in_padding=True,
             emptybins=True,
             snap=True,
             axes=None,
             ignore_binns=None,
             **kwargs):
    from rootpy.plotting.hist import _Hist
    from rootpy.plotting.graph import _Graph1DBase
    from rootpy.plotting.root2matplotlib import _set_bounds

    """
    Make a matplotlib errorbar plot from a ROOT histogram or graph
    or list of histograms and graphs.

    Parameters
    ----------

    hists : Hist, Graph or list of Hist and Graph
        The histogram(s) and/or Graph(s) to be plotted

    xerr : bool, optional (default=True)
        If True, x error bars will be displayed.

    yerr : bool or string, optional (default=True)
        If False, no y errors are displayed.  If True, an individual y
        error will be displayed for each hist in the stack.  If 'linear' or
        'quadratic', a single error bar will be displayed with either the
        linear or quadratic sum of the individual errors.

    xpadding : float or 2-tuple of floats, optional (default=0)
        Padding to add on the left and right sides of the plot as a fraction of
        the axes width after the padding has been added. Specify unique left
        and right padding with a 2-tuple.

    ypadding : float or 2-tuple of floats, optional (default=.1)
        Padding to add on the top and bottom of the plot as a fraction of
        the axes height after the padding has been added. Specify unique top
        and bottom padding with a 2-tuple.

    xerror_in_padding : bool, optional (default=True)
        If True then make the padding inclusive of the x errors otherwise
        only pad around the x values.

    yerror_in_padding : bool, optional (default=True)
        If True then make the padding inclusive of the y errors otherwise
        only pad around the y values.

    emptybins : bool, optional (default=True)
        If True (the default) then plot bins with zero content otherwise only
        show bins with nonzero content.

    snap : bool, optional (default=True)
        If True (the default) then the origin is an implicit lower bound of the
        histogram unless the histogram has both positive and negative bins.

    axes : matplotlib Axes instance, optional (default=None)
        The axes to plot on. If None then use the global current axes.

    kwargs : additional keyword arguments, optional
        All additional keyword arguments are passed to matplotlib's errorbar
        function.

    Returns
    -------

    The return value from matplotlib's errorbar function, or list of such
    return values if a list of histograms and/or graphs was plotted.

    """
    if axes is None:
        axes = plt.gca()
    curr_xlim = axes.get_xlim()
    curr_ylim = axes.get_ylim()
    was_empty = not axes.has_data()
    if isinstance(hists, (_Hist, _Graph1DBase)):
        # This is a single plottable object.
        returns = _duke__errorbar(
            hists, xerr, yerr,
            axes=axes, emptybins=emptybins,
            ignore_binns=ignore_binns,
             **kwargs)
        _set_bounds(hists, axes=axes,
                    was_empty=was_empty,
                    prev_ylim=curr_ylim,
                    xpadding=xpadding, ypadding=ypadding,
                    xerror_in_padding=xerror_in_padding,
                    yerror_in_padding=yerror_in_padding,
                    snap=snap)
    else:
        returns = []
        for h in hists:
            returns.append(duke_errorbar(
                h, xerr=xerr, yerr=yerr, axes=axes,
                xpadding=xpadding, ypadding=ypadding,
                xerror_in_padding=xerror_in_padding,
                yerror_in_padding=yerror_in_padding,
                snap=snap,
                emptybins=emptybins,
                **kwargs))
    return returns

# stolen from rootpy
# ignore_binns are hist if the bin values=0 the bins will be ignored
def _duke__errorbar(h, xerr, yerr, axes=None, emptybins=True, ignore_binns=None, zorder=None, **kwargs):
    from rootpy.plotting.root2matplotlib import _set_defaults
    import numpy as np
    if axes is None:
        axes = plt.gca()
    if zorder is None:
        zorder = 1
    _set_defaults(h, kwargs, ['common', 'errors', 'errorbar', 'marker'])
    if xerr:
        xerr = np.array([list(h.xerrl()), list(h.xerrh())])
    if yerr:
        yerr = np.array([list(h.yerrl()), list(h.yerrh())])
    x = np.array(list(h.x()))
    y = np.array(list(h.y()))
    remove=[]
    if ignore_binns is not None:
        for i in range(len(x)-1,-1,-1):
            ignore=False
            for ihist in ignore_binns:
                if ihist[i+1].value==0:
                    ignore=True
            if ignore:
                remove.append(i)
    x=np.delete(x, remove)
    y=np.delete(y, remove)
    if xerr is not False:
        xerr=np.delete(xerr, remove,1)
    if yerr is not False:
        yerr=np.delete(yerr, remove,1)

    if not emptybins:
        nonempty = y != 0
        x = x[nonempty]
        y = y[nonempty]
        if xerr is not False:
            xerr = xerr[:, nonempty]
        if yerr is not False:
            yerr = yerr[:, nonempty]
    return axes.errorbar(x, y, xerr=xerr, yerr=yerr, zorder=zorder, **kwargs)


def getRGBTColor(color):
    import ROOT
    if type(color)==type(""):
        if "+" in color:
            tcolor,modif=color.split("+")
            color=getattr(ROOT,tcolor)+int(modif)
        elif "-" in color:
            tcolor,modif=color.split("-")
            color=getattr(ROOT,tcolor)-int(modif)
        else:
            color=getattr(ROOT,color)
    col=ROOT.gROOT.GetColor(color)
    return (col.GetRed(),col.GetGreen(),col.GetBlue())

def makeTeXTable(dir_name,hist,histContainer):
    rnd=rounding()
    #make a nice tex table
    mtBinOutTeX=open(dir_name+"/"+hist.replace("/","")+"output.tex","w+")


    head=r"""\documentclass[a4paper,landscape]{article}
\usepackage[utf8]{inputenc}
\usepackage[landscape]{geometry}
\title{}
\date{\today}
\begin{document}


"""
    mtBinOutTeX.write(head)
    mtBinOutTeX.write(r"\begin{table}[h]")
    mtBinOutTeX.write(r"\begin{tabular}{"+" ".join(["l" for i in range(len(histContainer.bg.hists)+len(histContainer.sg.hists)+3)])+r"}"+"\n")



    line=["$M^{min}_{T}$"]
    for bg in histContainer.bg.hists:
        line.append(bg)
    line.append("Background")
    line.append("Data")
    for sg in histContainer.sg.hists:
        line.append(sg)
    mtBinOutTeX.write(" & ".join(line) + r"\\"+"\n" )
    emptylines=0
    for ibin in range(1,histContainer.getData().GetNbinsX()):
        #mt bin
        line=[]
        line.append( "{0}".format(histContainer.getData().xedgesl(ibin)))
        for bg in histContainer.bg.hists:
            line.append( rnd.latex(histContainer.bg.hists[bg].integral(xbin1=ibin)))
        integ,err=histContainer.bg.getAllAdded().integral(xbin1=ibin,error=True)
        line.append(rnd.latex(integ,err))
        line.append(rnd.latex(histContainer.getData().integral(xbin1=ibin)))
        for sg in histContainer.sg.hists:
            line.append( rnd.latex(histContainer.sg.hists[sg].integral(xbin1=ibin)))
        if (histContainer.getData().integral(xbin1=ibin)+histContainer.bg.getAllAdded().integral(xbin1=ibin)) < 0.01:
            emptylines+=1
        else:
            mtBinOutTeX.write(" & ".join(line) + r"\\"+"\n" )
        if emptylines>0:
            break
    mtBinOutTeX.write("\n\n")


    mtBinOutTeX.write(r"\end{tabular}"+"\n")
    mtBinOutTeX.write(r"\caption{}"+"\n")
    mtBinOutTeX.write(r"\label{tab:}"+"\n")
    mtBinOutTeX.write(r"\end{table}"+"\n")
    mtBinOutTeX.write(r"\begin{table}[h]"+"\n")
    mtBinOutTeX.write(r"\begin{tabular}{"+" ".join(["l" for i in range(len(histContainer.bg.hists)+len(histContainer.sg.hists)+3)])+r"}"+"\n")

    line=["$M_{T}$"]
    for bg in histContainer.bg.hists:
        line.append(bg)
    line.append("Background")
    line.append("Data")
    for sg in histContainer.sg.hists:
        line.append(sg)
    mtBinOutTeX.write(" & ".join(line) + r"\\"+"\n" )
    emptylines=0
    for ibin in range(1,histContainer.getData().GetNbinsX()):
        #mt bin
        line=[]
        line.append( "{0} - {1}".format(histContainer.getData().xedgesl(ibin),histContainer.getData().xedgesh(ibin)))
        for bg in histContainer.bg.hists:
            line.append( rnd.latex(histContainer.bg.hists[bg][ibin].value))
        line.append(rnd.latex(histContainer.bg.getAllAdded()[ibin].value,histContainer.bg.getAllAdded()[ibin].error))
        line.append(rnd.latex(histContainer.getData()[ibin].value))
        for sg in histContainer.sg.hists:
            line.append( rnd.latex(histContainer.sg.hists[sg][ibin].value))
        if (histContainer.getData()[ibin].value+histContainer.bg.getAllAdded()[ibin].value) < 0.01:
            emptylines+=1
        else:
            mtBinOutTeX.write(" & ".join(line) + r"\\"+"\n" )
        if emptylines>5:
            break
        tail="""
\end{tabular}
\caption{}
\label{tab:}
\end{table}


\end{document}
"""
    mtBinOutTeX.write(tail)
    mtBinOutTeX.close()

##@class HistSorageContainer Class to handle data, bg and sg HistStorages
#
# written by Klaas Padeken 2015
class HistStorageContainer():
    ## Init function
    #
    # In this function the default variables are set and initialized.
    # @param[in] bg HistSorage for bg (default=False)
    # @param[in] sg HistSorage for sg (default=False)
    # @param[in] data HistSorage for data (default=False)
    def __init__(self,bg=False,sg=False,data=False):
        self.bg=bg
        self.sg=sg
        self.data=data
        self.allStored=[]
        if bg is not False:
            self.allStored.append(self.bg)
            self.bg.initStyle(style="bg")
        if sg is not False:
            self.allStored.append(self.sg)
            self.sg.initStyle(style="sg")
        if data is not False:
            self.allStored.append(self.data)
            self.data.initStyle(style="data")

    ## Function to clear hists
    #
    # use this if you want to plot a new set of hists
    #def clearHists(self):
        #for stored in self.allStored:
            #stored.clearHists()

    ## Function get hists from files
    #
    # the hists ate added to .hists and joined if a joinList exist
    # the hist will be put in the corresponding HistStorage class
    # @param[in] hist string name of the hist
    def getHist(self, hist, backupHist = None, noScale = False):
        for stored in self.allStored:
            stored.getHist(hist, backupHist = backupHist, noScale = noScale)

    def getHistFromTree(self,binns,xmin,xmax,xtitle,cut,value,tree,weight=None):
        for stored in self.allStored:
            stored.getHistFromTree(binns,xmin,xmax,xtitle,cut,value,tree,weight=weight)


    def getHistFromTree2d(self,xbinns,xmin,xmax,xtitle,ybins,ymin,ymax,ytitle,cut,value,tree,weight=None):
        for stored in self.allStored:
            stored.getHistFromTree2d(xbinns,xmin,xmax,xtitle,ybins,ymin,ymax,ytitle,cut,value,tree,weight=weight)

    ## Function rebin the all hists
    #
    # @param[in] width try to rebin to a specific width
    # @param[in] factor rebin to with a factor
    # if both are given the width is used
    def rebin(self,width=0,factor=0,vector=None):
        for stored in self.allStored:
            stored.rebin(width=width,factor=factor,vector=vector)

    ## Function scale xaxis by factor
    #
    # @param[in] factor to rescale
    def scale_xaxis(self,factor):
        for stored in self.allStored:
            stored.scale_xaxis(factor)

    ## Function setTitle for the all hists
    #
    # @param[in] xtitle
    def setTitle(self,xtitle,ytitle=""):
        for stored in self.allStored:
            stored.setTitle(xtitle,ytitle)

    ## Function to set the style of the histograms
    #
    # @param[in] bgcolors a list/dict of colors that the hists should have
    # if colors is not specified the internal colorListis used if set
    # @param[in] sgcolors a list/dict of colors that the hists should have
    # if colors is not specified the internal colorListis used if set
    def setStyle(self,bgcolors=None,sgcolors=None):
        if self.bg is not False:
            self.bg.setStyle()
            if bgcolors!=None:
                self.bg.colorList=bgcolors
        if self.sg is not False:
            self.sg.setStyle()
            if sgcolors!=None:
                self.sg.colorList=sgcolors
        if self.data is not False:
            self.data.setStyle()

    ## Function getbg
    #
    # @param[out] bg.getHistList() list of bg hists
    def getBGList(self):
        return self.bg.getHistList()

    ## Function getsg
    #
    # @param[out] sg.getHistList() list of sg hists
    def getSGList(self):
        return self.sg.getHistList()

    ## Function getData
    #
    # @param[out] data.getHistList()[0] hist of Data
    def getData(self):
        if len(self.data.getHistList())>1:
            return sum(self.data.getHistList())
        return self.data.getHistList()[0]

    ## Function makeCumulative
    #
    # make the distribution cumulative
    def makeCumulative(self,width=False):
        for h in self.allStored:
            h.makeCumulative(width=width)








##@class HistStorage
# Class to handle histograms functions
#
# To store a list of hists, which are scaled, joined, rebinned and otherwise
# manipulated.
#
# written by Klaas Padeken 2015
class HistStorage(object):
    ## Init function
    #
    # In this function the default variables are set and initialized.
    # @param[in] xs is a xs ConfigObj needed for scaling
    # @param[in] lumi is the lumi in pb
    # @param[in] path is the default path of the files (default=None)
    # @param[in] isData is a switch  (default=None)
    def __init__(self,xs=1., lumi=1.,xstype="pythonConfig", path=None,isData=False,useRoot=False):
        self.views=OrderedDict()
        self.hists=OrderedDict()
        self.files=OrderedDict()
        self.colorList={}
        self.genNumber={}
        self.weight={}
        self.basepath=path
        self.verbosity=3
        self._scaled=False
        self.datadrivenHist=[]
        self.xs=xs
        self.configxs=xstype
        self.lumi=lumi
        self.isData=isData
        self._joinList=False
        self.style={}
        self._style_changed=False
        self.isCumulative=False
        self.forcedWidth=False
        self.Unit=""
        self.matplotlibStyle= (not useRoot)
        self.additionalWeight={}
        self.eventString="Events"

    ## del function
    #
    # This deletes the main objects nedded to not get a crash at the end!
    #def __del__(self):
        ##for name in self.views:
            ##del self.views[name]
        ##for name in self.hists:
            ##del self.hists[name]
        #for name in  self.files:
            #self.files[name].Close()
    ##------------------------------------------------------------------
    ## Private functions
    ##------------------------------------------------------------------
    ## Function to get the event numbers from h_counter in file!
    #
    # The function fills the dict genNumber with the event numbers.
    def _getGenNumbers(self):
        for name in self.files:
            if name in self.genNumber or name in self.datadrivenHist or self.isData:
                continue
            try:
                if self.configxs=="BSM3G":
                    for path, dirs, objects in self.files[name].walk():
                        if "Events" in objects:
                            evHist=self.files[name].Get(path+"/Events")
                            Nev=evHist[1].value
                            break
                else:
                    counter=self.files[name].Get("h_counters")
                    if(counter[2].value>0):
                        Nev=counter[2].value
                    else:
                        Nev=counter[1].value
            except:
                if self.verbosity==3:
                    log_plotlib.info("[Info] If you want to use Nev from the root file store 'h_counters'")
                    log_plotlib.info("will set to 1 for %s"%(name))
                Nev=1
            self.genNumber[name]=Nev

    ## Function get the unit from xaxis.
    #
    # returns the unit of the hist
    def _getUnit(self):
        if self.Unit!="":
            return self.Unit
        xtitle=list(self.hists.values())[0].GetXaxis().GetTitle()
        noUnit=["phi","eta","_{T}/E^{miss}"]
        for veto in noUnit:
            if veto in xtitle:
                return ""
        t=""
        # this should cover all the usual cases, like [unit] /unit or *eV.
        prossibleUnits=['[\[\(]\S*[\]\)]','/\S*\}','\SeV']
        for unit in prossibleUnits:
            m = re.search(unit, xtitle)
            if m != None:
                break
        if m != None:
            t=m.group(0).translate(None,"[\[\(\]\)]/\}")
        self.Unit=t
        return t

    ## Function to add files to a scaled view.
    #
    # The scaled view dict "views" now retruns all histograms scaled!
    def _addToScaledView(self):
        for name in self.files:
            if name in self.views:
                continue
            weight=self._getWeight(name)
            self.views[name]=ScaleView(self.files[name],weight)
            self._scaled=True

    def _getWeight(self,name):
        if name in self.weight:
            return self.weight[name]
        else:
            if name in self.datadrivenHist or self.isData:
                weight=1.
            else:
                if self.configxs=="pythonConfig":
                    if "xs" in self.xs[name]:
                        if "weight" in self.xs[name]:
                            weight=self.xs[name].as_float("xs")*self.xs[name].as_float("weight")*self.lumi/self.genNumber[name]
                        else:
                            weight=self.xs[name].as_float("xs")*self.lumi/self.genNumber[name]
                    else:
                        if "weight" in self.xs[name]:
                            weight=self.xs[name].as_float("crosssection")*self.xs[name].as_float("weight")*self.lumi/self.genNumber[name]
                        else:
                            weight=self.xs[name].as_float("crosssection")*self.lumi/self.genNumber[name]
                elif self.configxs=="music":
                    weight=self.xs.as_float("%s.XSec"%(name))*self.xs.as_float("%s.FilterEff"%(name))*self.xs.as_float("%s.kFactor"%(name))*self.xs.as_float("Lumi")/self.genNumber[name]
                elif self.configxs=="music_scaled":
                    raise Exception('not implemented music n_files getter yet!!')
                    weight=1/self._getNjobs()
                elif self.configxs=="BSM3G":
                    _name=name
                    if name not in self.xs:
                        _name=name+".root"
                    weight=self.xs[_name]["xs"]*self.xs[_name]["eff"]*self.lumi/self.genNumber[name]
                elif self.configxs==None:
                    weight=1.
            if name in self.additionalWeight:
                weight*=self.additionalWeight[name]
            self.weight[name]=weight
            return weight

    def __getitem__(self,item):
        return self.hists[item]


    ## Function to add all files in the given path
    #
    # Use setPath(path) to set the path if did not in the init.
    # @param[in] tag if regexpr is not used all *.root files containing the tag are added
    # @param[in] veto define a !list!! of veto strings not case sensitive
    # @param[in] regexpr use a regular expression to find the file names (need .root at the end if
    # you want to use root files!!
    # @param[in] joinName if specified all files matching the expressions above will be added to the list of files that should be joined.
    def addAllFiles(self,tag="",veto=None, regexpr=None, joinName=None):
        if self.basepath==None:
            raise RuntimeError("You must set a basepath to add all files from one directory!")
        if regexpr is not None:
            import re,os
            fileList = [f for f in os.listdir(self.basepath+"/") if re.search(r'%s'%(regexpr), f)]
        else:
            import glob
            fileList=glob.glob(self.basepath+"/*"+tag+"*.root")
        tmpList=[]
        for file in fileList:
            if veto is not None:
                vetoed=False
                for v in veto:
                    if v.lower() in file.split("/")[-1].lower():
                        vetoed=True
                if vetoed:
                    log_plotlib.info("vetoed file: %s"%file)
                    continue
            name=file.split("/")[-1].replace(".root","")
            self.files[name]=File(file, "read")
            tmpList.append(name)
        self._getGenNumbers()
        self._addToScaledView()
        if joinName is not None:
            if self._joinList is not False:
                self._joinList[joinName]=tmpList
            else:
                self._joinList=OrderedDict()
                self._joinList[joinName]=tmpList

    ## Function to add a single file
    #
    # Use setPath(path) to set the path if did not in the init.
    # @param[in] name the name of the file that should be added!
    def addFile(self,name):
        self.files[name]=File(self.basepath+"/"+name+".root", "read")
        self._getGenNumbers()
        self._addToScaledView()


    ## Function to remove a single file
    #
    # Use setPath(path) to set the path if did not in the init.
    # @param[in] name the name of the file that should be added!
    def removeFile(self,name):
        del self.files[name]
        del self.views[name]
        self._joinList = {key: value for key, value in list(self._joinList.items()) if name not in value}

    ## Function to add files specified as a list or (ordered)dict
    #
    # Use setPath(path) to set the path if did not in the init.
    # @param[in] fileList list or dict of the files you want to add
    # if the dict is used the files are joined to a single hist with this key
    def addFileList(self,fileList):
        if type(fileList)==type(list()):
            for file in fileList:
                self.files[file]=File(self.basepath+"/"+file+".root", "read")
        if isinstance(fileList,dict):
            import itertools
            useList=list(itertools.chain.from_iterable(list(fileList.values())))
            for file in useList:
                try:
                    self.files[file]=File(self.basepath+"/"+file+".root", "read")
                except:
                    yesno=input("file %s is not there continue? [y/n]"%(file))
                    if yesno!="y":
                        import sys
                        sys.exit(1)
                    fileList={key: value for key, value in list(fileList.items())
                        if file not in value}
            self._joinList=fileList
        self._getGenNumbers()
        self._addToScaledView()

    ## Function to add a dict to join files
    #
    # @param[in] joinList wich should be a (ordered)dict
    def addJoinList(self,joinList):
        self._joinList=joinList

    ## Function to apply a style to a single histogram
    #
    # @param[in] name string of the file
    # @param[in] kwargs dict all the styles can be set like fillstyle = 'solid'
    def applyStyle(self, name, **kwargs):
        if name in self.style:
            self.style[name].update(kwargs)
        else:
            self.style[name]=kwargs.copy()

    ## Function to apply a style to all histograms
    #
    # @param[in] kwargs all the styles can be set like fillstyle = 'solid'
    def applyStyleAll(self, **kwargs):
        iteratorList=self.files
        if self._joinList is not False:
            iteratorList=self._joinList

        for view in iteratorList:
            if view in self.style:
                self.style[view].update(kwargs)
            else:
                self.style[view]=kwargs.copy()

    ## Function to clear hists
    #
    # use this if you want to plot a new set of hists
    def clearHists(self):
        self.Unit=""
        self.forcedWidth=False
        self.isCumulative=False
        del self.hists
        self.hists=OrderedDict()

    ## Function to get a hist that is the sum of hists in the storage
    #
    # handy if you want only a subgroup as a hist
    # @param[in] name add only files that contain the name (default="")
    # @param[in] ignoreScale if you want to add hists that are not scaled (default=False)
    # @param[out] Hist
    def getAdded(self,name="",ignoreScale=False):
        if not self._scaled and not ignoreScale:
            raise RuntimeError("Add all histograms without scaling. I think not!")
        temp = []
        for key in self.hists:
            if name in key:
                temp.append( self.hists[key] )
        self.setStyle()
        return sum(temp)

    ## Function to get a hist that is the sum of all hists in the storage
    #
    # same as getAdded() perhaps faster
    # @param[in] ignoreScale if you want to add hists that are not scaled (default=False)
    # @param[out] Hist
    def getAllAdded(self,ignoreScale=False):
        if not self._scaled and not ignoreScale:
            raise RuntimeError("Add all histograms without scaling. I think not!")
        self.setStyle()
        return sum(self.hists.values())

    ## Function get hists from files
    #
    # the hists are added to .hists and joined if a joinList exist
    # @param[in] hist string of the hist in the files
    def getHist(self, hist, backupHist = None, noScale = False):
        self.clearHists()
        for f in self.views:
            try:
                if not "/eff_" in hist and not noScale:
                    self.hists[f]=self.views[f].Get(hist)
                    self.hists[f].Sumw2()
                else:
                    self.hists[f]=self.files[f].Get(hist)
            except Exception as e:
                if backupHist != None:
                    log_plotlib.warning( "No %s in %s, using %s instead"%(hist,f,backupHist))
                    self.hists[f]=self.views[f].Get(backupHist)
                    self.hists[f].Sumw2()
                else:
                    log_plotlib.warning( "No %s in %s (error:%s)"%(hist,f,e))
                    if len(self.hists)>0:
                        self.hists[f]=self.hists.values()[0].clone()
                        self.hists[f].Reset()
                    else:
                        self.hists[f]=Hist(100,0,100)

        if self._joinList is not False:
            self.joinList(self._joinList)
        for hist in self.hists:
            if hist in self.style:
                self.hists[hist].decorate(**self.style[hist])



    ## Function get hists via trees from files
    #
    # the hists ate added to .hists and joined if a joinList exist
    # @param[in] hist string of the hist in the files
    def getHistFromTree(self,bins,xmin,xmax,xtitle,cut,value,tree,weight=None):
        self.clearHists()
        for f in self.files:
            try:
                _tree=self.files[f].Get(tree)
            except AttributeError as e:
                log_plotlib.warning( "No %s in %s"%(tree,f))
                log_plotlib.warning( "Will try without %s, and add an empty hist."%f)
                log_plotlib.warning( e)
                self.hists[f]=Hist(binns,xmin,xmax)
                continue
            self.hists[f]=Hist(bins,xmin,xmax)
            self.hists[f].GetXaxis().SetTitle(xtitle)
            try:
                if weight is None:
                    _tree.Draw(value,selection=cut,hist=self.hists[f])
                else:
                    #_tree.Draw(value,selection="(%s)*(%s)"%(cut,weight),hist=self.hists[f])
                    tmpFile=File("/tmp/tmp.root", "recreate")
                    #sel_tree=_tree.copy_tree(selection=cut)
                    sel_tree=asrootpy(_tree.CopyTree(cut))
                    ##print weight
                    sel_tree.Draw(value,selection=weight,hist=self.hists[f])
                    tmpFile.Close()
            except Exception as e:
                log_plotlib.info( "error:%s"%(e))
                log_plotlib.info( "file :%s"%(f))
                log_plotlib.info( "Perhaps try this one:")
                for i in _tree.glob("*"):
                    log_plotlib.info( i)
                raise RuntimeError("Will stop here!")
            self.hists[f].Scale(self._getWeight(f))
            self.hists[f].Sumw2()
        if self._joinList is not False:
            self.joinList(self._joinList)
        for hist in self.hists:
            if hist in self.style:
                self.hists[hist].decorate(**self.style[hist])


    ## Function get hists via trees from files
    #
    # the hists ate added to .hists and joined if a joinList exist
    # @param[in] hist string of the hist in the files
    def getHistFromTree2d(self,xbins,xmin,xmax,xtitle,ybins,ymin,ymax,ytitle,cut,value,tree,weight=None):
        self.clearHists()
        for f in self.files:
            try:
                _tree=self.files[f].Get(tree)
            except AttributeError as e:
                log_plotlib.warning( "No %s in %s"%(tree,f))
                log_plotlib.warning( "Will try without %s, and add an empty hist."%f)
                log_plotlib.warning( e)
                self.hists[f]=Hist(binns,xmin,xmax)
                continue
                #TH2F(const char* name, const char* title, Int_t nbinsx, const Float_t* xbins, Int_t nbinsy, const Float_t* ybins)
                #TH2F(const char* name, const char* title, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup)
            self.hists[f]=Hist2D(xbins,xmin,xmax,ybins,ymin,ymax)
            self.hists[f].GetXaxis().SetTitle(xtitle)
            self.hists[f].GetXaxis().SetTitle(ytitle)
            try:
            #if weight is None:
                _tree.Draw(value,selection=cut,hist=self.hists[f])
            #else:
                #tmpFile=File("tmp.root", "recreate")
                #sel_tree=_tree.copy_tree(selection=cut)
                #print weight
                #sel_tree.Draw(value,selection=weight,hist=self.hists[f])
            except:
                log_plotlib.info( "Perhaps try this one:")
                for i in _tree.glob("*"):
                    log_plotlib.info( i)
                raise RuntimeError("Will stop here!")
            self.hists[f].Scale(self._getWeight(f))
            self.hists[f].Sumw2()
        if self._joinList is not False:
            self.joinList(self._joinList)
        for hist in self.hists:
            if hist in self.style:
                self.hists[hist].decorate(**self.style[hist])

    ## Function to get the hists as a list
    #
    # @param[out] list of all stored hists
    def getHistList(self):
        self.setStyle()
        return list(self.hists.values())

    ## Function join files containing name to one nameed label
    #
    # @param[in] name add all files containing name
    # @param[in] label name of the resulting new hist
    def join(self,name,label):
        if name == "*" or name == "":
            # join all bg:
            joined = self.getAllAdded()
            self.hists=OrderedDict()
            self.hists[label]=joined
        else:
            # join hist filtered by name:
            joined = self.getAdded(name)
            for key in self.hists:
                if name in key:
                    self.hists.pop(key)
            self.hists[label]=joined

    ## Function join files via a (ordered)dict
    #
    # @param[in] joinList add all files that are in the (ordered)dict to one hist with the name of the key
    def joinList(self,joinList):
        #import ROOT

        for name in joinList:
            self.hists[name]=self.hists[joinList[name][0]]
            self.hists.pop(joinList[name][0])
            for h in joinList[name][1:]:
                try:
                    self.hists[name]+=(self.hists[h])
                except:
                    self.hists[name].Add(self.hists[h])
                self.hists.pop(h)

    ## Function rebin the all hists
    #
    # @param[in] width float try to rebin to a specific width
    # @param[in] factor float rebin to with a factor
    # @param[in] vector list of all binns
    # if both are given the width is used
    def rebin(self,width=0,factor=0,vector=None):
        if vector != None:
            for name in self.hists:
                rebinnedHist=self.hists[name].rebinned(vector)
                rebinnedHist.xaxis.SetTitle(self.hists[name].xaxis.GetTitle())
                self.hists[name]=rebinnedHist
                if name in self.style:
                    self.hists[name].decorate(**self.style[name])
                if width!=0:
                    for ibin in self.hists[name]:
                        ibin.value=ibin.value/(ibin.x.width/width)
                        ibin.error=ibin.error/(ibin.x.width/width)
                    self.forcedWidth=width
        else:
            if width!=0:
                factor=int(width/list(self.hists.values())[-1].xwidth(1)+0.5)
            for name in self.hists:
                self.hists[name].Rebin(factor)


    ## Function scale xaxis by factor
    #
    # @param[in] factor to rescale
    def scale_xaxis(self,factor):
        for name in self.hists:
            scale_xaxis(self.hists[name],factor)


    ## Function to make the hist cumulative
    #
    # @param[in] width if specified the bins are specified the bins are corrected for the width
    def makeCumulative(self,width=False):
        for hist in self.hists:
            for ibin in self.hists[hist].bins():
                if width is not False:
                    #hmm does one realy want this I think the error is wrong
                    ibin.value,ibin.error=self.hists[hist].integral(xbin1=ibin.idx, error=True)
                    ibin.value/=(float(width)/self.hists[hist].xwidth(ibin.idx))
                    ibin.error/=(float(width)/self.hists[hist].xwidth(ibin.idx))
                else:
                    ibin.value,ibin.error=self.hists[hist].integral(xbin1=ibin.idx, error=True)
        self.isCumulative=True
        self.forcedWidth=width

    ## Function to set the datadiven name flag
    #
    # @param[in] ddhist the name of the datadriven hist
    def setDataDriven(self,ddhist):
        self.datadrivenHist.append(ddhist)

    ## Function to set path of the hists
    #
    # @param[in] path
    def setPath(self,path):
        self.basepath=path

    ## Function to set x axis title for all hists
    #
    # @param[in] xtitle
    def setTitle(self,xtitle,ytitle=""):
        for h in list(self.hists.values()):
            h.GetXaxis().SetTitle(xtitle)
            if ytitle!="":
                h.GetYaxis().SetTitle(ytitle)


    ## Function to set the style of the histograms
    #
    # sets the axis labels and titles
    def setStyle(self):
        for key in self.hists:
            try:
                if self.matplotlibStyle:
                    if "$" not in self.hists[key].xaxis.GetTitle():
                        self.hists[key].xaxis.SetTitle("$\\mathrm{"+self.hists[key].xaxis.GetTitle().replace("#","\\")+"}$")
                        #self.hists[key].xaxis.SetTitle("${"+self.hists[key].xaxis.GetTitle().replace("#","\\")+"}$")
                if self.isCumulative and ">" not in self.eventString:
                    self.eventString+=">%s"%(self.hists[key].xaxis.GetTitle().translate(None,self._getUnit()+"[]/()"))
                if self.forcedWidth is not False:
                    width=self.forcedWidth
                else:
                    width=list(self.hists.values())[-1].xwidth(2)
                self.hists[key].yaxis.SetTitle("%s/(%s %s)"%(self.eventString,rnd.latex(width),self._getUnit()))
                if self.isData:
                    self.hists[key].SetTitle("Data")
                else:
                    self.hists[key].SetTitle(key)
            except:
                log_plotlib.warning("Could not change the axis title")

    ## Function to init the style of the histograms
    #
    # @param[in] style "bg" and "sg" posible (default="bg")
    # @param[in] colors a list/dict of colors that the hists should have
    # if colors is not specified the internal colorListis used if set
    def initStyle(self,style="bg",colors=None):
        if style=="bg":
            self.applyStyleAll(fillstyle = 'solid',linewidth = 0)
        if style=="sg":
            self.applyStyleAll(fillstyle = '0',linewidth = 2 )
        #nothing to do here yet
        if style=="data":
            self.applyStyleAll(markersize=0.5)
            pass

        iteratorList=self.views
        if self._joinList is not False:
            iteratorList=self._joinList

        for key in iteratorList:
            if colors!=None:
                if isinstance(colors, (list)):
                    usecolor=colors.pop()
                    self.applyStyle(key,fillcolor = usecolor,linecolor = usecolor)
                if isinstance(colors, (dict)):
                    if key in colors:
                        self.applyStyle(key,fillcolor = colors[key],linecolor = colors[key])
            else:
                if isinstance(self.colorList, list):
                    if len(self.colorList)>0:
                        usecolor=self.colorList.pop()
                        self.applyStyle(key,fillcolor = usecolor,linecolor = usecolor)
                if isinstance(self.colorList, dict):
                    if key in self.colorList:
                        self.applyStyle(key,fillcolor = self.colorList[key],linecolor = self.colorList[key])


