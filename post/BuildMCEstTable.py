#!/usr/bin/env python
#coding=utf-8

import sys
import os
import numpy as np
import matplotlib.pylab as plt
import statsmodels.api as sm
from matplotlib.ticker import FuncFormatter

from TableTools import *
from StatedepTools import *

################################################################################
# HARD-CODING
################################################################################
FNOUT = 'TableMCEst.tex'
NBASE = 3435
NPLOTLIST = [1,2] # Smallest to largest, numbered from 0
NPLOTCOLORLIST = ['gray', 'black']
ARROWPERCENTBUFFERY = .1 # ???
ARROWPERCENTBUFFERX = .3
ARROWSIZE = 8
PERCENTBUFFERY = .05
PLOTGRIDDIM = 1000
PLOTLINEWIDTH = 2.0

################################################################################
################################################################################
################################################################################
CodeDir = os.path.dirname(os.path.abspath(sys.argv[0]))
ResultsDir = os.path.abspath(sys.argv[1])
if not os.path.isdir(ResultsDir):
    print ('Could not find directory ' + ResultsDir)
    sys.exit()

NDIRLIST = [dirname for dirname in os.listdir(ResultsDir)
            if os.path.isdir(os.path.join(ResultsDir, dirname))]
NDIRLIST.sort()
NLIST = [float(dirname.lstrip('N')) for dirname in NDIRLIST]

# Load data
LB = {}
UB = {}
MinCrit = {}
FirstFlag = True
for n in NDIRLIST:
    LB[n] = np.genfromtxt(   os.path.join(ResultsDir, n, FNLB),\
                                delimiter=',', names=True)
    UB[n] = np.genfromtxt(   os.path.join(ResultsDir, n, FNUB),\
                                delimiter=',', names=True)
    assert LB[n].dtype.names == UB[n].dtype.names
    if FirstFlag:
        (TrueLB, TrueUB) = createBoundsDataFrame(\
                [os.path.join(ResultsDir, n, FNTRUEBOUNDS)])
        ParamNames = tuple(TrueLB.index.values)
        FirstFlag = False
    else:
        assert LB[n].dtype.names == ParamNames

    MinCrit[n] = np.loadtxt(os.path.join(ResultsDir, n, FNMINCRITERION))

# Table specification and header
fout = open(os.path.join(ResultsDir, FNOUT), 'w')
colspec = 'cl' + (1 + 2*len(NLIST))*'c'
startTable(fout, colspec)
insertTopRule(fout)

row = 5*['']
row[2] = '\multicolumn{' + '%d' % len(NLIST) + '}{c}' + \
            '{$\hat{\\theta}{}^{\star}_{\\text{lb}}$}'
row[4] = '\multicolumn{' + '%d' % len(NLIST) + '}{c}' + \
            '{$\hat{\\theta}{}^{\star}_{\\text{ub}}$}'
writeRow(fout, row)
insertCMidRule(fout, 3, 3 + len(NLIST) - 1, SPEC='lr')
insertCMidRule(fout, 3 + 1 + len(NLIST), 3 + 1 + 2*len(NLIST) - 1, SPEC='l')

row = (len(colspec)-1)*['']
row[0] = '\multicolumn{2}{r}{sample size}'
count = 0
for n in NLIST:
    ntotal = '$%d$' % round(n*NBASE)
    row[1 + count] = ntotal
    row[1 + 1 + len(NLIST) + count] = ntotal
    count = count + 1
writeRow(fout, row)
insertMidRule(fout)

row = len(colspec)*['']
FirstFlag = True
for p in PARAMUNIVERSE.keys():
    if p in ParamNames:
        if FirstFlag:
            FirstFlag = False
        else:
            insertCMidRule(fout, 2, len(colspec), 'l')

        row[0] = '\multirow{6}{*}{' + PARAMUNIVERSE[p] + '}'
        row[1] = 'true'
        for c in range(0,len(NLIST)):
            row[2 + c] = formatNum(TrueLB.loc[p][0])
            row[2 + 1 + len(NLIST) + c] = formatNum(TrueUB.loc[p][0])
        writeRow(fout, row)

        row[0] = ''
        row[1] = 'mean'

        for i, n in enumerate(NDIRLIST):
            row[2 + i] = formatNum(LB[n][p].mean())
            row[2 + 1 + len(NLIST) + i] = formatNum(UB[n][p].mean())
        writeRow(fout, row)

        row[0] = ''
        row[1] = 'std'
        for i, n in enumerate(NDIRLIST):
            row[2 + i] = formatNum(LB[n][p].std())
            row[2 + 1 + len(NLIST) + i] = formatNum(UB[n][p].std())
        writeRow(fout, row)

        row[0] = ''
        row[1] = 'rmse'
        for i, n in enumerate(NDIRLIST):
            row[2 + i] = formatNum(\
                        np.sqrt(\
                              (LB[n][p].mean() - TrueLB.loc[p][0])**2 \
                            + LB[n][p].var()\
                        ))
            row[2 + 1 + len(NLIST) + i] = formatNum(\
                        np.sqrt(\
                              (UB[n][p].mean() - TrueUB.loc[p][0])**2 \
                            + UB[n][p].var()\
                        ))
        writeRow(fout, row)

        row[0] = ''
        row[1] = '5/95\%'
        for i, n in enumerate(NDIRLIST):
            row[2 + i] = formatNum(percentile(LB[n][p], 5))
            row[2 + 1 + len(NLIST) + i] = formatNum(percentile(UB[n][p], 95))
        writeRow(fout, row)

        row[0] = ''
        row[1] = 'min/max'
        for i, n in enumerate(NDIRLIST):
            row[2 + i] = formatNum(min(LB[n][p]))
            row[2 + 1 + len(NLIST) + i] = formatNum(max(UB[n][p]))
        writeRow(fout, row)


insertMidRule(fout)
row = (len(colspec) - 1)*['']
row[0] = '\multicolumn{2}{r}{$\mathbb{P}[\Theta^{\star}=\emptyset \\text{ in sample}]$}'
for i, n in enumerate(NDIRLIST):
    row[1 + i] = formatNum(float(np.count_nonzero(MinCrit[n]))/len(MinCrit[n]))
    row[1 + 1 + len(NLIST) + i] = '--'
writeRow(fout, row)

insertBottomRule(fout)
endTable(fout)
createTableViewerAndCompile(os.path.join(CodeDir, FNVIEWTEMPLATE),
                            FNOUT, ResultsDir)

################################################################################
# Lets also make some plots while we're at it
################################################################################
for p in PARAMUNIVERSE.keys():
    if p in ParamNames:
        LBComb = np.vstack((LB[NDIRLIST[n]][p] for n in NPLOTLIST))
        LeftLB = np.amin(LBComb)
        LeftUB = np.amax(LBComb)
        UBComb = np.vstack((UB[NDIRLIST[n]][p] for n in NPLOTLIST))
        RightLB = np.amin(UBComb)
        RightUB = np.amax(UBComb)

        LeftGrid = np.linspace( LeftLB, LeftUB, PLOTGRIDDIM)
        RightGrid = np.linspace(RightLB, RightUB, PLOTGRIDDIM)

        fig,(axleft,axright) = plt.subplots(1,2, sharey=True)
        kwargs = dict(var_type='c', bw='normal_reference')

        Height = 0
        LegendPlots = []
        LegendLabels = []
        for i, n in enumerate(NPLOTLIST):
            LabelPDF = 'n = %d' % (NLIST[n]*NBASE)
            LegendLabels.append(LabelPDF)

            LeftDens = sm.nonparametric.KDEMultivariate(
                    data=LB[NDIRLIST[n]][p], **kwargs)
            LeftPlot = LeftDens.pdf(LeftGrid)
            pl, = axleft.plot(LeftGrid, LeftPlot, color=NPLOTCOLORLIST[i],
                    linewidth=PLOTLINEWIDTH)
            LegendPlots.append(pl)

            RightDens = sm.nonparametric.KDEMultivariate(
                    data=UB[NDIRLIST[n]][p], **kwargs)
            RightPlot = RightDens.pdf(RightGrid)
            pr, = axright.plot(RightGrid, RightPlot, color=NPLOTCOLORLIST[i],
                    linewidth=PLOTLINEWIDTH)

            Height = max(Height, np.amax(np.vstack([LeftPlot, RightPlot])))

        ### Cosmetic aspects of the plot
        axleft.set_xlim(LeftLB, LeftUB)
        axright.set_xlim(RightLB, RightUB)
        axleft.yaxis.set_ticks_position('none')
        axright.yaxis.set_ticks_position('none')
        axleft.xaxis.set_ticks_position('bottom')
        axright.xaxis.set_ticks_position('bottom')
        axleft.get_yaxis().set_ticks([])
        axright.get_yaxis().set_ticks([])

        YTop = (1 + PERCENTBUFFERY)*Height
        axes = plt.gca()
        axes.set_ylim([0, YTop])
        PlotLeftID = axleft.plot(   TrueLB.loc[p][0]*np.ones(PLOTGRIDDIM), \
                                    np.linspace(0, YTop, PLOTGRIDDIM))
        PlotRightID = axright.plot(  TrueUB.loc[p][0]*np.ones(PLOTGRIDDIM), \
                                    np.linspace(0, YTop, PLOTGRIDDIM))
        plt.setp(   [PlotLeftID, PlotRightID],
                    color='gray',
                    linewidth=1.5,
                    linestyle='--')
        xticks = [('%4.3f' % s) for s in \
                    [   LeftLB, LeftUB,
                        RightLB, RightUB,
                        TrueLB.loc[p][0], TrueUB.loc[p][0]]]
        axleft.xaxis.set_ticks([LeftLB, LeftUB, TrueLB.loc[p][0]])
        axright.xaxis.set_ticks([RightLB, RightUB, TrueUB.loc[p][0]])
        majorFormatter = FuncFormatter(removeLeadingZero)
        axleft.xaxis.set_major_formatter(majorFormatter)
        axright.xaxis.set_major_formatter(majorFormatter)

        axleft.set_title('Lower bound')
        axright.set_title('Upper bound')
        fig.legend( LegendPlots, LegendLabels,
                    loc='lower center',
                    ncol=2,
                    frameon=False)
        plt.subplots_adjust(bottom=0.15)

        plt.savefig(os.path.join(ResultsDir, 'MCDensityPlot_' + p))
