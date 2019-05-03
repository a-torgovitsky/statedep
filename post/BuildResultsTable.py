#!/usr/bin/env python
#coding=utf-8

import sys
import os
import math
import pandas as pd
import numpy
import re
import itertools
import copy
import csv
import jinja2
import operator
from collections import OrderedDict
from shutil import copyfile

from TableTools import *
from StatedepTools import *

################################################################################
# HARDCODED VARIABLES
################################################################################
FNOUT = "TableResults.tex"
CRFONTSIZE = '\\tiny'
ASSUMPTIONLABEL = OrderedDict(( \
    ('mST', '$\\text{ST}(m)$'),\
    ('SigmaST', '$\;\;\sigma$'),\
    ('qMTS', '$\\text{MTS}(q)$'),\
    ('MTR', 'MTR'),\
    ('MATR', 'MATR'),\
    ('TIV', 'TIV'),\
    ('DSC', 'DSC'),\
    ('MIV', 'MIV'),\
))
def TABLECHECKMARK(x):
    if x:
        return '\\checkmark'
    else:
        return ' '

def TABLEINTEGER(x):
    if (x >= 0):
        return  ('%d' % x)
    else:
        return ' '

def TABLEFLOAT(x):
    if (x >= 0):
        return ('%.2f' % x)
    else:
        return ' '

ASSUMPTIONFORMAT = dict(( \
    ('mST', TABLEINTEGER),\
    ('SigmaST', TABLEFLOAT),\
    ('qMTS', TABLEINTEGER),\
    ('TIV', TABLECHECKMARK),\
    ('MTR', TABLECHECKMARK),\
    ('MATR', TABLECHECKMARK),\
    ('DSC', TABLECHECKMARK),\
    ('MIV', TABLECHECKMARK),\
))

################################################################################
CodeDir = os.path.dirname(os.path.abspath(sys.argv[0]))
ResultsDir = os.path.abspath(sys.argv[1])
if not os.path.isdir(ResultsDir):
    print ('Could not find directory ' + ResultsDir)

SimNames = getSimNames(ResultsDir)
copyfile(os.path.join(CodeDir, FNVIEWTEMPLATE),
         os.path.join(ResultsDir, FNVIEWTEMPLATE))

fout = open(os.path.join(ResultsDir, FNOUT), 'w')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Column spec
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Create bounds dataframe
FileListBounds = [os.path.join(ResultsDir, SimNames[d], FNBOUNDS) \
    for d in range(0,len(SimNames))]
(lb, ub) = createBoundsDataFrame(FileListBounds)

FlagCR = False
try:
    FileListCRBounds = [os.path.join(ResultsDir, SimNames[d], FNCRBOUNDS) \
        for d in range(0,len(SimNames))]
    (crlb, crub) = createBoundsDataFrame(FileListCRBounds)
    FlagCR = True
except:
    crlb = lb.copy(deep=True)
    crub = ub.copy(deep=True)
    crlb[:] = float('nan')
    crub[:] = float('nan')

# Create constraints dataframe
FileListAssumptions = [os.path.join(ResultsDir, SimNames[d], \
    FNASSUMPTIONS) for d in range(0,len(SimNames))]
assert(len(FileListBounds) == len(FileListAssumptions))
assumptions = createDataFrame(FileListAssumptions)

# Left align for row labels
colspec = 'l' + len(lb.columns)*'c'
numcols = len(colspec)
startTable(fout, colspec)
insertTopRule(fout)

# For each simulation determine if its pdbr or not
ListPDBR = [os.path.isfile(os.path.join(ResultsDir, SimNames[d], FNPDBR)) \
        for d in range(0,len(SimNames))]
# Make sure PDBR's come at the end
for i in range(1,len(ListPDBR)):
    assert(not(ListPDBR[i-1] and not ListPDBR[i]))

################################################################################
# Column header
################################################################################
row = 2*['']
row[0] = '\t'
# sim types
row[1] = '\multicolumn{%d}{c}{\\textbf{DPO}}' % (len(ListPDBR) - sum(ListPDBR))
if sum(ListPDBR) > 0:
    row.append('\multicolumn{%d}{c}{\\textbf{PDBR}}' % sum(ListPDBR))
writeRow(fout, row)

# dividing lines
insertCMidRule(fout, 2, (len(ListPDBR) - sum(ListPDBR)) + 1, 'lr')
if sum(ListPDBR) > 0:
    insertCMidRule(fout, (len(ListPDBR) - sum(ListPDBR)) + 2, numcols, 'l')

# numbers
row = numcols*['']
for i in range(0, len(lb.columns)):
    row[1+i] = '(' + lb.columns[i].lstrip('0') + ')'
    row[1+i] = '\\textbf{' + row[1+i] + '}'

writeRow(fout, row)


################################################################################
# Assumptions
################################################################################
insertMidRule(fout)
subtitlerow = generateSubTitleRow(numcols, 'Assumptions')
writeRow(fout, subtitlerow)
insertMidRule(fout)

for a in ASSUMPTIONLABEL.keys():
    if a in assumptions.index and any(assumptions.loc[a,:].values > 0):
        row[0] = ASSUMPTIONLABEL[a]
        datastring = assumptions.loc[a,:].values
        f = ASSUMPTIONFORMAT.get(a)
        for i in range(0,len(datastring)):
            row[1 + i] = f(datastring[i])
        writeRow(fout, row)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Specification
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
insertMidRule(fout)
subtitlerow = generateSubTitleRow(numcols, 'Misspecification')
writeRow(fout, subtitlerow)
insertMidRule(fout)

FileListMinCriterion = [os.path.join(ResultsDir, SimNames[d], \
    FNMINCRITERION) for d in range(0,len(SimNames))]
assert(len(FileListBounds) == len(FileListMinCriterion))
mincrit = createDataFrame(FileListMinCriterion)
toprow = ['$\Theta^{\star} = \emptyset$']
for (pos, i) in enumerate(mincrit.iloc[0,:].values):
    s = '\multirow{2}{*}{'
    if ListPDBR[pos]:
        s = s + ''
    else:
        if i <= 0:
            s = s + 'No'
        else:
            s = s + 'Yes'
    s = s + '}'
    toprow.extend([s])
writeRow(fout, toprow)
bottomrow = ['in sample']
bottomrow.extend(['' for i in mincrit.iloc[0,:].values])
writeRow(fout, bottomrow, SKIPPT=5)

try:
    FileListMisspecification = [os.path.join(ResultsDir, SimNames[d], \
        FNMISSPECIFICATION) for d in range(0,len(SimNames))]
    assert(len(FileListBounds) == len(FileListMisspecification))
    misspec = createDataFrame(FileListMisspecification)
except:
    misspec = mincrit.copy(deep=True)
    misspec[:] = float('nan')

toprow = ['p-value for']
for (pos, i) in enumerate(misspec.iloc[0,:].values):
    if i < 1 and not ListPDBR[pos]:
        s = formatNum(i)
    else:
        s = ''
    toprow.extend(['\multirow{2}{*}{' + s + '}'])
writeRow(fout, toprow)
bottomrow = ['$H_{0}: \Theta^{\star} \\neq \emptyset$']
bottomrow.extend(['' for i in misspec.iloc[0,:].values])
writeRow(fout, bottomrow)

################################################################################
# Bounds
################################################################################
insertMidRule(fout)
if FlagCR:
    s = 'Bounds and 95\% Confidence Intervals'
else:
    s = 'Bounds'
subtitlerow = generateSubTitleRow(numcols, s)
writeRow(fout, subtitlerow)
insertMidRule(fout)

first = True

for p in PARAMUNIVERSE.keys():
    if p in lb.index:
        if first:
            first = False
        else:
            insertCMidRule(fout, 2, numcols, 'l')

        rowcrlb = numcols*['']
        rowlb = numcols*['']
        rowub = numcols*['']
        rowcrub = numcols*['']

        if FlagCR:
            rowcrlb[0] = '\t\multirow{4}{*}{' + PARAMUNIVERSE[p] + '}'
        else:
            rowlb[0] = '\t\multirow{2}{*}{' + PARAMUNIVERSE[p] + '}'

        crlbp = (crlb.loc[p,:].values).astype(numpy.float)
        crubp = (crub.loc[p,:].values).astype(numpy.float)
        lbp = (lb.loc[p,:].values).astype(numpy.float)
        ubp = (ub.loc[p,:].values).astype(numpy.float)

        pointid = [     (numpy.isfinite(lbp[i])) \
                    and (ubp[i] - lbp[i] == 0) \
                    for i in range(0,len(lbp))]

        for i in range(0,len(lbp)):
            if not pointid[i]:
                rowcrlb[1+i] = formatNum(crlbp[i], CRFONTSIZE)
                rowlb[1+i] = formatNum(lbp[i])
                rowub[1+i] = formatNum(ubp[i])
                rowcrub[1+i] = formatNum(crubp[i], CRFONTSIZE)
            elif pointid[i]:
                rowcrlb[1+i] = formatNum(crlbp[i], CRFONTSIZE)
                rowlb[1+i] = '\multirow{2}{*}{' + formatNum(lbp[i]) + '}'
                rowcrub[1+i] = formatNum(crubp[i], CRFONTSIZE)
            else:
                rowlb[1+i] = '---'

        if FlagCR:
            writeRow(fout, rowcrlb)
        writeRow(fout, rowlb)
        writeRow(fout, rowub)
        if FlagCR:
            writeRow(fout, rowcrub)

insertBottomRule(fout)
endTable(fout)
fout.close()
createTableViewerAndCompile(os.path.join(CodeDir, FNVIEWTEMPLATE),
                            FNOUT, ResultsDir)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Compile plots if multiple sigma estimates in this set
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(assumptions.loc['SigmaST',:].values > 0):
    order = assumptions.columns.argsort()

    with open(os.path.join(ResultsDir, 'BoundsSigma.csv'), mode='w') as outf:
        outfwriter = csv.writer(outf, delimiter=',')

        headerrow = ['Sigma']
        for p in lb.index:
            headerrow.extend([p + ' LB', p + ' UB'])
        outfwriter.writerow(headerrow)

        # for r in lb.columns:
        for r in order:
            datarow = []
            datarow.append(assumptions.ix['SigmaST', r])
            for p in lb.index:
                datarow.extend([lb.ix[p,r], ub.ix[p,r]])
            outfwriter.writerow(datarow)

    # reader = csv.reader(open(os.path.join(ResultsDir, 'BoundsSigma.csv')),
                        # delimiter=',')
    # sortedreader = sorted(reader, key=operator.itemgetter(0))

    with open(os.path.join(ResultsDir, 'CIsSigma.csv'), mode='w') as outf:
        outfwriter = csv.writer(outf, delimiter=',')

        headerrow = ['Sigma']
        for p in crlb.index:
            headerrow.extend([p + ' LB', p + ' UB'])
        outfwriter.writerow(headerrow)

        for r in order:
            datarow = []
            datarow.append(assumptions.ix['SigmaST', r])
            for p in crlb.index:
                datarow.extend([crlb.ix[p,r], crub.ix[p,r]])
            outfwriter.writerow(datarow)

    # Initialize Jinja templating
    latex_jinja_env = jinja2.Environment(
            block_start_string = '\BLOCK{',
            block_end_string = '}',
            variable_start_string = '\VAR{',
            variable_end_string = '}',
            comment_start_string = '\#{',
            comment_end_string = '}',
            trim_blocks = True,
            autoescape = False,
            loader = jinja2.FileSystemLoader(os.path.abspath('.'))
    )
    templatesigma = latex_jinja_env.get_template(FNTEMPLATESIGMA)

    idxcount = 1
    c = {}
    for p in lb.index:
        c['collb'] = idxcount
        c['colub'] = idxcount + 1
        idxcount = idxcount + 2

        c['parameter'] = PARAMUNIVERSE[p]

        fn = 'SigmaPlot' + p + '.tex'
        with open(os.path.join(ResultsDir, fn), 'w') as f:
            f.write(templatesigma.render(c))

        cwd = os.getcwd()
        os.chdir(ResultsDir)
        callLatexQuietly(fn)
        os.chdir(cwd)
