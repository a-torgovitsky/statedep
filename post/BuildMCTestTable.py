#!/usr/bin/env python
#coding=utf-8

import sys
import os
import numpy as np

from TableTools import *
from StatedepTools import *

################################################################################
# HARD-CODING
################################################################################
LEVELS = [1, 5, 10]
TESTS = ['SS', 'CNS']
TOL = 1e-3
FNOUT = 'TableMCTest.tex'

################################################################################
################################################################################
################################################################################
################################################################################
CodeDir = os.path.dirname(os.path.abspath(sys.argv[0]))
ResultsDir = os.path.abspath(sys.argv[1])
if not os.path.isdir(ResultsDir):
    print ('Could not find directory ' + ResultsDir)

# Load data
RejectProb = {}
for l in LEVELS:
    for t in TESTS:
        fn = os.path.join(ResultsDir, FNREJECTMASK % (l,t))
        R = np.loadtxt(fn)
        RejectProb[l,t] = R.mean(axis=0)

TrueBounds = np.loadtxt(os.path.join(ResultsDir, FNTRUEBOUNDS), usecols=(1,2))
TestPoints = np.loadtxt(os.path.join(ResultsDir, FNTESTPOINTS))

# Table specification and header
fout = open(os.path.join(ResultsDir, FNOUT), 'w')
colspec = (2 + len(TestPoints))*'c'
startTable(fout, colspec)
insertTopRule(fout)

row = 3*['']
row[2] = '\multicolumn{' + '%d' % len(TestPoints) + '}{c}' + \
        '{rejection probability of' + \
        ' $H_{0}: t \\in' + \
        ' \Theta^{\star}$ for $t = \ldots$}'

writeRow(fout, row)

row = len(colspec)*['']
row[0] = 'level'
row[1] = 'test'
row[2:] = [formatNum(t) \
           if (t < TrueBounds[0] - TOL or t > TrueBounds[1] + TOL) else \
           ('\\fbox{' + formatNum(t) + '}') for t in TestPoints]
writeRow(fout, row)
insertMidRule(fout)

for l in LEVELS:
    firsttest = True
    for t in TESTS:
        if firsttest:
            row[0] = ('\multirow{%d' % len(TESTS)) + \
                     '}{*}{' + ('%.2f' % (.01*l)).lstrip('0') + '}'
            firsttest = False
        else:
            row[0] = ''
        row[1] = t
        row[2:] = [formatNum(p) for p in RejectProb[l,t]]
        writeRow(fout, row)

insertBottomRule(fout)
endTable(fout)

createTableViewerAndCompile(os.path.join(CodeDir, FNVIEWTEMPLATE),
                            FNOUT, ResultsDir)
