#!/usr/bin/env python
#coding=utf-8

import sys
import os
import numpy as np
import scipy.stats
import statsmodels.api as sm

from TableTools import *
from StatedepTools import *

################################################################################
# HARD-CODING
################################################################################
PANELTWOFORMAT = '%.2f'
NODATASTR = '--'
################################################################################
################################################################################

def main(DataFile, DestDir):
    CheckInputs(DataFile, DestDir)

    data = np.loadtxt(DataFile)
    fnout = os.path.join(\
                DestDir,\
                os.path.splitext(os.path.basename(DataFile))[0] + '.tex'\
            )
    fout = open(fnout, 'w')

# Table specification and header
    colspec = 'r' + (data.shape[1] + 1)*'c'
    startTable(fout, colspec)
    insertTopRule(fout)
    row = 3*['']
    row[1] = '\multicolumn{' + str(data.shape[1]) + '}{c}{time period $t$}'
    writeRow(fout, row)

    row = len(colspec)*['']
    for t in range(0,data.shape[1]):
        row[t+1] = '%d' % t
    writeRow(fout, row)

    insertCMidRule(fout, 2, len(colspec) - 1, 'lr')

# Probability of employment
    row2 = len(row)*['']
    emp_m = data.mean(axis=0);
    emp_se = scipy.stats.sem(data, axis=0);
    row[0] = '$\mathbb{P}[Y_{it} = 1]$'
    row2[0] = ''
    for t in range(0,data.shape[1]):
        row[t+1] = formatNum(emp_m[t])
        row2[t+1] = '(' + formatNum(emp_se[t], '\scriptsize') + ')'
    writeRow(fout, row)
    writeRow(fout, row2)

# Probability of transition
    transition = np.zeros((data.shape[0],data.shape[1] - 1));
    for t in range(1,data.shape[1]):
        for i in range(0,data.shape[0]):
            if data[i,t] != data[i,t-1]:
                transition[i,t-1] = 1;
            else:
                transition[i,t-1] = 0;

    trans_m = transition.mean(axis=0);
    trans_se = scipy.stats.sem(transition, axis=0);

    row[0] = '$\mathbb{P}[Y_{it} \\neq Y_{i(t-1)}]$'
    row[1] = '\multirow{2}{*}{' + NODATASTR + '}'
    row2[0] = ''
    row2[1] = ''
    for t in range(1,data.shape[1]):
        row[t+1] = formatNum(trans_m[t-1])
        row2[t+1] = '(' + formatNum(trans_se[t-1], '\scriptsize') + ')'
    writeRow(fout, row)
    writeRow(fout, row2)

# Probability of remaining unemployed
    p0g0 = np.zeros((data.shape[1]-1,2));
    tests = np.array([[1, 0]]);
    for t in range(1,data.shape[1]):
        Y = np.subtract(np.ones_like(data[:,t]), data[:,t]);
        X = data[:,t-1];
        X = sm.add_constant(X);
        model = sm.OLS(Y,X);
        results = model.fit(cov_type='HC0');
        results_lincom = results.t_test(tests);
        p0g0[t-1,0] = results_lincom.effect[0]
        p0g0[t-1,1] = results_lincom.sd[0]

    row[0] = '$\mathbb{P}[Y_{it} = 0 \\vert Y_{i(t-1)} = 0]$'
    for t in range(1,data.shape[1]):
        row[t+1] = formatNum(p0g0[t-1, 0])
        row2[t+1] = '(' + formatNum(p0g0[t-1, 1], '\scriptsize') + ')'
    writeRow(fout, row)
    writeRow(fout, row2)

# Probability of remaining employed
    p1g1 = np.zeros((data.shape[1]-1,2));
    naiveate = np.zeros((data.shape[1]-1,2));
    tests = np.array([[1, 1], [0, 1]]);
    for t in range(1,data.shape[1]):
        Y = data[:,t];
        X = data[:,t-1];
        X = sm.add_constant(X);
        model = sm.OLS(Y,X);
        results = model.fit(cov_type='HC0');
        results_lincom = results.t_test(tests);
        p1g1[t-1,0] = results_lincom.effect[0]
        p1g1[t-1,1] = results_lincom.sd[0]
        naiveate[t-1,0] = results_lincom.effect[1]
        naiveate[t-1,1] = results_lincom.sd[1]

    row[0] = '$\mathbb{P}[Y_{it} = 1 \\vert Y_{i(t-1)} = 1]$'
    for t in range(1,data.shape[1]):
        row[t+1] = formatNum(p1g1[t-1, 0])
        row2[t+1] = '(' + formatNum(p1g1[t-1, 1], '\scriptsize') + ')'
    writeRow(fout, row)
    writeRow(fout, row2)

# Naive ATE (already computed above)
    row[0] = 'naive ATE'
    for t in range(1,data.shape[1]):
        row[t+1] = formatNum(naiveate[t-1, 0])
        row2[t+1] = '(' + formatNum(naiveate[t-1, 1], '\scriptsize') + ')'
    writeRow(fout, row)
    writeRow(fout, row2)

################################################################################
# Second panel
################################################################################
# Header
    insertMidRule(fout)
    row = 3*['']
    row[1] = '\multicolumn{' + \
             str(data.shape[1]) + \
             '}{c}{percentage of agents with \ldots}'
    writeRow(fout, row)

    row = len(colspec)*['']
    for t in range(0,data.shape[1] + 1):
        row[t+1] = '%d' % t
    writeRow(fout, row)

    insertCMidRule(fout, 2, len(colspec), 'lr')

# Periods of unemployment
    totalemp = data.sum(axis=1)
    empcount = np.zeros(data.shape[1] + 1)
    for t in range(0, data.shape[1] + 1):
        a = t*np.ones(totalemp.shape);
        empcount[t] = \
            data.shape[0] - np.count_nonzero(np.subtract(totalemp, a))
    empcount = empcount[::-1] # Flip to be periods of unemployment

    row[0] = 'periods of unemployment'
    for t in range(0,data.shape[1] + 1):
        row[t+1] = PANELTWOFORMAT % (100*empcount[t]/data.shape[0])
    writeRow(fout, row)

# Unemployment spells
    spells = np.zeros(data.shape[0])
    for i in range(0,data.shape[0]):
        for t in range(0, data.shape[1]):
            if t == 0:
                if data[i,t] == 0:
                    spells[i] = spells[i] + 1
            elif t > 0:
                if data[i,t-1] == 1 and data[i,t] == 0:
                    spells[i] = spells[i] + 1
    spellcount = np.zeros(int(np.floor(data.shape[1]/2) + 2))
    for s in range(0,spellcount.shape[0]):
        a = s*np.ones(spells.shape)
        spellcount[s] = \
            data.shape[0] - np.count_nonzero(np.subtract(spells, a))

    row[0] = 'unemployment spells'
    for t in range(0, data.shape[1] + 1):
        if t < len(spellcount):
            row[t+1] = PANELTWOFORMAT % (100*spellcount[t]/data.shape[0])
        else:
            row[t+1] = NODATASTR
    writeRow(fout, row)

# Transitions
    totaltrans = transition.sum(axis=1);
    transcount = numpy.zeros(data.shape);
    for r in range(0,data.shape[1]):
        a = r*numpy.ones(totaltrans.shape);
        transcount[:,r] = (totaltrans == a);
    tottrans_m = transcount.mean(axis=0);
    tottrans_se = scipy.stats.sem(transcount, axis=0);

    row[0] = 'transitions'
    for t in range(0, data.shape[1]):
        row[t+1] = PANELTWOFORMAT % (100*tottrans_m[t])
    row[-1] = NODATASTR
    writeRow(fout, row)

# Build tex file
    insertBottomRule(fout)
    endTable(fout)
    CodeDir = os.path.dirname(os.path.abspath(sys.argv[0]))
    createTableViewerAndCompile(os.path.join(CodeDir, FNVIEWTEMPLATE),
                                os.path.basename(fnout),
                                DestDir)

if __name__ == '__main__':
    DataFile = os.path.abspath(sys.argv[1])
    DestDir = os.path.abspath(sys.argv[2])

    main(DataFile, DestDir)
