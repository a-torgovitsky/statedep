#!/usr/bin/env python
#coding=utf-8

import os
import math
import numpy
from numpy import *
import pandas as pd
import itertools
import re
import subprocess
import shutil
import jinja2

def formatNum(x, fontsize=''):
    SMALL = .001

    if (math.isnan(x)):
        s = "--"
    if (x < 0):
        s = "-" + ("%4.3f" % x).lstrip('-0')
    if (x > SMALL and x < 1 - SMALL):
        s = ("%4.3f" % x).lstrip('0')
    if (x < SMALL):
        s = ".000"
    if (x > 1 - SMALL and x <= 1):
        s = "1.00"
    if (x > 1):
        s = ("%3.2f" % x)

    if fontsize:
        s = '{' + fontsize + ' ' + s + '}'

    return s

def removeLeadingZero(x, pos):
    return ('%4.3f' % x).lstrip('0')

def createBoundsDataFrame(FileListBounds):
    for count, f in enumerate(FileListBounds, start=0):
        # Pull off the bottom most directory name as the title
        bounds = pd.read_csv(f, sep=None, index_col=0, \
                    skipinitialspace=True, header=None, engine='python')
        # 10/12/18: Deprecated command -- but also not needed???
        # bounds.convert_objects(convert_numeric=True)
        if (count == 0):
            lb = pd.concat([bounds.iloc[:, [0]]], axis=1);
            ub = pd.concat([bounds.iloc[:, [1]]], axis=1);
        else:
            lb = pd.concat([lb, bounds.iloc[:, [0]]], axis=1)
            ub = pd.concat([ub, bounds.iloc[:, [1]]], axis=1)
        simid = os.path.basename(os.path.dirname(f))
        lb.rename(columns={1: simid}, inplace=True);
        ub.rename(columns={2: simid}, inplace=True);
    return (lb, ub)

def createDataFrame(FileList ):
    for count, f in enumerate(FileList, start=0):
        df_in = pd.read_csv(f, sep=None, index_col=0, \
                    skipinitialspace=True, header=None, engine='python')
        if (count == 0):
            df = df_in;
        else:
            df = pd.concat([df, df_in.iloc[:,[0]]], axis=1)
        simid = os.path.basename(os.path.dirname(f))
        df.rename(columns={1: simid}, inplace=True);
    return df

def getSimNames(ResultsDir):
    # Note that d[1] is all subdirectories in d
    # The simulation folders don't have any subdirectories
    # So if not d[1] evaluates to True for these folders
    # Need this conditioning because os.walk includes the top folder
    # which I don't want to include.
    SimDirs = [d[0] for d in os.walk(ResultsDir) if not d[1]]
    # Each directory is a column
    SimNames = [os.path.basename(d) for d in SimDirs]
    SimNames.sort()
    return SimNames

def startTable(outfile, array):
    tablespec = '\\begin{tabular}{@{}'
    for i in array:
        assert(i in {'l', 'c', 'r'})
        tablespec = tablespec + i
    tablespec = tablespec + '@{}}\n'
    outfile.write(tablespec);

def writeRow(fout, inrow, SKIPPT=0):
    l = len(inrow)
    outrow = '\t'
    for i in range(0,l):
        outrow = outrow + inrow[i]
        if i < l-1:
            outrow = outrow + ' & '
        else:
            outrow = outrow + ' \\\\'
            if SKIPPT:
                outrow = outrow + '[%dpt]' % SKIPPT
            outrow = outrow + '\n'
    fout.write(outrow)
    return outrow

def insertTopRule(outfile):
    outfile.write('\t\\toprule\n')

def insertBottomRule(outfile):
    outfile.write('\t\\bottomrule\n')

def insertMidRule(outfile):
    outfile.write('\t\\midrule\n')

def insertCMidRule(outfile, start, stop, SPEC=''):
    s = '\t\cmidrule'
    if SPEC:
        s = s + '(' + SPEC + ')'
    s = s + '{' + str(start) + '-' + str(stop) + '}\n'
    outfile.write(s)

def endTable(outfile):
    tableend = '\end{tabular}';
    outfile.write(tableend);
    outfile.close();

def callLatexQuietly(texfile, outputdir=''):
    subprocess.call('pdflatex -halt-on-error ' + texfile \
        + ' | grep -a3 ^!', shell=True)
    if outputdir:
        shutil.copy(os.path.splitext(texfile)[0] + '.pdf', outputdir)
    s = subprocess.Popen(['latexmk', '-c'], \
      stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]

def createTableViewerAndCompile(fnviewtemplate, fntable, destination):
    # Initialize Jinja templating
    templatedir = os.path.dirname(os.path.abspath(fnviewtemplate))
    latex_jinja_env = jinja2.Environment(
            block_start_string = '\BLOCK{',
            block_end_string = '}',
            variable_start_string = '\VAR{',
            variable_end_string = '}',
            comment_start_string = '\#{',
            comment_end_string = '}',
            trim_blocks = True,
            autoescape = False,
            loader = jinja2.FileSystemLoader(templatedir)
    )
    templateview = latex_jinja_env.get_template(
            os.path.basename(fnviewtemplate))
    c = {}
    c['tablefn'] = fntable
    fnview = 'View-' + fntable

    with open(os.path.join(destination, fnview), 'w') as f:
        f.write(templateview.render(c))

    cwd = os.getcwd()
    os.chdir(destination)
    callLatexQuietly(fnview, '')
    os.chdir(cwd)

def generateSubTitleRow(numcols, text):
    subtitlerow = 2*['']
    subtitlerow[0] = '\t'
    subtitlerow[1] = '\multicolumn{' + \
                     ('%d' % (numcols-1)) + \
                     '}{c}{\\textbf{' + text + '}}'
    return subtitlerow
