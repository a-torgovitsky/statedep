#!/usr/bin/env python
#coding=utf-8

from collections import OrderedDict
import sys
import os

################################################################################
# HARDCODING -- Filenames
################################################################################
FNLB = 'EstimatedLB.out'
FNUB = 'EstimatedUB.out'
FNTRUEBOUNDS = 'TrueBounds.out'
FNVIEWTEMPLATE = 'ViewTableTemplate.tex'
FNREJECTMASK = 'Reject_A%d_%s.out'
FNTESTPOINTS = 'TestPoints.out'
FNTRUEBOUNDS = 'TrueBounds.out'
FNCRBOUNDS = 'ConfidenceRegions_A5_CNS.out'
FNBOUNDS = 'Bounds.out'
FNASSUMPTIONS = 'Assumptions.out'
FNMISSPECIFICATION = 'Misspecification.out'
FNMINCRITERION = 'MinCriterion.out'
FNTEMPLATESIGMA = 'TemplateSigma.tex'
FNPDBR = 'PDBR.out'
MCTESTINGTMAX = 3

################################################################################
# HARDCODING -- Parameter dictionary with Latex representation
################################################################################
PARAMUNIVERSE = OrderedDict(( \
    ('ATE', '$\\text{ATE}_{t}$'), \
    ('ATE_Avg', '$\overline{\\text{ATE}}$'), \
    ('PSD', '$\\text{SD}_{\\text{avg}}^{+}$'), \
    ('PSD_Avg', '$\overline{\\text{SD}}^{+}$'), \
    ('PSD_G0', '$\\text{SD}_{\\text{avg}}^{+}(\cdot \\vert 0)$'), \
    ('PSD_G0_Avg', '$\overline{\\text{SD}}^{+}(\cdot \\vert 0)$'), \
    ('PSD_G00', '$\\text{SD}_{\\text{avg}}^{+}(\cdot \\vert 00)$'), \
    ('PSD_G00_Avg', '$\overline{\\text{SD}}^{+}(\cdot \\vert 00)$'), \
    ('PSD_G000' , '$\\text{SD}_{t}^{+}(\cdot \\vert 000)$'), \
    ('PSD_G000_Avg' , '$\overline{\\text{SD}}^{+}(\cdot \\vert 000)$'), \
    ('PSD_G1' , '$\\text{SD}_{\\text{avg}}^{+}(\cdot \\vert 1)$'), \
    ('PSD_G1_Avg' , '$\overline{\\text{SD}}^{+}(\cdot \\vert 1)$'), \
    ('PSD_G11' , '$\\text{SD}_{\\text{avg}}^{+}(\cdot \\vert 11)$'), \
    ('PSD_G11_Avg' , '$\overline{\\text{SD}}^{+}(\cdot \\vert 11)$'), \
    ('PSD_G111' , '$\\text{SD}_{t}^{+}(\cdot \\vert 111)$'), \
    ('PSD_G111_Avg' , '$\overline{\\text{SD}}^{+}(\cdot \\vert 111)$'), \
    ('NSD' , '$\\text{SD}_{t}^{-}$'), \
    ('NSD_Avg' , '$\overline{\\text{SD}}^{-}$'), \
    ('TSD' , '$\\text{SD}_{\\text{avg}}$'), \
    ('TSD_Avg' , '$\overline{\\text{SD}}$') \
))

def CheckInputs(fn, dirname):
    if not os.path.isfile(fn):
        print ('Could not find data file ' + fn)
        sys.exit()
    if not os.path.isdir(dirname):
        print ('Could not find destination directory ' + \
                dirname + ', so making it.')
        try:
            os.mkdir(dirname)
        except OSError:
            print ("Creating %s returned an error." % dirname)
            sys.exit()
