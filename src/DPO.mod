################################################################################
# DPO.mod
#
# A set/param is in ALL CAPS if and only it's value is passed from MATLAB.
################################################################################

# Sample size variables---only used here when loading the data
param N integer, >= 1; # Sample size
param SQRTN = sqrt(N);

# Time is indexed as t=0,1,...,T
param T integer >= 1 default 4; # Overall length of data
set DomU = 1..(2^(1 + 2*T));
set YSEQS {t in 1..(T+1)}; # Set of sequences of length t

################################################################################
# OBSERVED DATA
################################################################################
# YHAT contains all Y sequences that have positive probability in the data
# It's desirable for the way I handle resampling to not automatically
#   have YHAT declared as {y : Q[y] > 0}.
# This is because I don't want YHAT to change when I redraw samples.
# So YHAT gets declared by MATLAB rather than automatically.
#
# Q[y] is the observed probability of y = (Y_{0},Y_{1},...,Y_{T})
set YHAT ordered by Integers within YSEQS[T+1];
param Q {y in YHAT} in [0,1] default 0;

# Saved quantities from the sample
param Q_Sample {y in YHAT} default Q[y];

################################################################################
# UNOBSERVED OBJECTS
################################################################################
################################################################################
# UHAT contains all sequences of u that could have positive probability given
# the observed data.
#
# Using this set (vs. 1..SU) helps on memory usage.
#
# P[u] is the unobserved probability of
#   U = (U_{0}, U_{1}(0),...,U_{T}(0),
#               U_{1}(1),...,U_{T}(1))
# This set is only indexed over UHAT for computational speed.
# H[u] is the local deviation version used in CNS
################################################################################
set UHAT within DomU;
var P {u in UHAT} in [0,1];
var H {u in UHAT} default 0;

param r_ZeroOne_CNS >= 0, default 0;
subject to ZeroOne_CNS {u in UHAT}:
    1 - r_ZeroOne_CNS >= P[u] + H[u]/SQRTN >= r_ZeroOne_CNS;

#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
# PARAMETERS TO OPTIMIZE
# AND THEIR MIN/MAX OBJECTIVE FUNCTIONS
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
# What time period to evaluate the t-specific parameters at
# Note that period T+1 indicates the average over the other periods
param OptPeriod in {1..T+1} default T+1;

################################################################################
# POSITIVE STATE DEPENDENCE
#
# U_PSD[t] is the set of u such that
# Pr[U_{t}(0) = 0, U_{t}(1) = 1] = sum_{u in U_PSD} P(u);
# This set is determined in MATLAB.
################################################################################
set U_PSD {t in 1..T} within DomU;
var PSD {t in 1..T+1} =
    if (t <= T) then (
        sum {u in (U_PSD[t] inter UHAT)} P[u])
    else (
        (1/T)*(sum {tt in 1..T} PSD[tt])
    );
var PSD_H {t in 1..T+1} =
    if (t <= T) then (
        sum {u in (U_PSD[t] inter UHAT)} H[u])
    else (
        (1/T)*(sum {tt in 1..T} PSD_H[tt])
    );
minimize minPSD : PSD[OptPeriod];
maximize maxPSD : PSD[OptPeriod];

# NSD = negative state dependence
set U_NSD {t in 1..T} within DomU;
var NSD {t in 1..T+1} =
    if (t <= T) then (
        sum {u in (U_NSD[t] inter UHAT)} P[u])
    else (
        (1/T)*(sum {tt in 1..T} NSD[tt])
    );
var NSD_H {t in 1..T+1} =
    if (t <= T) then (
        sum {u in (U_NSD[t] inter UHAT)} H[u]
    ) else (
        (1/T)*(sum {tt in 1..T} NSD_H[tt])
    );
minimize minNSD : NSD[OptPeriod];
maximize maxNSD : NSD[OptPeriod];

# ATE = average treatment effect
var ATE {t in 1..T+1} = PSD[t] - NSD[t];
minimize minATE : ATE[OptPeriod];
maximize maxATE : ATE[OptPeriod];

# TSD = total state dependence
var TSD {t in 1..T+1} = PSD[t] + NSD[t];
var TSD_H {t in 1..T+1} = PSD_H[t] + NSD_H[t];
minimize minTSD : TSD[OptPeriod];
maximize maxTSD : TSD[OptPeriod];

################################################################################
################################################################################
################################################################################
# CONDITIONAL PARAMETERS
################################################################################
################################################################################
################################################################################
################################################################################
# PSD_G0 is PSD given Y_{t} = 0
#
# Y_G0 is all Y sequences with Y_{t} = 0
# G0 is the probability that Y_{t} = 0
# PSD_G0 is used for computing identified sets directly and is linear
# since G0 is taken directly from the observed data and hence a parameter.
#   (This is equivalent as long as the identified set is non-empty.)
# The numerator and denominator variables are used for inference,
# as are their H counterparts (for CNS).
################################################################################
set Y_G0 {t in 1..T};
param g0 {t in 1..T} = sum {y in (Y_G0[t] inter YHAT)} Q[y];
param g0_Sample {t in 1..T} = sum {y in (Y_G0[t] inter YHAT)} Q_Sample[y];

set U_PSD_G0_NUM {t in 1..T} within DomU;
var PSD_G0_Num {t in 1..T} = sum {u in (U_PSD_G0_NUM[t] inter UHAT)} P[u];
var PSD_G0_Num_H {t in 1..T} = sum {u in (U_PSD_G0_NUM[t] inter UHAT)} H[u];

var PSD_G0 {t in 1..T+1} =
    if (t <= T) then (
        if (g0[t] > 0) then (
            PSD_G0_Num[t]/g0[t]
        ) else (
            Infinity
        )
    ) else (
        (1/T)*(sum {tt in 1..T : g0[tt] > 0} PSD_G0[tt])
    );
var PSD_G0_Sample {t in 1..T+1} =
    if (t <= T) then (
        if (g0_Sample[t] > 0) then (
            PSD_G0_Num[t]/g0_Sample[t]
        ) else (
            Infinity
        )
    ) else (
        (1/T)*(sum {tt in 1..T : g0_Sample[tt] > 0} PSD_G0_Sample[tt])
    );
var PSD_G0_Sample_H {t in 1..T+1} =
    if (t <= T) then (
        if (g0_Sample[t] > 0) then (
            PSD_G0_Num_H[t]/g0_Sample[t]
        ) else (
            Infinity
        )
    ) else (
        (1/T)*(sum {tt in 1..T : g0_Sample[tt] > 0} PSD_G0_Sample_H[tt])
    );

minimize minPSD_G0 : PSD_G0[OptPeriod];
maximize maxPSD_G0 : PSD_G0[OptPeriod];
subject to PSD_G0_Bound {t in 1..T+1}:
    (if (t > T) then 1 else
        (if (g0[t] > 0) then 1 else Infinity)
    ) >= PSD_G0[t] >= 0;

################################################################################
# PSD_G00 is PSD given Y_{t} = 0, Y_{t-1} = 0
################################################################################
set Y_G00 {t in 1..T};
param g00 {t in 1..T} = sum {y in (Y_G00[t] inter YHAT)} Q[y];
param g00_Sample {t in 1..T} = sum {y in (Y_G00[t] inter YHAT)} Q_Sample[y];

set U_PSD_G00_NUM {t in 1..T} within DomU;
var PSD_G00_Num {t in 1..T} = sum {u in (U_PSD_G00_NUM[t] inter UHAT)} P[u];
var PSD_G00_Num_H {t in 1..T} = sum {u in (U_PSD_G00_NUM[t] inter UHAT)} H[u];

var PSD_G00 {t in 1..T+1} =
    if (t <= T) then (
        if (g00[t] > 0) then (
            PSD_G00_Num[t]/g00[t]
        ) else (
            Infinity
        )
    ) else (
        (1/T)*(sum {tt in 1..T : g00[tt] > 0} PSD_G00[tt])
    );
var PSD_G00_Sample {t in 1..T+1} =
    if (t <= T) then (
        if (g00_Sample[t] > 0) then (
            PSD_G00_Num[t]/g00_Sample[t]
        ) else (
            Infinity
        )
    ) else (
        (1/T)*(sum {tt in 1..T : g00_Sample[tt] > 0} PSD_G00_Sample[tt])
    );
var PSD_G00_Sample_H {t in 1..T+1} =
    if (t <= T) then (
        if (g00_Sample[t] > 0) then (
            PSD_G00_Num_H[t]/g00_Sample[t]
        ) else (
            Infinity
        )
    ) else (
        (1/T)*(sum {tt in 1..T : g00_Sample[tt] > 0} PSD_G00_Sample_H[tt])
    );

minimize minPSD_G00 : PSD_G00[OptPeriod];
maximize maxPSD_G00 : PSD_G00[OptPeriod];
subject to PSD_G00_Bound {t in 1..T+1}:
    (if (t > T) then 1 else
        (if (g00[t] > 0) then 1 else Infinity)
    ) >= PSD_G00[t] >= 0;

################################################################################
# PSD_G1 is PSD given Y_{t} = 1
################################################################################
set Y_G1 {t in 1..T};
param g1 {t in 1..T} = sum {y in (Y_G1[t] inter YHAT)} Q[y];
param g1_Sample {t in 1..T} = sum {y in (Y_G1[t] inter YHAT)} Q_Sample[y];

set U_PSD_G1_NUM {t in 1..T} within DomU;
var PSD_G1_Num {t in 1..T} = sum {u in (U_PSD_G1_NUM[t] inter UHAT)} P[u];
var PSD_G1_Num_H {t in 1..T} = sum {u in (U_PSD_G1_NUM[t] inter UHAT)} H[u];

var PSD_G1 {t in 1..T+1} =
    if (t <= T) then (
        if (g1[t] > 0) then (
            PSD_G1_Num[t]/g1[t]
        ) else (
            Infinity
        )
    ) else (
        (1/T)*(sum {tt in 1..T : g1[tt] > 0} PSD_G1[tt])
    );
var PSD_G1_Sample {t in 1..T+1} =
    if (t <= T) then (
        if (g1_Sample[t] > 0) then (
            PSD_G1_Num[t]/g1_Sample[t]
        ) else (
            Infinity
        )
    ) else (
        (1/T)*(sum {tt in 1..T : g1_Sample[tt] > 0} PSD_G1_Sample[tt])
    );
var PSD_G1_Sample_H {t in 1..T+1} =
    if (t <= T) then (
        if (g1_Sample[t] > 0) then (
            PSD_G1_Num_H[t]/g1_Sample[t]
        ) else (
            Infinity
        )
    ) else (
        (1/T)*(sum {tt in 1..T : g1_Sample[tt] > 0} PSD_G1_Sample_H[tt])
    );

minimize minPSD_G1 : PSD_G1[OptPeriod];
maximize maxPSD_G1 : PSD_G1[OptPeriod];
subject to PSD_G1_Bound {t in 1..T+1}:
    (if (t > T) then 1 else
        (if (g1[t] > 0) then 1 else Infinity)
    ) >= PSD_G1[t] >= 0;

################################################################################
# PSD_G11 is PSD given Y_{t} = 1, Y_{t-1} = 1
################################################################################
set Y_G11 {t in 1..T};
param g11 {t in 1..T} = sum {y in (Y_G11[t] inter YHAT)} Q[y];
param g11_Sample {t in 1..T} = sum {y in (Y_G11[t] inter YHAT)} Q_Sample[y];

set U_PSD_G11_NUM {t in 1..T} within DomU;
var PSD_G11_Num {t in 1..T} = sum {u in (U_PSD_G11_NUM[t] inter UHAT)} P[u];
var PSD_G11_Num_H {t in 1..T} = sum {u in (U_PSD_G11_NUM[t] inter UHAT)} H[u];

var PSD_G11 {t in 1..T+1} =
    if (t <= T) then (
        if (g11[t] > 0) then (
            PSD_G11_Num[t]/g11[t]
        ) else (
            Infinity
        )
    ) else (
        (1/T)*(sum {tt in 1..T : g11[tt] > 0} PSD_G11[tt])
    );
var PSD_G11_Sample {t in 1..T+1} =
    if (t <= T) then (
        if (g11_Sample[t] > 0) then (
            PSD_G11_Num[t]/g11_Sample[t]
        ) else (
            Infinity
        )
    ) else (
        (1/T)*(sum {tt in 1..T : g11_Sample[tt] > 0} PSD_G11_Sample[tt])
    );
var PSD_G11_Sample_H {t in 1..T+1} =
    if (t <= T) then (
        if (g11_Sample[t] > 0) then (
            PSD_G11_Num_H[t]/g11_Sample[t]
        ) else (
            Infinity
        )
    ) else (
        (1/T)*(sum {tt in 1..T : g11_Sample[tt] > 0} PSD_G11_Sample_H[tt])
    );

minimize minPSD_G11 : PSD_G11[OptPeriod];
maximize maxPSD_G11 : PSD_G11[OptPeriod];
subject to PSD_G11_Bound {t in 1..T+1}:
    (if (t > T) then 1 else
        (if (g11[t] > 0) then 1 else Infinity)
    ) >= PSD_G11[t] >= 0;

################################################################################
################################################################################
################################################################################
# PARAMETER CONSTRAINTS IN INFERENCE PROCEDURES
#   Note that choosing "MS" does not actually impose any constraint, since
#   this sets variable FixGap = 0, so constraint FixParameter is 0 = 0.
#   This is a convenient way to re-use code for the misspecification test.
################################################################################
################################################################################
################################################################################
set InferenceParams = { 'TSD', 'PSD',
                        'PSD_G0', 'PSD_G00',
                        'PSD_G1', 'PSD_G11', 'MS'};
param ActiveParam symbolic within InferenceParams, default 'PSD';
param Fix default .5;

var FixGap {p in InferenceParams} =
    (if (p == 'TSD') then (
        TSD[OptPeriod] - Fix
    ) else 0)
    +
    (if (p == 'PSD') then (
        PSD[OptPeriod] - Fix
    ) else 0)
    +
    (if (p == 'PSD_G0') then (
        PSD_G0_Sample[OptPeriod] - Fix
    ) else 0)
    +
    (if (p == 'PSD_G00') then (
        PSD_G00_Sample[OptPeriod] - Fix
    ) else 0)
    +
    (if (p == 'PSD_G1') then (
        PSD_G1_Sample[OptPeriod] - Fix
    ) else 0)
    +
    (if (p == 'PSD_G11') then (
        PSD_G11_Sample[OptPeriod] - Fix
    ) else 0)
    +
    (if (p == 'MS') then (
        0
    ) else 0);
subject to FixParameter: FixGap[ActiveParam] = 0;

# Note that the form is different for PSD since it is separable in Fix
var FixGap_CNS {p in InferenceParams} =
    (if (p == 'TSD') then (
        TSD_H[OptPeriod]
    ) else 0)
    +
    (if (p == 'PSD') then (
        PSD_H[OptPeriod]
    ) else 0)
    +
    (if (p == 'PSD_G0') then (
        PSD_G0_Sample_H[OptPeriod]
    ) else 0)
    +
    (if (p == 'PSD_G00') then (
        PSD_G00_Sample_H[OptPeriod]
    ) else 0)
    +
    (if (p == 'PSD_G1') then (
        PSD_G1_Sample_H[OptPeriod]
    ) else 0)
    +
    (if (p == 'PSD_G11') then (
        PSD_G11_Sample_H[OptPeriod]
    ) else 0)
    +
    (if (p == 'MS') then (
        0
    ) else 0);
subject to FixParameter_CNS: FixGap_CNS[ActiveParam] = 0;

################################################################################
################################################################################
################################################################################
# REQUIRED CONSTRAINTS
################################################################################
################################################################################
################################################################################

################################################################################
# PROPER PMF
#
# Make sure P[ ] is a proper pmf for
# U = (U_{0}, U_{1}(0),...,U_{T}(0),
#             U_{1}(1),...,U_{T}(1)))
################################################################################
subject to Probability : sum {u in UHAT} P[u] = 1;
subject to Probability_CNS : sum {u in UHAT} H[u] = 0;

################################################################################
################################################################################
################################################################################
# OBSERVATIONAL EQUIVALENCE
#
# Defining this as a "gap" variable is useful to then define a criterion
# functions for estimating.
#
# ObsEq_Gap_Sample is useful for saving data, which is needed in CNS
################################################################################
################################################################################
################################################################################
set U_OEQ {y in YHAT} within UHAT;
var ObsEq_Pred {y in YHAT} = sum {u in U_OEQ[y]} P[u];
var ObsEq_Gap {y in YHAT} = Q[y] - ObsEq_Pred[y];
var ObsEq_Gap_Sample {y in YHAT} = Q_Sample[y] - ObsEq_Pred[y];
subject to ObsEq {y in YHAT}: ObsEq_Gap[y] = 0;

################################################################################
# The equality criterion is the scaled sum of
#
#   | ObsEq_Gap[y] |
#
# The absolute value is turned into a linear function by introducing an
# auxiliary variable ObsEq_Gap_Abs that is constrained to
# be greater than \pm ObsEq_Gap.
# Since the criterion is being minimized, and since ObsEq_Gap_Abs appears
# nowhere else in the program, the optimal solution will be to force
# ObsEq_Gap down until it hits the larger of \pm ObsEq_Gap,
# i.e. until it hits the absolute value of ObsEq_Gap.
################################################################################
var ObsEq_Gap_Abs {y in YHAT} >= 0;
subject to ObsEq_Gap_Abs_Pos {y in YHAT}: ObsEq_Gap_Abs[y] >= ObsEq_Gap[y];
subject to ObsEq_Gap_Abs_Neg {y in YHAT}: ObsEq_Gap_Abs[y] >= -1*ObsEq_Gap[y];
var EqCriterion = SQRTN*(sum {y in YHAT} ObsEq_Gap_Abs[y]);

################################################################################
# The "sample" equality criterion is used in the CNS constraint---
# see EstimatedIDSet below.
#
# The constraint is: EqCriterion_Sample <= (1+tau)CriterionHat,
# where CriterionHat and tau are parameters.
# EqCriterion_Sample is the scaled sum of
#
#   | ObsEq_Gap_Sample[y] |
#
# This gets turned into a linear constraint by introducing positive auxiliary
# variables ObsEq_Gap_Sample_Pos and ObsEq_Gap_Sample_Neg, replacing
#
# | ObsEq_Gap_Sample | by ObsEq_Gap_Sample_Pos + ObsEq_Gap_Sample_Neg
#
# and including the constraint that
#
# ObsEq_Gap_Sample = ObsEq_Gap_Sample_Pos - ObsEq_Gap_Sample_Neg.
#
# It can be shown that
#
#   | ObsEq_Gap_Sample[y] | <= c
#
# if and only if such auxiliary variables exist.
#
# Note that this differs from the reformulation strategy for the equality
# criterion above because this is a constraint (vs. an objective).
# Indeed, it is not in general the case that for a feasible solution
#
#   | ObsEq_Gap_Sample | = ObsEq_Gap_Sample_Pos + ObsEq_Gap_Sample_Neg
#
# since the statement is there EXIST such auxiliary variables, not that
# they are unique.
################################################################################
var ObsEq_Gap_Sample_Pos {y in YHAT} >= 0;
var ObsEq_Gap_Sample_Neg {y in YHAT} >= 0;
var EqCriterion_Sample =
    SQRTN*(sum {y in YHAT} (ObsEq_Gap_Sample_Pos[y] + ObsEq_Gap_Sample_Neg[y]));
subject to ObsEq_Gap_Sample_Abs {y in YHAT}:
    ObsEq_Gap_Sample[y] = ObsEq_Gap_Sample_Pos[y] - ObsEq_Gap_Sample_Neg[y];

################################################################################
# For the observational equivalence conditions, Vstar reduces to
#    (\Prob_{n}^{*}[Y = j] - \Prob_{n}[Y = j])
#   /(\Prob_{n}[Y = j](1-\Prob_{n}[Y=j]))
#
# EqCriterion_CNS is reformulated analogously to EqCriterion above.
################################################################################
param VstarObsEq {y in YHAT} = (Q[y] - Q_Sample[y]);
var ObsEq_Gap_CNS {y in YHAT} =
    VstarObsEq[y] - (1/SQRTN)*(sum {u in U_OEQ[y]} H[u]);
var ObsEq_Gap_CNS_Abs {y in YHAT} >= 0;
subject to ObsEq_Gap_CNS_Abs_Pos {y in YHAT}:
    ObsEq_Gap_CNS_Abs[y] >= ObsEq_Gap_CNS[y];
subject to ObsEq_Gap_CNS_Abs_Neg {y in YHAT}:
    ObsEq_Gap_CNS_Abs[y] >= -1*ObsEq_Gap_CNS[y];
var EqCriterion_CNS = SQRTN*(sum {y in YHAT} ObsEq_Gap_CNS_Abs[y]);

################################################################################
################################################################################
################################################################################
# CONSTRAINTS FOR ADDITIONAL ASSUMPTIONS
################################################################################
################################################################################
################################################################################

################################################################################
# Stationarity (ST)
#
# meaning that:
#       P[U_{t}^{m}(0) = u^{m}(0) and U_{t}^{m}(1) = u^{m}(1)]
#    =  P[U_{1}^{m}(0) = u^{m}(0) and U_{1}^{m}(1) = u^{m}(1)]
# for all t >= 2 and <= T - m.
#
# where
#   U_{t}^{m}(d) = (U_{t}(d),U_{t+1}(d),...,U_{t+m}(d)) for d = 0,1
#   m --> DIMST is an integer chosen by the analyst
#   0 <= DIMST <= (T-2) if T >= 2
#   This definition only makes sense for t <= T - m
#
# Note that a realization of U_{t}^{m} = (U_{t}^{m}(0), U_{t}^{m}(1))
#   has length LUST = 2*(DIMST + 1).
# So there are 2^LUST total realizations of each U_{t}^{m}
################################################################################
param DIMST integer >= 0, <= max(0,(T-2)), default 0;
set U_ST_EQUATE;
set U_ST {t in 1..(T-DIMST), u in U_ST_EQUATE} within DomU;

subject to ST {t in 2..(T-DIMST), u in U_ST_EQUATE}:
    sum {uu in (U_ST[t,u] inter UHAT)} P[uu]
    =
    sum {uu in (U_ST[t-1,u] inter UHAT)} P[uu];

subject to ST_CNS {t in 2..(T-DIMST), u in U_ST_EQUATE}:
    sum {uu in (U_ST[t,u] inter UHAT)} H[uu]
    =
    sum {uu in (U_ST[t-1,u] inter UHAT)} H[uu];

param SIGMAST {t in 1..(T-DIMST), tt in 1..(T-DIMST)} >= 0, default 0;

subject to STR {t in 1..(T-DIMST), tt in 1..(T-DIMST), u in U_ST_EQUATE}:
    sum {uu in (U_ST[t,u] inter UHAT)} P[uu]
    <=
    if (SIGMAST[t,tt] == Infinity) then
        Infinity
    else
        (1 + SIGMAST[t,tt])*(sum {uu in (U_ST[tt,u] inter UHAT)} P[uu]);

subject to STR_CNS {t in 1..(T-DIMST), tt in 1..(T-DIMST), u in U_ST_EQUATE}:
    sum {uu in (U_ST[t,u] inter UHAT)}(P[uu] + H[uu]/SQRTN)
    <=
    if (SIGMAST[t,tt] == Infinity) then
        Infinity
    else
        (1 + SIGMAST[t,tt])*
        sum {uu in (U_ST[tt,u] inter UHAT)} (P[uu] + H[uu]/SQRTN);

subject to STL {t in 1..(T-DIMST), tt in 1..(T-DIMST), u in U_ST_EQUATE}:
    sum {uu in (U_ST[t,u] inter UHAT)} P[uu]
    >=
    if (SIGMAST[t,tt] == Infinity) then
        -Infinity
    else
        (1 - SIGMAST[t,tt])*
        (sum {uu in (U_ST[tt,u] inter UHAT)} P[uu]);

subject to STL_CNS {t in 1..(T-DIMST), tt in 1..(T-DIMST), u in U_ST_EQUATE}:
    sum {uu in (U_ST[t,u] inter UHAT)} (P[uu] + H[uu]/SQRTN)
    >=
    if (SIGMAST[t,tt] == Infinity) then
        -Infinity
    else
        (1 - SIGMAST[t,tt])*
        sum {uu in (U_ST[tt,u] inter UHAT)} (P[uu] + H[uu]/SQRTN);

################################################################################
# DSC (diminishing serial correlation)
#
# P[U_{t}(d) = 1, U_{t+r}(d) = 1] is increasing in |r|.
#
# This was easiest to do as the following two conditions:
#
# P[U_{t}(d) = 1, U_{t+r}(d) = 1] >= P[U_{t}(d) = 1, U_{t+r+1}(d) = 1]
#   for t = 1,...,T-2
#       r > t, but with r+1 <= T
#   and d = 0, 1
#
# and
#
# P[U_{t}(d) = 1, U_{t-r}(d) = 1] >= P[U_{t}(d) = 1, U_{t-r-1}(d) = 1]
#   for t = 3,...,T
#       r < t, but with r-1 >= 1
#   and d = 0, 1
################################################################################
set U_DSC {s1 in 1..T, s2 in 1..T, d in 0..1 : s1 <> s2} within DomU;

subject to DSC_Forw {t in 1..(T-2), r in (t+1)..(T-1), d in 0..1}:
        (sum {u in (U_DSC[t,r,d] inter UHAT)} P[u])
        >=
        (sum {u in (U_DSC[t,r+1,d] inter UHAT)} P[u]);

subject to DSC_Back {t in 3..T, r in (t-1)..2 by -1, d in 0..1}:
        (sum {u in (U_DSC[t,r,d] inter UHAT)} P[u])
        >=
        (sum {u in (U_DSC[t,r-1,d] inter UHAT)} P[u]);

################################################################################
# MTR (monotone treatment response)
#
# i.e. P[U_{t}(0) = 1, U_{t}(1) = 0] = NSD[t] = 0 for all t,
################################################################################
subject to MTR {t in 1..T}: NSD[t] = 0;
subject to MTR_CNS {t in 1..T}: NSD_H[t] = 0;

################################################################################
# MATR (monotone average treatment response)
#
# i.e. P[U_{t}(1) = 1] >= P[U_{t}(0) = 1] for all t,
################################################################################
set U_AE0 {t in 1..T} within DomU;
var AE0 {t in 1..T} = sum {u in (U_AE0[t] inter UHAT)} P[u];

set U_AE1 {t in 1..T} within DomU;
var AE1 {t in 1..T} = sum {u in (U_AE1[t] inter UHAT)} P[u];

subject to MATR {t in 1..T}: AE1[t] >= AE0[t];

################################################################################
# MTS (monotone treatment selection)
#
# The condition is
# P[U_{t}(d) = 1 | Y_{t-1} = 1, Y_{t-2} = y_{2}, ..., Y_{t-q} = y_{q}]
#   >=
# P[U_{t}(d) = 1 | Y_{t-1} = 0, Y_{t-2} = y_{2}, ..., Y_{t-q} = y_{q}]
#   for t = 2..T, d = 0,1 and all (y_{2},...,y_{q}) sequences
# where q is "DimMTS" in the Matlab code.
#
# Y_MTS_DENOM_EQUATE[t,q] is all sequences (Y_{t-2},...,Y_{t-q})
#   for all q = 2,...,DimMTS (or up to 0, whichever happens first)
# Y_MTS_DENOM_SUM[t,ytm1,y] is the set of sequences Y to sum to get
#   (Y_{t-1}, Y_{t-2},...,Y_{t-q}) = (ytm1, y)
# MTSDenom[t,ytm1,y] is the probability of
#   (Y_{t-1}, Y_{t-2},...,Y_{t-q}) = (ytm1, y)
################################################################################
param Y_MTS_LENMAX {t in 2..T} in 1..(t-1);

set MTSDenomIndex = {t in 2..T, ytm1 in 0..1, q in 1..Y_MTS_LENMAX[t],
    y in YSEQS[q]};
set Y_MTS_DENOM_SUM {(t, ytm1, q, y) in MTSDenomIndex} within YSEQS[T + 1];
param MTSDenom {(t,ytm1,q,y) in MTSDenomIndex} =
    (sum {yy in (Y_MTS_DENOM_SUM[t,ytm1,q,y] inter YHAT)} Q[yy]);
param MTSDenom_Sample {(t,ytm1,q,y) in MTSDenomIndex} =
    (sum {yy in (Y_MTS_DENOM_SUM[t,ytm1,q,y] inter YHAT)} Q_Sample[yy]);

set MTSNumerIndex = {t in 2..T, d in 0..1, ytm1 in 0..1,
    q in 1..Y_MTS_LENMAX[t], y in YSEQS[q]};
set U_MTS_NUMER {(t,d,ytm1,q,y) in MTSNumerIndex} within DomU;
var MTSNumer {(t,d,ytm1,q,y) in MTSNumerIndex} =
    (sum {u in (U_MTS_NUMER[t,d,ytm1,q,y] inter UHAT)} P[u]);
var MTSNumer_H {(t,d,ytm1,q,y) in MTSNumerIndex} =
    (sum {u in (U_MTS_NUMER[t,d,ytm1,q,y] inter UHAT)} H[u]);

set MTSGapIndex = {t in 2..T, d in 0..1, q in 1..Y_MTS_LENMAX[t],
    y in YSEQS[q]};
var MTS_Gap {(t,d,q,y) in MTSGapIndex} =
    (if (
                (MTSDenom[t,0,q,y] > 0)
            and (MTSDenom[t,1,q,y] > 0)
    ) then 1 else 0)
    *(
          MTSNumer[t,d,1,q,y]*MTSDenom[t,0,q,y]
        - MTSNumer[t,d,0,q,y]*MTSDenom[t,1,q,y]
    );
var MTS_Gap_Sample {(t,d,q,y) in MTSGapIndex} =
    (if (
                (MTSDenom_Sample[t,0,q,y] > 0)
            and (MTSDenom_Sample[t,1,q,y] > 0)
    ) then 1 else 0)
    *(
          MTSNumer[t,d,1,q,y]*MTSDenom_Sample[t,0,q,y]
        - MTSNumer[t,d,0,q,y]*MTSDenom_Sample[t,1,q,y]
    );
var MTS_Gap_H {(t,d,q,y) in MTSGapIndex} =
    (if (
                (MTSDenom_Sample[t,0,q,y] > 0)
            and (MTSDenom_Sample[t,1,q,y] > 0)
    ) then 1 else 0)
    *(
          MTSNumer_H[t,d,1,q,y]*MTSDenom_Sample[t,0,q,y]
        - MTSNumer_H[t,d,0,q,y]*MTSDenom_Sample[t,1,q,y]
    );

subject to MTS {(t,d,q,y) in MTSGapIndex}: MTS_Gap[t,d,q,y] >= 0;

param Flag_MTS_In_Objective binary default 0; # Useful for statistical inference

################################################################################
# Convert moment inequality to moment equality by adding a slack variable.
# i.e. MTS_Gap >= 0 if and only if there exists a Slack >= 0 such that
#   MTS_Gap - Slack = 0
# Then treat MTS_Gap - Slack = 0 as a moment equality.
################################################################################
var Slack {(t,d,q,y) in MTSGapIndex} >= 0;
var MTS_Gap_Eq {(t,d,q,y) in MTSGapIndex} = MTS_Gap[t,d,q,y] - Slack[t,d,q,y];

################################################################################
# Reformulate the absolute value of MTS_Gap in the same way as
# what was done for ObsEq_Gap in EqCriterion --- see above
################################################################################
var MTS_Gap_Eq_Abs {(t,d,q,y) in MTSGapIndex} >= 0;
subject to MTS_Gap_Eq_Abs_Pos {(t,d,q,y) in MTSGapIndex}:
    MTS_Gap_Eq_Abs[t,d,q,y] >= MTS_Gap_Eq[t,d,q,y];
subject to MTS_Gap_Eq_Abs_Neg {(t,d,q,y) in MTSGapIndex}:
    MTS_Gap_Eq_Abs[t,d,q,y] >= -1*MTS_Gap_Eq[t,d,q,y];
var IneqCriterion =
    if (Flag_MTS_In_Objective == 1) then
        SQRTN*(sum {(t,d,q,y) in MTSGapIndex} MTS_Gap_Eq_Abs[t,d,q,y])
    else 0;

################################################################################
# Reformulate IneqCriterion_Sample --- which gets used in the constraint
# EstimatedIDSet below --- in the same way as EqCriterion_Sample --- see above
################################################################################
var MTS_Gap_Eq_Sample {(t,d,q,y) in MTSGapIndex}
    = MTS_Gap_Sample[t,d,q,y] - Slack[t,d,q,y];
var MTS_Gap_Eq_Sample_Pos {(t,d,q,y) in MTSGapIndex} >= 0;
var MTS_Gap_Eq_Sample_Neg {(t,d,q,y) in MTSGapIndex} >= 0;
var IneqCriterion_Sample =
    if (Flag_MTS_In_Objective == 1) then SQRTN*(
        sum {(t,d,q,y) in MTSGapIndex} (
              MTS_Gap_Eq_Sample_Pos[t,d,q,y]
            + MTS_Gap_Eq_Sample_Neg[t,d,q,y]
        )
    ) else 0;
subject to MTS_Gap_Eq_Sample_Abs
    {(t,d,q,y) in MTSGapIndex}:
        MTS_Gap_Eq_Sample[t,d,q,y]
        =
        MTS_Gap_Eq_Sample_Pos[t,d,q,y] - MTS_Gap_Eq_Sample_Neg[t,d,q,y];

# For the MTS conditions, Vstar reduces to the following.
var VstarMTS {(t,d,q,y) in MTSGapIndex} =
    (if (
                (MTSDenom[t,0,q,y] > 0)
            and (MTSDenom[t,1,q,y] > 0)
            and (MTSDenom_Sample[t,0,q,y] > 0)
            and (MTSDenom_Sample[t,1,q,y] > 0)
    ) then 1 else 0)
    *(
         MTSNumer[t,d,1,q,y]
        *(MTSDenom[t,0,q,y] - MTSDenom_Sample[t,0,q,y])
        -
         MTSNumer[t,d,0,q,y]
        *(MTSDenom[t,1,q,y] - MTSDenom_Sample[t,1,q,y])
    );
var Slack_H {(t,d,q,y) in MTSGapIndex} default 0;
var MTS_Gap_Eq_CNS {(t,d,q,y) in MTSGapIndex}
    = VstarMTS[t,d,q,y] +  (1/SQRTN)*(MTS_Gap_H[t,d,q,y] - Slack_H[t,d,q,y]);
var MTS_Gap_Eq_CNS_Abs {(t,d,q,y) in MTSGapIndex} >= 0;
subject to MTS_Gap_Eq_CNS_Abs_Pos {(t,d,q,y) in MTSGapIndex}:
    MTS_Gap_Eq_CNS_Abs[t,d,q,y] >= MTS_Gap_Eq_CNS[t,d,q,y];
subject to MTS_Gap_Eq_CNS_Abs_Neg {(t,d,q,y) in MTSGapIndex}:
    MTS_Gap_Eq_CNS_Abs[t,d,q,y] >= -1*MTS_Gap_Eq_CNS[t,d,q,y];

var IneqCriterion_CNS =
    if (Flag_MTS_In_Objective == 1) then
        SQRTN*(sum {(t,d,q,y) in MTSGapIndex} MTS_Gap_Eq_CNS_Abs[t,d,q,y])
    else 0;

param r_MTS_CNS >= 0, default 0;
subject to MTS_CNS {(t,d,q,y) in MTSGapIndex}:
    Slack[t,d,q,y] + Slack_H[t,d,q,y]/SQRTN >= r_MTS_CNS;

################################################################################
# TIV
#
# An implication of the TIV condition in CFHN is that:
#
# P[U_{t} = u, (Y_{0},...,Y_{s-1} = y] = P[U_{s} = u, (Y_{0},...,Y_{s-1} = y]
#
# for all t > s > 1, all u in {0,1}^{2} and all y in {0,1}^{s}
#
# The set U_TIV[t,r,u,y] is then all u such that
#   [U_{t} = u, (Y_{0},Y_{1},...,Y_{r}) = y]
################################################################################
set U_TIV {t in 1..T, r in 0..(t-1), u0 in 0..1, u1 in 0..1, y in YSEQS[r+1]};

subject to TIV {t in 2..T, s in 1..(t-1),
    u0 in 0..1, u1 in 0..1, y in YSEQS[s]}:
        (sum {u in (U_TIV[t,s-1,u0,u1,y] inter UHAT)} P[u])
        =
        (sum {u in (U_TIV[s,s-1,u0,u1,y] inter UHAT)} P[u]);

subject to TIV_CNS {t in 2..T, s in 1..(t-1),
    u0 in 0..1, u1 in 0..1, y in YSEQS[s]}:
        (sum {u in (U_TIV[t,s-1,u0,u1,y] inter UHAT)} H[u])
        =
        (sum {u in (U_TIV[s,s-1,u0,u1,y] inter UHAT)} H[u]);

################################################################################
# Criterion functions
################################################################################
var Criterion = EqCriterion + IneqCriterion;
minimize minCriterion: Criterion;

var Criterion_CNS = EqCriterion_CNS + IneqCriterion_CNS;
minimize minCriterion_CNS: Criterion_CNS;

param CriterionHat default +Infinity;
param Tau default 0;
subject to EstimatedIDSet: # Used when estimating the identified set
    Criterion <= CriterionHat*(1 + Tau);

var Criterion_Sample = EqCriterion_Sample + IneqCriterion_Sample;
subject to EstimatedIDSet_Sample: # Used in CNS to keep sample probabilities
    Criterion_Sample <= CriterionHat*(1 + Tau);

maximize ChangeInPSD {t in 1..T, tt in 1..T}: PSD[t] - PSD[tt];
