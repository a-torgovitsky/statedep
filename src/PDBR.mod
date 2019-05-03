################################################################################
# PDBR.mod
#
# Dynamic random effects probit
################################################################################
load amplgsl_64.dll;
function gsl_cdf_ugaussian_P;
function gsl_ran_ugaussian_pdf;
function gsl_cdf_ugaussian_Pinv;
param pi = 4 *atan(1); # AMPL does not have pi, but does have arctan...ok

param N integer > 0; # Cross-section
param T integer >= 1; # Time
param J integer >= 2; # Number of nodes for integration
param K integer >= 0; # Covariates

param x {i in 1..N, t in 1..T, k in 1..K}; # Covariates
param y {i in 1..N, t in 0..T} binary; # Binary outcomes
param ghnodes {j in 1..J}; # Gauss-Hermite nodes
param ghweights {j in 1..J} >= 0; # Gauss-Hermite weights

var gamma, default 0;
var beta {k in 1..K}, default 0;
var lambda, default 0;
var sigma >= 0, default 1;

var Idx {d in 0..1, i in 1..N, j in 1..J, t in 1..T} =
      gamma*d
    + (sum {k in 1..K} beta[k]*x[i,t,k])
    + lambda*y[i,0]
    + sqrt(2)*sigma*ghnodes[j];

var LHij {i in 1..N, j in 1..J} =
    prod {t in 1..T} (
        y[i,t]*gsl_cdf_ugaussian_P(Idx[y[i,t-1],i,j,t])
        +
        (1-y[i,t])*(1 - gsl_cdf_ugaussian_P(Idx[y[i,t-1],i,j,t]))
    );

var LHi {i in 1..N} =
    (1/sqrt(pi))
    *
    sum {j in 1..J} (
        ghweights[j]
        *LHij[i,j]
    );

maximize LLH:
    sum {i in 1..N}(
        log(LHi[i])
    );

################################################################################
# TARGET PARAMETERS
################################################################################
# Average over draws of the random effect for each (i,t)
var ASFit {d in 0..1, i in 1..N, t in 1..T} =
    (1/sqrt(pi))*
    (sum {j in 1..J}
        ghweights[j]*gsl_cdf_ugaussian_P(Idx[d,i,j,t])
    );

var ATE {t in 1..T+1} =
    if (t <= T) then (
        (1/N)*sum{i in 1..N} (
            ASFit[1,i,t] - ASFit[0,i,t]
        )
    ) else (
        (1/T)*(sum {tt in 1..T} ATE[tt])
    );

var PSD {t in 1..T+1} =
    if (t <= T) then (
        if (gamma <= 0) then 0 else ATE[t]
    ) else (
        (1/T)*(sum {tt in 1..T} PSD[tt])
    );

var NSD {t in 1..T+1} =
    if (t <= T) then (
        if (gamma >= 0) then 0 else ATE[t]
    ) else (
        (1/T)*(sum {tt in 1..T} NSD[tt])
    );

var TSD {t in 1..T+1} =
    if (t <= T) then (
        PSD[t] + NSD[t]
    ) else (
        (1/T)*(sum {tt in 1..T} TSD[tt])
    );

# Conditional probabilities of Y used below
param PrYt {t in 0..T} = (1/N)*(sum {i in 1..N} y[i,t]);
param PrYt_given_Ytm1 {t in 1..T, d in 0..1} =
    (1/N)*(sum {i in 1..N} (
          y[i,t]
        * (if (y[i,t-1] == d) then 1 else 0)
    ))
    /
    (if (d == 1) then PrYt[t-1] else (1 - PrYt[t-1]));

# Parameters conditional on Y_{t} = 0
var PrUd_and_Yt0 {d in 0..1, t in 1..T} =
    (1/N)*(sum {i in 1..N} (
        (if (gamma*d > gamma*y[i,t-1]) then 1 else 0)
        *
        (ASFit[d,i,t] - ASFit[y[i,t-1],i,t])
    ));

var PrUd_given_Yt0 {d in 0..1, t in 1..T} =
    PrUd_and_Yt0[d,t] / (1 - PrYt[t]);

var ATE_G0 {t in 1..T+1} =
    if (t <= T) then (
        PrUd_given_Yt0[1,t] - PrUd_given_Yt0[0,t]
    ) else (
        (1/T)*sum{tt in 1..T} ATE_G0[tt]
    );

var PSD_G0 {t in 1..T+1} =
    if (t <= T) then (
        if (gamma <= 0) then 0 else ATE_G0[t]
    ) else (
        (1/T)*sum{tt in 1..T} PSD_G0[tt]
    );

var NSD_G0 {t in 1..T+1} =
    if (t <= T) then (
        if (gamma >= 0) then 0 else ATE_G0[t]
    ) else (
        (1/T)*sum{tt in 1..T} NSD_G0[tt]
    );

# Parameters conditional on Y_{t} = 0, Y_{t-1} = 0
var PrUd_and_Yt0_and_Ytm10 {d in 0..1, t in 1..T} =
    (1/N)*(sum {i in 1..N} (
          (if (gamma*d > gamma*y[i,t-1]) then 1 else 0)
        * (if (y[i,t-1] == 0) then 1 else 0)
        * (ASFit[d,i,t] - ASFit[y[i,t-1],i,t])
    ));

var PrUd_given_Yt0_and_Ytm10 {d in 0..1, t in 1..T} =
    PrUd_and_Yt0_and_Ytm10[d,t]
    /
    (
          (1 - PrYt_given_Ytm1[t,0])
        * (1 - PrYt[t-1])
    );

var ATE_G00 {t in 1..T+1} =
    if (t <= T) then (
        PrUd_given_Yt0_and_Ytm10[1,t] - PrUd_given_Yt0_and_Ytm10[0,t]
    ) else (
        (1/T)*sum{tt in 1..T} ATE_G00[tt]
    );

var PSD_G00 {t in 1..T+1} =
    if (t <= T) then (
        if (gamma <= 0) then 0 else ATE_G00[t]
    ) else (
        (1/T)*sum{tt in 1..T} PSD_G00[tt]
    );

var NSD_G00 {t in 1..T+1} =
    if (t <= T) then (
        if (gamma >= 0) then 0 else ATE_G00[t]
    ) else (
        (1/T)*sum{tt in 1..T} NSD_G00[tt]
    );

# Parameters conditional on Y_{t} = 1
var PrUd_and_Yt1 {d in 0..1, t in 1..T} =
    (1/N)*(sum {i in 1..N} (
        if (gamma >= 0)
        then ASFit[min(d, y[i,t-1]),i,t]
        else ASFit[max(d, y[i,t-1]),i,t]
    ));

var PrUd_given_Yt1 {d in 0..1, t in 1..T} =
    PrUd_and_Yt1[d,t]/PrYt[t];

var ATE_G1 {t in 1..T+1} =
    if (t <= T) then (
        PrUd_given_Yt1[1,t] - PrUd_given_Yt1[0,t]
    ) else (
        (1/T)*sum{tt in 1..T} ATE_G1[tt]
    );

var PSD_G1 {t in 1..T+1} =
    if (t <= T) then (
        if (gamma <= 0) then 0 else ATE_G1[t]
    ) else (
        (1/T)*sum{tt in 1..T} PSD_G1[tt]
    );

var NSD_G1 {t in 1..T+1} =
    if (t <= T) then (
        if (gamma >= 0) then 0 else ATE_G1[t]
    ) else (
        (1/T)*sum{tt in 1..T} NSD_G1[tt]
    );

# Parameters conditional on Y_{t} = 1, Y_{t-1} = 1
var PrUd_and_Yt1_and_Ytm11 {d in 0..1, t in 1..T} =
    (1/N)*(sum {i in 1..N} (
        (if (gamma >= 0)
         then ASFit[min(d, y[i,t-1]),i,t]
         else ASFit[max(d, y[i,t-1]),i,t])
        * (if (y[i,t-1] == 1) then 1 else 0)
    ));

var PrUd_given_Yt1_and_Ytm11 {d in 0..1, t in 1..T} =
    PrUd_and_Yt1_and_Ytm11[d,t]
    /
    (PrYt_given_Ytm1[t,1]*PrYt[t-1]);

var ATE_G11 {t in 1..T+1} =
    if (t <= T) then (
        PrUd_given_Yt1_and_Ytm11[1,t] - PrUd_given_Yt1_and_Ytm11[0,t]
    ) else (
        (1/T)*sum {tt in 1..T} ATE_G11[tt]
    );

var PSD_G11 {t in 1..T+1} =
    if (t <= T) then (
        if (gamma <= 0) then 0 else ATE_G11[t]
    ) else (
        (1/T) * sum {tt in 1..T} PSD_G11[tt]
    );

var NSD_G11 {t in 1..T+1} =
    if (t <= T) then (
        if (gamma >= 0) then 0 else ATE_G11[t]
    ) else (
        (1/T) * sum {tt in 1..T} NSD_G11[tt]
    );
