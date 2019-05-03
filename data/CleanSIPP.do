********************************************************************************
* CleanSIPP.do
********************************************************************************
clear all
set more off
cap log close
log using CleanSIPP.log, append
********************************************************************************
global initial_wave = 8
global final_wave = 14
********************************************************************************
* For each wave:
*   Keep variables of potential use
*   Create unique unit numeric id
*   Keep only final reference month
*   Remove type z observations
********************************************************************************
forvalues j = $initial_wave/$final_wave {
    display "Loading in raw data set for wave `j'"

    local filename_in = "sippl08puw`j'.dta"
    use `filename_in', clear

    keep ssuid epppnum swave srefmon rhcalmn eppintvw ///
        tage esex erace eorigin edisprev eeducate ///
        renroll eafnow eptwrk rmesr tfipsst tpearn ///
        tpyrate* apyrate* ejobcntr ///
        epayhr1 apayhr1 ///
        eptwrk aptwrk ///
        ehrsall ahrsall ///

    * Create a unique sipp id
    egen sippid = concat(ssuid epppnum)
    drop ssuid epppnum
    order sippid swave

    * Only keep the final reference month to mitigate seam bias
    keep if srefmon == 4
    drop srefmon

    * Get rid of type z observations (those that are highly imputed)
    * This also appears to get rid of children but I will drop them explicitly
    * in the next file anyway
    keep if eppintvw <= 2
    drop eppintvw

    local filename_out = "sipp08_w`j'.dta"
    save `filename_out', replace
}

********************************************************************************
* Create a merged dataset with all waves
* Select by age
* Create employment variables
********************************************************************************
local filename_in = "sipp08_w" + "$initial_wave" + ".dta"
use `filename_in', clear

local i = $initial_wave + 1
while `i' <= $final_wave {
    display "Loading in wave `i'"
    local filename_in = "sipp08_w`i'.dta"
    qui append using `filename_in'
    local i = `i' + 1
}

sort sippid swave
* Convenient for the unit id to be numeric and small
egen id = group(sippid)
drop sippid
order id swave
sort id swave

* Keep only those between 18 and 55 during the initial wave
by id: egen initage = min(tage)
tab initage
keep if (initage >= 18)*(initage <= 55)

* Labor force status variable based on
* rmesr (Employment status recode)
********************************************************************************
/*1. With a job entire month, worked all weeks*/
/*2. With a job entire month, absent from work without pay 1+ weeks,
absence not due to layoff*/
/*3. With a job entire month, absent from work without pay 1+ weeks, absence*/
/*due to layoff*/
/*4. With a job at least 1 but not all weeks, no time on layoff and no time*/
/*looking for work*/
/*5. With a job at least 1 but not all weeks, remaining weeks on layoff or*/
/*looking for work*/
/*6. No job all month, on layoff or looking for work all weeks*/
/*7. No job all month, at least one but not all weeks on layoff or looking for*/
/*work*/
/*8 .No job all month, no time on layoff and no time looking for work. **/
********************************************************************************
* 10/19/18: Chetty (2008) appears to classify employed as rmesr = 1 or 2.
* I am going to follow him on that.
gen emp = 1*(rmesr <= 2) + 2*(rmesr == 4 | rmesr == 8)
label define emplabel ///
    2 "Not participating" ///
    0 "Unemployed" ///
    1 "Employed"
label values emp emplabel

local filename_out = "sipp08_merged.dta"
save `filename_out', replace

*###############################################################################
* Create a balanced panel
* Do some additional cleaning based on time-varying criteria:
*   Remove those who were disabled, in school or in the military at any
*   time
*   Remove those whose education level changed over the time
*   (which indicates schooling or measurement error)
*   Select on education
*   Create some variables about employment transitions
*###############################################################################
* Balance the panel
sort id swave
by id: gen nwaves = _N
by id: gen nvals = (_n == 1)
display "Cross-sectional observations before balancing:"
count if nvals
keep if nwaves == $final_wave - $initial_wave + 1
display "Cross-sectional observations after balancing:"
count if nvals

* Remove those whose sex variable changed
by id (esex), sort: gen sexchg = (esex[1] != esex[_N])
tab sexchg if swave==$initial_wave
drop if sexchg

* Keep only men
tab esex if swave==$initial_wave
keep if (esex == 1)

* Keep only men who were in the labor force over all periods
sort id swave
by id: egen periodsnotpartic = sum((emp==2))
tab periodsnotpartic if swave==$initial_wave
drop if (periodsnotpartic > 0)

* Remove men who ever reported a work-preventing disability
by id: egen ndisabled = total((edisprev == 1))
tab ndisabled if swave==$initial_wave
drop if ndisabled > 0
drop ndisabled edisprev

* Remove men who were ever enrolled in school
* (should have been taken care of by labor force nonparticipation)
by id: egen nschool = total((renroll < 3))
tab nschool if swave==$initial_wave
drop if nschool > 0
drop nschool

* Remove those who were ever in the military
* (might have been taken care of by labor force nonparticipation?)
by id: egen naf = total((eafnow == 1))
tab naf if swave==$initial_wave
drop if naf > 0
drop eafnow naf

* Keep men with a high school education
* Also remove people whose education level changed through the sample
by id: egen maxeducate = max(eeducate)
by id: egen mineducate = min(eeducate)
drop if (maxeducate != mineducate)
tab maxeducate if swave==$initial_wave
/*38          "12th grade, no diploma"        */
/*39          "High School Graduate - (diploma"*/
/*40          "Some college, but no degree"   */
/*41          "Diploma or certificate from a" */
/*43          "Associate (2-yr) college degree"*/
/*44          "Bachelor's degree (for example:"*/
/*45          "Master's degree (For example: MA,"*/
drop if (maxeducate < 39) | (maxeducate >= 43)
drop maxeducate mineducate

* Count number of cross-sectional observations remaining
display "Cross-sectional observations remaining:"
count if nvals

* Create some variables about employment transitions
xtset id swave
by id: gen demp = D.emp
* Did your employment status change?
gen cemp = abs(demp)
order id swave emp demp cemp
* Total number of transitions, i.e. emp != L.emp
by id: gen ntrans = sum(cemp)

local filename_out = "sipp08_clean.dta"
save `filename_out', replace
********************************************************************************
********************************************************************************
********************************************************************************
* Create tab-delimited datasets for reading into Matlab
********************************************************************************
********************************************************************************
********************************************************************************
keep emp id swave initage tage
replace swave = swave - $initial_wave // Start time at 0
label drop _all // Make all variables numeric

preserve
by id: gen nvals = (_n == 1)
display "Cross-sectional observations in full sample:"
count if nvals
display "Cross-sectional observations in young sample:"
count if nvals & initage <= 40

outsheet using "sipp08.tsv", replace
outsheet using "sipp08-young.tsv" if initage <= 40, replace
restore

drop tage initage
reshape wide emp, i(id) j(swave) // wide form
label drop _all
outsheet emp* using "sipp08-wide.tsv", nonames replace

log close
exit
