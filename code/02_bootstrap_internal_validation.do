*** Bootstrap Internal Validation***
timer on 1
foreach a of numlist 1/1000 {
di "SAMPLE" `a'
qui{
use temp_analysis_file.dta, clear

bsample, strata(WHORegion male ageg eth_group)
append using temp_analysis_file.dta, gen(orig)

xi: mfp, alpha(0.05) select(0.157) center(no) all: regress lnDD_ffm_IAEACons htm weight Age BA SA AO Other i.male (WHO_Afr WHO_EMed WHO_Eur WHO_SEAsia WHO_WPac) ht_Age	ht_male		wt_Age	wt_male if orig==0
local rsq_apparent=e(r2)
local fp = e(fp_fvl)
predict pred
qui regress lnDD_ffm_IAEACons pred if orig==0
matrix A = r(table)
local calib = A[1,1]
local cil = A[1,2]

regress lnDD_ffm_IAEACons `fp' if orig == 1
local rsq_test=e(r2)
predict pred2 
qui regress lnDD_ffm_IAEACons pred2 if orig==1
matrix A = r(table)
local calib_ = A[1,1]
local cil_ = A[1,2]

local calib_opt = `calib' - `calib_'
local cil_opt = `cil' - `cil_'
local rsq_diff = `rsq_test' - `rsq_apparent'

di `calib'
di `rsq_apparent'
di `rsq_test'
di `rsq_diff'
di `calib_opt'
di `cil_opt'


preserve
clear
use "opt_r2_final model", clear

capture gen opt_r2 = .
capture gen calib = .
capture gen calib_ = .
capture gen opt_calib = .

capture gen cil = .
capture gen cil_ = .
capture gen opt_cil = .


replace opt_r2 = `rsq_diff' in `a'
replace calib = `calib' in `a'
replace calib_ = `calib_' in `a'
replace opt_calib = `calib_opt' in `a'

replace cil = `cil' in `a'
replace cil_ = `cil_' in `a'
replace opt_cil = `cil_opt' in `a'

save "opt_r2_final model",replace emptyok
restore
}
}


*Optimism in apparent performance stats
preserve
use "opt_r2_final model", clear
summ opt_r2 
matrix opt_r2 = r(mean) * 100
summ opt_calib
matrix opt_calib = r(mean)
summ calib
matrix calib = r(mean)
summ opt_cil
matrix opt_cil = r(mean)

matrix bootstrap = opt_r2 \ opt_calib \ opt_cil
mat rownames bootstrap = opt_r2 opt_calib opt_cil
restore

