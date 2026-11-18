***External Validation - NAMIBIA
use "Namibia.dta", clear
rename *, lower
summ 
drop if fm_d2o>weight
	gen ffm_d2o_newhy = tbw / 0.732
	gen fm_d2o_newhy = weight - ffm_d2o_newhy

gen htm =  height/100
rename age Age
gen male = 0
gen BA = 1
gen SA = 0
gen AO = 0
gen Other=0

gen ageg = 1 if Age>=16 & Age <25
replace ageg = 2 if Age >=25 & Age < 40
replace ageg = 3 if Age >=40 & Age < 65
replace ageg = 4 if Age >=65 & Age < 80
replace ageg = 5 if Age >=80 & !missing(Age)


*Summary Stats of Age, ht, wt, fm, ffm, male, Ethnicity overall
foreach x of varlist Age htm weight ffm_d2o fm_d2o {
qui summ `x' ,d
mat `x' = r(p50), r(p25), r(p75)
}

count if male==1
scalar male_n = r(N)
count
scalar n = r(N)
mat male = male_n, male_n/n *100,.

gen eth_group = 2

foreach z of numlist 1/5 {
	count if eth_group==`z'
	scalar eth = r(N)
	count
	scalar n = r(N)
	mat eth`z' = eth, eth/n *100,.
}
mat ethnic = eth1\eth2\eth3\eth4\eth5

mat Tab2b = Age \ htm \ weight \ ffm_d2o \ fm_d2o \ male \ ethnic


gen pred_lnffm = 0.4621*htm - 322*(weight^-2) - 3.164*(weight^-0.5) + 0.004578*weight  - 0.00000979*((weight)^2) - 0.0011*Age + 0.000006*(weight*Age) + 0.01*BA + 0.0113*SA - 0.0142*AO +0.0184*Other + 361*((htm*Age)^-2) + 0.0262*((ht*Age/100)^-2)*ln(ht*Age/100) + 3.3122 if male==1

replace pred_lnffm = 0.5392*htm - 322*(weight^-2) + 0.004578*weight - 0.0011*Age + 0.000006*(weight*Age) + 0.01*BA + 0.0113*SA - 0.0142*AO +0.0184*Other + 361*(htm*Age)^-2 + 0.0262*((ht*Age/100)^-2)*ln(ht*Age/100) + 2.6045 if male==0


gen pred_ffm = exp(pred_lnffm)
gen pred_fm = weight - pred_ffm

keep if !missing(pred_fm) 


*overall 
regress fm_d2o pred_fm
matrix temp = r(table)'
matrix slope_fm = temp[1,1], temp[1,5], temp[1,6]

constraint 1 pred_fm = 1
cnsreg fm_d2o pred_fm, constraints(1)
matrix temp = r(table)'
matrix cil_fm = temp[2,1], temp[2,5], temp[2,6]

pmcalplot pred_fm fm_d2o, cont nohist legend(off) xlab(0(20)100) ylab(0(20)100) xtitle(Predicted FM (kg)) ytitle(DD Assessed FM (kg)) title(Namibia) range(1 100) nostat

*R2
corr fm_d2o pred_fm 
scalar r2_fm = r(rho)^2 
scalar r2_fm_var = (4*r2_fm*((1-r2_fm)^2))/r(N)
matrix r2_fm = r2_fm*100 , (r2_fm  - 1.96*(sqrt(r2_fm_var)))*100 , (r2_fm  + 1.96*(sqrt(r2_fm_var)))*100 

*RMSE
capture drop diff
gen diff = (fm_d2o -  pred_fm)^2
summ diff
di sqrt((r(sum) / r(N)))
matrix rmse_fm =  sqrt((r(sum) / r(N)))

mat stats_fm = r2_fm , slope_fm , cil_fm , rmse_fm
mat colnames stats_fm = r2_fm  r2_fm_low r2_fm_up slope_fm  slope_fm_low slope_fm_up cil_fm  cil_fm_low cil_fm_up rmse_fm


**RECALIB OF INTERCEPT
gen lnffm_d2o = ln(ffm_d2o)

gen pred_lnffm_recal =.
constraint 1 pred_lnffm = 1 if male==0
cnsreg lnffm_d2o pred_lnffm if male==0, constraints(1)
di _b[_cons]

replace pred_lnffm_recal = pred_lnffm - 0.08942182

gen pred_fm_recal = weight - exp(pred_lnffm_recal)
pmcalplot pred_fm_recal fm_d2o, cont nohist legend(off) xlab(0(20)100) ylab(0(20)100) xtitle(Predicted FM (kg) - AFTER RECALIBRATION) ytitle(DD Assessed FM (kg)) title(Namibia - after recal)

corr fm_d2o pred_fm_recal 
scalar r2_fm = r(rho)^2 
scalar r2_fm_var = (4*r2_fm*((1-r2_fm)^2))/r(N)
matrix r2_fm = r2_fm*100 , (r2_fm  - 1.96*(sqrt(r2_fm_var)))*100 , (r2_fm  + 1.96*(sqrt(r2_fm_var)))*100 

regress fm_d2o pred_fm_recal
constraint 1 pred_fm_recal = 1
cnsreg fm_d2o pred_fm_recal, constraints(1)

capture drop diff
gen diff = (fm_d2o -  pred_fm_recal)^2
summ diff
di sqrt((r(sum) / r(N)))
matrix rmse_fm =  sqrt((r(sum) / r(N)))
