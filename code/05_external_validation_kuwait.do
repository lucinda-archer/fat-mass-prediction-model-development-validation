***External Validation - Kuwait Wellness Study DXA
use "Kuwait.dta", clear

gen htm = height/100
gen fm_dxa =  Total_fat_mass_g/1000
gen fmp_dxa = fm_dxa/weight *100
gen male = 1 if Sex=="M"
	replace male = 0 if Sex=="F"
gen Other = 1 
gen BA= 0 
gen SA = 0 
gen AO = 0	
gen ageg = 1 if Age>=16 & Age <25
replace ageg = 2 if Age >=25 & Age < 40
replace ageg = 3 if Age >=40 & Age < 65
replace ageg = 4 if Age >=65 & Age < 80
replace ageg = 5 if Age >=80 & !missing(Age)

*Exclusions
drop if missing(fm_dxa)
drop if weight < 10
drop if height<110
drop if fm_dxa>weight
drop if fmp_dxa>75

gen ffm_dxa = weight - fm_dxa

*Summary Stats of Age, ht, wt, fm, ffm, male, Ethnicity overall
foreach x of varlist Age htm weight ffm_dxa fm_dxa {
qui summ `x' ,d
mat `x' = r(p50), r(p25), r(p75)
}

count if male==1
scalar male_n = r(N)
count
scalar n = r(N)
mat male = male_n, male_n/n *100,.

gen eth_group=5
foreach z of numlist 1/5 {
	count if eth_group==`z'
	scalar eth = r(N)
	count
	scalar n = r(N)
	mat eth`z' = eth, eth/n *100,.
}
mat ethnic = eth1\eth2\eth3\eth4\eth5

mat Tab2c = Age \ htm \ weight_kg \ ffm_dxa \ fm_dxa \ male \ ethnic

gen pred_lnffm = 0.4621*htm - 322*(weight^-2) - 3.164*(weight^-0.5) + 0.004578*weight  - 0.00000979*((weight)^2) - 0.0011*Age + 0.000006*(weight*Age) + 0.01*BA + 0.0113*SA - 0.0142*AO +0.0184*Other + 361*((htm*Age)^-2) + 0.0262*((ht*Age/100)^-2)*ln(ht*Age/100) + 3.2535 if male==1

replace pred_lnffm = 0.5392*htm - 322*(weight^-2) + 0.004578*weight - 0.0011*Age + 0.000006*(weight*Age) + 0.01*BA + 0.0113*SA - 0.0142*AO +0.0184*Other + 361*(htm*Age)^-2 + 0.0262*((ht*Age/100)^-2)*ln(ht*Age/100) + 2.5458 if male==0


gen pred_ffm = exp(pred_lnffm)
gen pred_fm = weight - pred_ffm


*overall 
regress fm_dxa pred_fm
matrix temp = r(table)'
matrix slope_fm = temp[1,1], temp[1,5], temp[1,6]

constraint 1 pred_fm = 1
cnsreg fm_dxa pred_fm, constraints(1)
matrix temp = r(table)'
matrix cil_fm = temp[2,1], temp[2,5], temp[2,6]

pmcalplot pred_fm fm_dxa, cont nohist legend(off) xlab(0(20)100) ylab(0(20)100) xtitle(Predicted FM (kg)) ytitle(DXA Assessed FM (kg)) title(Kuwait) nostat range(0 100)

*R2
corr fm_dxa pred_fm 
scalar r2_fm = r(rho)^2 
scalar r2_fm_var = (4*r2_fm*((1-r2_fm)^2))/r(N)
matrix r2_fm = r2_fm*100 , (r2_fm  - 1.96*(sqrt(r2_fm_var)))*100 , (r2_fm  + 1.96*(sqrt(r2_fm_var)))*100 

*RMSE
capture drop diff
gen diff = (fm_dxa -  pred_fm)^2
summ diff
di sqrt((r(sum) / r(N)))
matrix rmse_fm =  sqrt((r(sum) / r(N)))

mat stats_fm = r2_fm , slope_fm , cil_fm , rmse_fm
mat colnames stats_fm = r2_fm  r2_fm_low r2_fm_up slope_fm  slope_fm_low slope_fm_up cil_fm  cil_fm_low cil_fm_up rmse_fm




**By Sex
*males 
preserve
keep if male==1
regress fm_dxa pred_fm
matrix temp = r(table)'
matrix slope_fm = temp[1,1], temp[1,5], temp[1,6]

constraint 1 pred_fm = 1
cnsreg fm_dxa pred_fm, constraints(1)
matrix temp = r(table)'
matrix cil_fm = temp[2,1], temp[2,5], temp[2,6]

pmcalplot pred_fm fm_dxa, cont nohist legend(off) xtitle(Predicted FM (kg)) ytitle(DXA Assessed FM (kg)) title(Males) nostat

*R2
corr fm_dxa pred_fm 
scalar r2_fm = r(rho)^2 
scalar r2_fm_var = (4*r2_fm*((1-r2_fm)^2))/r(N)
matrix r2_fm = r2_fm*100 , (r2_fm  - 1.96*(sqrt(r2_fm_var)))*100 , (r2_fm  + 1.96*(sqrt(r2_fm_var)))*100 

*RMSE
capture drop diff
gen diff = (fm_dxa -  pred_fm)^2
summ diff
di sqrt((r(sum) / r(N)))
matrix rmse_fm =  sqrt((r(sum) / r(N)))

mat stats_fm_m = r2_fm , slope_fm , cil_fm , rmse_fm
mat colnames stats_fm_m = r2_fm  r2_fm_low r2_fm_up slope_fm  slope_fm_low slope_fm_up cil_fm  cil_fm_low cil_fm_up rmse_fm
restore

*females 
preserve
keep if male==0
regress fm_dxa pred_fm
matrix temp = r(table)'
matrix slope_fm = temp[1,1], temp[1,5], temp[1,6]

constraint 1 pred_fm = 1
cnsreg fm_dxa pred_fm, constraints(1)
matrix temp = r(table)'
matrix cil_fm = temp[2,1], temp[2,5], temp[2,6]

pmcalplot pred_fm fm_dxa, cont nohist legend(off) xtitle(Predicted FM (kg)) ytitle(DXA Assessed FM (kg)) title(Females) nostat

*R2
corr fm_dxa pred_fm 
scalar r2_fm = r(rho)^2 
scalar r2_fm_var = (4*r2_fm*((1-r2_fm)^2))/r(N)
matrix r2_fm = r2_fm*100 , (r2_fm  - 1.96*(sqrt(r2_fm_var)))*100 , (r2_fm  + 1.96*(sqrt(r2_fm_var)))*100 

*RMSE
capture drop diff
gen diff = (fm_dxa -  pred_fm)^2
summ diff
di sqrt((r(sum) / r(N)))
matrix rmse_fm =  sqrt((r(sum) / r(N)))

mat stats_fm_f = r2_fm , slope_fm , cil_fm , rmse_fm
mat colnames stats_fm_f = r2_fm  r2_fm_low r2_fm_up slope_fm  slope_fm_low slope_fm_up cil_fm  cil_fm_low cil_fm_up rmse_fm
restore

gen sex = "Female" if male==0
replace sex = "Male" if male==1
pmcalplot pred_fm fm_dxa, by(sex) cont nohist legend(off) xtitle(Predicted FM (kg)) ytitle(DXA Assessed FM (kg)) nostat range(0 100)


**By Ageg
foreach x of numlist 1/4 {
pmcalplot pred_fm fm_dxa if ageg==`x', cont nohist legend(off) xtitle(Predicted FM (kg)) ytitle(DXA Assessed FM (kg)) title(Age Group `x') nostat
}

gen ageg2 = "16-24" if ageg==1
replace ageg2 = "25-39" if ageg==2
replace ageg2 = "40-64" if ageg==3

pmcalplot pred_fm fm_dxa , by(ageg2, col(3)) cont nohist legend(off) xtitle(Predicted FM (kg)) ytitle(DXA Assessed FM (kg)) nostat range(0 100)

*RMSE
foreach x of numlist 1/3 {
preserve
keep if ageg==`x'
capture drop diff
gen diff = (fm_dxa -  pred_fm)^2
summ diff
di sqrt((r(sum) / r(N)))
matrix rmse_fm_a`x' =  sqrt((r(sum) / r(N)))
restore
}
mat rmse_fm_allageg = rmse_fm_a1 \ rmse_fm_a2 \ rmse_fm_a3




**by age
foreach x of numlist 1/4 {
preserve
keep if ageg==`x'
regress fm_dxa pred_fm
matrix temp = r(table)'
matrix slope_fm = temp[1,1], temp[1,5], temp[1,6]

constraint 1 pred_fm = 1
cnsreg fm_dxa pred_fm, constraints(1)
matrix temp = r(table)'
matrix cil_fm = temp[2,1], temp[2,5], temp[2,6]

*R2
corr fm_dxa pred_fm 
scalar r2_fm = r(rho)^2 
scalar r2_fm_var = (4*r2_fm*((1-r2_fm)^2))/r(N)
matrix r2_fm = r2_fm*100 , (r2_fm  - 1.96*(sqrt(r2_fm_var)))*100 , (r2_fm  + 1.96*(sqrt(r2_fm_var)))*100 

*RMSE
capture drop diff
gen diff = (fm_dxa -  pred_fm)^2
summ diff
di sqrt((r(sum) / r(N)))
matrix rmse_fm =  sqrt((r(sum) / r(N)))

mat stats_fm_a`x' = r2_fm , slope_fm , cil_fm , rmse_fm
mat colnames stats_fm_a`x' = r2_fm  r2_fm_low r2_fm_up slope_fm  slope_fm_low slope_fm_up cil_fm  cil_fm_low cil_fm_up rmse_fm
restore
}
