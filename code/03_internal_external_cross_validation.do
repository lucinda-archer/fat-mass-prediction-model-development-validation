***INTERNAL-EXTERNAL VALIDATION***
*By REGION
xi: mfp, alpha(0.05) select(0.157) center(no) all: regress lnDD_ffm_IAEACons htm weight Age BA SA AO Other i.male (WHO_Afr WHO_EMed WHO_Eur WHO_SEAsia WHO_WPac) ht_Age	ht_male		wt_Age	wt_male
	
capture drop pred_*
foreach R of numlist 1/6 {
di `R'

regress lnDD_ffm_IAEACons htm  Iweig__1 Iweig__2 Age BA SA AO Other male WHO_Afr WHO_EMed WHO_Eur WHO_SEAsia WHO_WPac 
			Iht_A__1 Iht_A__2	ht_male wt_Age  Iwt_m__1 Iwt_m__2	if WHOR_c !=`R'
				 
matrix beta`R' = e(b)
scalar shrink`R' = (-2*(e(ll_0)-e(ll)) - 10) / (-2*(e(ll_0)-e(ll)))
di shrink`R'   
matrix beta2`R' = shrink`R'*beta`R'

gen pred_lnffm_R`R' = beta2`R'[1,1]*htm + beta2`R'[1,2]*Iweig__1 + beta2`R'[1,3]*Iweig__2 + beta2`R'[1,4]*Age + beta2`R'[1,5]*BA + beta2`R'[1,6]*SA + beta2`R'[1,7]*AO + 		beta2`R'[1,8]*Other + beta2`R'[1,9]*male + beta2`R'[1,10]*WHO_Afr + beta2`R'[1,11]*WHO_EMed + beta2`R'[1,12]*WHO_Eur + beta2`R'[1,13]*WHO_SEAsia 	+ 				beta2`R'[1,14]*WHO_WPac	+ beta2`R'[1,15]*Iht_A__1 + beta2`R'[1,16]*Iht_A__2 + beta2`R'[1,17]*ht_male + beta2`R'[1,18]*wt_Age + beta2`R'[1,19]*Iwt_m__1 + 				beta2`R'[1,20]*Iwt_m__2 + beta2`R'[1,21]
				
				
gen pred_ffm_R`R' = exp(pred_lnffm_R`R') if WHOR_c==`R'
gen pred_fm_R`R' = weight - pred_ffm_R`R' if WHOR_c==`R'


*OVERALL PERFORMANCE
regress DD_fm_IAEACons pred_fm_R`R'  if WHOR_c==`R'
matrix temp = r(table)'
matrix calib_fm`R' = temp[1,1], temp[1,5], temp[1,6]

constraint 1 pred_fm_R`R' = 1
cnsreg DD_fm_IAEACons pred_fm_R`R' if WHOR_c==`R', constraints(1)
matrix cil_fm`R' = temp[2,1], temp[2,5], temp[2,6]

qui corr pred_fm_R`R' DD_fm_IAEACons if WHOR_c==`R'
scalar r2 = r(rho)^2
scalar r2_var = (4*r2*((1-r2)^2))/r(N)
scalar r = r(rho)
scalar n = r(N)
matrix r2_fm`R' = r2, r2_var
matrix z_fm`R' = 0.5 * (ln((1 + r) / (1 - r))) , 1 / sqrt(n - 3)


capture drop diff
gen diff = (DD_fm_IAEACons -  pred_fm_R`R')^2  if WHOR_c==`R'
summ diff  if WHOR_c==`R'
matrix rmse_fm`R' =  sqrt((r(sum) / r(N)))  , . , . , .



*BY SEX
*MALES
regress DD_fm_IAEACons pred_fm_R`R'  if WHOR_c==`R' & male==1
matrix temp = r(table)'
matrix calib_fm`R'_male = temp[1,1], temp[1,5], temp[1,6]

constraint 1 pred_fm_R`R' = 1
cnsreg DD_fm_IAEACons pred_fm_R`R' if WHOR_c==`R' & male==1, constraints(1)
matrix cil_fm`R'_male = temp[2,1], temp[2,5], temp[2,6] 

qui corr pred_fm_R`R' DD_fm_IAEACons if WHOR_c==`R' & male==1
scalar r2 = r(rho)^2
scalar r2_var = (4*r2*((1-r2)^2))/r(N)
matrix r2_fm`R'_male = r2, r2_var
scalar r = r(rho)
scalar n = r(N)
matrix z_fm`R'_male = 0.5 * (ln((1 + r) / (1 - r))) , 1 / sqrt(n - 3)

capture drop diff
gen diff = (DD_fm_IAEACons -  pred_fm_R`R')^2  if WHOR_c==`R' & male==1
summ diff  if WHOR_c==`R' & male==1
matrix rmse_fm`R'_male =  sqrt((r(sum) / r(N)))  , . , . , .

*FEMALES
regress DD_fm_IAEACons pred_fm_R`R'  if WHOR_c==`R' & male==0
matrix temp = r(table)'
matrix calib_fm`R'_female = temp[1,1], temp[1,5], temp[1,6]

constraint 1 pred_fm_R`R' = 1
cnsreg DD_fm_IAEACons pred_fm_R`R' if WHOR_c==`R' & male==0, constraints(1)
matrix cil_fm`R'_female = temp[2,1], temp[2,5], temp[2,6] 

qui corr pred_fm_R`R' DD_fm_IAEACons if WHOR_c==`R' & male==0
scalar r2 = r(rho)^2
scalar r2_var = (4*r2*((1-r2)^2))/r(N)
matrix r2_fm`R'_female = r2, r2_var
scalar r = r(rho)
scalar n = r(N)
matrix z_fm`R'_female = 0.5 * (ln((1 + r) / (1 - r))) , 1 / sqrt(n - 3)

capture drop diff
gen diff = (DD_fm_IAEACons -  pred_fm_R`R')^2  if WHOR_c==`R' & male==0
summ diff  if WHOR_c==`R' & male==0
matrix rmse_fm`R'_female =  sqrt((r(sum) / r(N)))  , . , . , .


*Age Group
foreach a of numlist 1/5 {

count if WHOR_c==`R' & ageg==`a'
	if r(N) > 10 {
	regress DD_fm_IAEACons pred_fm_R`R'  if WHOR_c==`R' & ageg==`a'
	matrix temp = r(table)'
	matrix calib_fm`R'_ageg`a' = temp[1,1], temp[1,5], temp[1,6]

	constraint 1 pred_fm_R`R' = 1
	cnsreg DD_fm_IAEACons pred_fm_R`R' if WHOR_c==`R' & ageg==`a', constraints(1)
	matrix cil_fm`R'_ageg`a' = temp[2,1], temp[2,5], temp[2,6] 

	qui corr pred_fm_R`R' DD_fm_IAEACons if WHOR_c==`R' & ageg==`a'
	scalar r2 = r(rho)^2
	scalar r2_var = (4*r2*((1-r2)^2))/r(N)
	matrix r2_fm`R'_ageg`a' = r2, r2_var
	scalar r = r(rho)
	scalar n = r(N)
	matrix z_fm`R'_ageg`a' = 0.5 * (ln((1 + r) / (1 - r))) , 1 / sqrt(n - 3)

	capture drop diff
	gen diff = (DD_fm_IAEACons -  pred_fm_R`R')^2  if WHOR_c==`R' & ageg==`a'
	summ diff  if WHOR_c==`R' & ageg==`a'
	matrix rmse_fm`R'_ageg`a' =  sqrt((r(sum) / r(N)))  , . , . , .
	}

	else {
	matrix calib_fm`R'_ageg`a' = .,.,.
	matrix cil_fm`R'_ageg`a' = .,.,.
	matrix r2_fm`R'_ageg`a' = .,.
	matrix z_fm`R'_ageg`a' = .,.
	matrix rmse_fm`R'_ageg`a' = .,.,.,.
	}

}



*Eth Group
foreach a of numlist 1/5 {

count if WHOR_c==`R' & eth_group==`a'
	if r(N) > 10 {
	regress DD_fm_IAEACons pred_fm_R`R'  if WHOR_c==`R' & eth_group==`a'
	matrix temp = r(table)'
	matrix calib_fm`R'_eth`a' = temp[1,1], temp[1,5], temp[1,6]

	constraint 1 pred_fm_R`R' = 1
	cnsreg DD_fm_IAEACons pred_fm_R`R' if WHOR_c==`R' & eth_group==`a', constraints(1)
	matrix cil_fm`R'_eth`a' = temp[2,1], temp[2,5], temp[2,6] 

	qui corr pred_fm_R`R' DD_fm_IAEACons if WHOR_c==`R' & eth_group==`a'
	scalar r2 = r(rho)^2
	scalar r2_var = (4*r2*((1-r2)^2))/r(N)
	matrix r2_fm`R'_eth`a' = r2, r2_var
	scalar r = r(rho)
	scalar n = r(N)
	matrix z_fm`R'_eth`a' = 0.5 * (ln((1 + r) / (1 - r))) , 1 / sqrt(n - 3)

	capture drop diff
	gen diff = (DD_fm_IAEACons -  pred_fm_R`R')^2  if WHOR_c==`R' & eth_group==`a'
	summ diff  if WHOR_c==`R' & eth_group==`a'
	matrix rmse_fm`R'_eth`a' =  sqrt((r(sum) / r(N)))  , . , . , .
	}

	else {
	matrix calib_fm`R'_eth`a' = .,.,.
	matrix cil_fm`R'_eth`a' = .,.,.
	matrix r2_fm`R'_eth`a' = .,.
	matrix z_fm`R'_eth`a' = .,.
	matrix rmse_fm`R'_eth`a' = .,.,.,.
	}

}
}

mat r2_fm_intext_ALL = r2_fm1 \ r2_fm2 \ r2_fm3 \ r2_fm4 \ r2_fm5 \ r2_fm6
mat z_fm_intext_ALL = z_fm1 \ z_fm2 \ z_fm3 \ z_fm4 \ z_fm5 \ z_fm6
mat calib_fm_intext_ALL = calib_fm1 \ calib_fm2 \ calib_fm3 \ calib_fm4 \ calib_fm5 \ calib_fm6
mat cil_fm_intext_ALL = cil_fm1 \ cil_fm2 \ cil_fm3 \ cil_fm4 \ cil_fm5 \ cil_fm6
mat rmse_fm_intext_ALL = rmse_fm1 \ rmse_fm2 \ rmse_fm3 \ rmse_fm4 \ rmse_fm5 \ rmse_fm6


foreach x in male female ageg1 ageg2 ageg3 ageg4 ageg5 eth1 eth2 eth3 eth4 eth5 {
mat r2_fm_intext_ALL_`x' = r2_fm1_`x' \ r2_fm2_`x' \ r2_fm3_`x' \ r2_fm4_`x' \ r2_fm5_`x' \ r2_fm6_`x'
mat z_fm_intext_ALL_`x' = z_fm1_`x' \ z_fm2_`x' \ z_fm3_`x' \ z_fm4_`x' \ z_fm5_`x' \ z_fm6_`x'
mat calib_fm_intext_ALL_`x' = calib_fm1_`x' \ calib_fm2_`x' \ calib_fm3_`x' \ calib_fm4_`x' \ calib_fm5_`x'  \ calib_fm6_`x'
mat cil_fm_intext_ALL_`x' = cil_fm1_`x' \ cil_fm2_`x' \ cil_fm3_`x' \ cil_fm4_`x' \ cil_fm5_`x'  \ cil_fm6_`x'
mat rmse_fm_intext_ALL_`x' = rmse_fm1_`x' \ rmse_fm2_`x' \ rmse_fm3_`x' \ rmse_fm4_`x' \ rmse_fm5_`x'  \ rmse_fm6_`x'
}


**POOLING PERFORMANCE ACROSS OMITTED REGIONS
preserve
clear
set obs 7
gen Region = "African Region" in 1
replace Region = "Eastern Mediterranean Region" in 2
replace Region = "European Region" in 3
replace Region = "Region of the Americas" in 4
replace Region = "South-East Asia Region" in 5
replace Region = "Western Pacific Region" in 6
replace Region = "POOLED" in 7

svmat r2_fm_intext_ALL
meta set r2_fm_intext_ALL1 r2_fm_intext_ALL2 if Region!="POOLED", random(reml) studylabel(Region)
meta summarize if Region!="POOLED", random se(hknapp) predinterval(95)
replace r2_fm_intext_ALL1 = round(r(theta)*100,0.01) if Region=="POOLED"
replace _meta_cil  = round(r(ci_lb)*100,0.01) if Region=="POOLED"
replace _meta_ciu = round(r(ci_ub)*100,0.01) if Region=="POOLED"
meta forestplot, random se(hknapp) esrefline xline(1, lwidth(thin) lcol(gs8)) xlab(0.5(0.1)1.0) xtitle(R{superscript:2})

svmat z_fm_intext_ALL
meta set z_fm_intext_ALL1 z_fm_intext_ALL2 if Region!="POOLED", random(reml) studylabel(Region)
meta summarize if Region!="POOLED", random predinterval(95) transform(corr)
meta forestplot, random esrefline xline(1, lwidth(thin) lcol(gs8))  xtitle(R{superscript:2}) transform(corr)


svmat calib_fm_intext_ALL
meta set calib_fm_intext_ALL1 calib_fm_intext_ALL2 calib_fm_intext_ALL3 if Region!="POOLED", random(reml) studylabel(Region)
meta summarize if Region!="POOLED", random se(hknapp) predinterval(95)
replace calib_fm_intext_ALL1 = r(theta) if Region=="POOLED"
replace calib_fm_intext_ALL2  = r(ci_lb) if Region=="POOLED"
replace calib_fm_intext_ALL3 = r(ci_ub) if Region=="POOLED"
meta forestplot if Region!="POOLED", random se(hknapp) esrefline xline(1, lwidth(thin) lcol(gs8)) xlab(0.5(0.2)1.5) xtitle(Calibration Slope)


svmat cil_fm_intext_ALL
meta set cil_fm_intext_ALL1 cil_fm_intext_ALL2 cil_fm_intext_ALL3 if Region!="POOLED", random(reml) studylabel(Region)
meta summarize if Region!="POOLED", random se(hknapp) predinterval(95)
replace cil_fm_intext_ALL1 = r(theta) if Region=="POOLED"
replace cil_fm_intext_ALL2  = r(ci_lb) if Region=="POOLED"
replace cil_fm_intext_ALL3 = r(ci_ub) if Region=="POOLED"
meta forestplot if Region!="POOLED", random se(hknapp) esrefline xline(0, lwidth(thin) lcol(gs8)) xlab(-14(2)10, angle(vertical)) xtitle(CITL)
restore




*SEX-SPECIFIC PERFORMANCE ACROSS OMITTED REGIONS
preserve
clear
mat r2_fm_intext_ALL_malefem = r2_fm_intext_ALL_female\r2_fm_intext_ALL_male
mat z_fm_intext_ALL_malefem = z_fm_intext_ALL_female\z_fm_intext_ALL_male
mat calib_fm_intext_ALL_malefem = calib_fm_intext_ALL_female\calib_fm_intext_ALL_male
mat cil_fm_intext_ALL_malefem = cil_fm_intext_ALL_female\cil_fm_intext_ALL_male

set obs 12
gen Region = "African Region" in 1/2
replace Region = "Eastern Mediterranean Region" in 3/4
replace Region = "European Region" in 5/6
replace Region = "Region of the Americas" in 7/8
replace Region = "South-East Asia Region" in 9/10
replace Region = "Western Pacific Region" in 11/12

gen Sex = "Male" in 1
	foreach x of numlist 3 5 7 9 11 {
	replace Sex = "Male" in `x'
	}
	replace Sex = "Female" if missing(Sex)

sort Sex Region 
svmat r2_fm_intext_ALL_malefem
svmat z_fm_intext_ALL_malefem
svmat calib_fm_intext_ALL_malefem
svmat cil_fm_intext_ALL_malefem

set obs 18
replace Region = "African Region" in 13
replace Region = "Eastern Mediterranean Region" in 14
replace Region = "European Region" in 15
replace Region = "Region of the Americas" in 16
replace Region = "South-East Asia Region" in 17
replace Region = "Western Pacific Region" in 18
replace Sex = "POOLED" if missing(Sex)

encode Region, gen(region)


foreach x of numlist 1/6 {
meta set r2_fm_intext_ALL_malefem1 r2_fm_intext_ALL_malefem2 if region==`x', random(reml) studylabel(Sex)
meta summarize if region==`x', random predinterval(95)
replace r2_fm_intext_ALL_malefem1 = round(r(theta)*100,0.01) if Sex=="POOLED" & region==`x'
replace _meta_cil  = round(r(ci_lb)*100,0.01) if Sex=="POOLED" & region==`x'
replace _meta_ciu = round(r(ci_ub)*100,0.01) if Sex=="POOLED" & region==`x'
meta forestplot if region==`x', random esrefline xline(1, lwidth(thin) lcol(gs8)) xlab(0.5(0.1)1.0) title(Region `x')

meta set calib_fm_intext_ALL_malefem1 calib_fm_intext_ALL_malefem2 calib_fm_intext_ALL_malefem3  if region==`x', random(reml) studylabel(Sex)
meta summarize if region==`x', random predinterval(95)
replace calib_fm_intext_ALL_malefem1 = r(theta) if Sex=="POOLED" & region==`x'
replace calib_fm_intext_ALL_malefem2  = r(ci_lb) if Sex=="POOLED" & region==`x'
replace calib_fm_intext_ALL_malefem3 = r(ci_ub) if Sex=="POOLED" & region==`x'
meta forestplot  if region==`x', random esrefline xline(1, lwidth(thin) lcol(gs8)) xlab(0.5(0.2)1.5) title(Region `x')

meta set cil_fm_intext_ALL_malefem1 cil_fm_intext_ALL_malefem2 cil_fm_intext_ALL_malefem3  if region==`x', random(reml) studylabel(Sex)
meta summarize  if region==`x', random predinterval(95)
replace cil_fm_intext_ALL_malefem1 = r(theta) if Sex=="POOLED" & region==`x'
replace cil_fm_intext_ALL_malefem2  = r(ci_lb) if Sex=="POOLED" & region==`x'
replace cil_fm_intext_ALL_malefem3 = r(ci_ub) if Sex=="POOLED" & region==`x'
meta forestplot  if region==`x', random esrefline xline(0, lwidth(thin) lcol(gs8)) xlab(-14(2)10, angle(vertical)) title(Region `x')
}
restore




**AGE-SPECIFIC PERFORMANCE ACROSS OMITTED REGIONS
preserve
clear
mat r2_fm_intext_ALL_ages = r2_fm_intext_ALL_ageg1\r2_fm_intext_ALL_ageg2\r2_fm_intext_ALL_ageg3\r2_fm_intext_ALL_ageg4\r2_fm_intext_ALL_ageg5
mat calib_fm_intext_ALL_ages = calib_fm_intext_ALL_ageg1\calib_fm_intext_ALL_ageg2\calib_fm_intext_ALL_ageg3\calib_fm_intext_ALL_ageg4\calib_fm_intext_ALL_ageg5
mat cil_fm_intext_ALL_ages = cil_fm_intext_ALL_ageg1\cil_fm_intext_ALL_ageg2\cil_fm_intext_ALL_ageg3\cil_fm_intext_ALL_ageg4\cil_fm_intext_ALL_ageg5

set obs 30
gen Region = "African Region" in 1/5
replace Region = "Eastern Mediterranean Region" in 6/10
replace Region = "European Region" in 11/15
replace Region = "Region of the Americas" in 16/20
replace Region = "South-East Asia Region" in 21/25
replace Region = "Western Pacific Region" in 26/30

gen Age = ""
	foreach a of numlist 1 6 11 16 21 26 {
	replace Age = "16-24 years" in `a'
local b = `a'+1	
	replace Age = "25-39 years" in `b'
local c = `a'+2	
	replace Age = "40-64 years" in `c'
local d = `a'+3	
	replace Age = "65-79 years" in `d'
local e = `a'+4	
	replace Age = "80+ years" in `e'
	}

sort Age Region 
svmat r2_fm_intext_ALL_ages
svmat calib_fm_intext_ALL_ages
svmat cil_fm_intext_ALL_ages

set obs 36
replace Region = "African Region" in 31
replace Region = "Eastern Mediterranean Region" in 32
replace Region = "European Region" in 33
replace Region = "Region of the Americas" in 34
replace Region = "South-East Asia Region" in 35
replace Region = "Western Pacific Region" in 36
replace Age = "POOLED" if missing(Age)

encode Region, gen(region)

foreach x of numlist 1 3/6 {
meta set r2_fm_intext_ALL_ages1 r2_fm_intext_ALL_ages2 if region==`x', random(reml) studylabel(Age)
meta summarize if region==`x', random  predinterval(95)
replace r2_fm_intext_ALL_ages1 = round(r(theta)*100,0.01) if Age=="POOLED" & region==`x'
replace _meta_cil  = round(r(ci_lb)*100,0.01) if Age=="POOLED" & region==`x'
replace _meta_ciu = round(r(ci_ub)*100,0.01) if Age=="POOLED" & region==`x'
meta forestplot if region==`x', random esrefline xline(1, lwidth(thin) lcol(gs8)) xlab(0.5(0.1)1.0) title(Region `x')

meta set calib_fm_intext_ALL_ages1 calib_fm_intext_ALL_ages2 calib_fm_intext_ALL_ages3  if region==`x', random(reml) studylabel(Age)
meta summarize if region==`x', random  predinterval(95)
replace calib_fm_intext_ALL_ages1 = r(theta) if Age=="POOLED" & region==`x'
replace calib_fm_intext_ALL_ages2  = r(ci_lb) if Age=="POOLED" & region==`x'
replace calib_fm_intext_ALL_ages3 = r(ci_ub) if Age=="POOLED" & region==`x'
meta forestplot  if region==`x', random  esrefline xline(1, lwidth(thin) lcol(gs8)) xlab(0.5(0.2)1.5) title(Region `x')

meta set cil_fm_intext_ALL_ages1 cil_fm_intext_ALL_ages2 cil_fm_intext_ALL_ages3  if region==`x', random(reml) studylabel(Age)
meta summarize  if region==`x', random  predinterval(95)
replace cil_fm_intext_ALL_ages1 = r(theta) if Age=="POOLED" & region==`x'
replace cil_fm_intext_ALL_ages2  = r(ci_lb) if Age=="POOLED" & region==`x'
replace cil_fm_intext_ALL_ages3 = r(ci_ub) if Age=="POOLED" & region==`x'
meta forestplot  if region==`x', random  esrefline xline(0, lwidth(thin) lcol(gs8)) xlab(-14(2)10, angle(vertical)) title(Region `x')
}
restore




**Eth-SPECIFIC PERFORMANCE ACROSS OMITTED REGIONS
preserve
clear
mat r2_fm_intext_ALL_eth = r2_fm_intext_ALL_eth1\r2_fm_intext_ALL_eth2\r2_fm_intext_ALL_eth3\r2_fm_intext_ALL_eth4\r2_fm_intext_ALL_eth5
mat calib_fm_intext_ALL_eth = calib_fm_intext_ALL_eth1\calib_fm_intext_ALL_eth2\calib_fm_intext_ALL_eth3\calib_fm_intext_ALL_eth4\calib_fm_intext_ALL_eth5
mat cil_fm_intext_ALL_eth = cil_fm_intext_ALL_eth1\cil_fm_intext_ALL_eth2\cil_fm_intext_ALL_eth3\cil_fm_intext_ALL_eth4\cil_fm_intext_ALL_eth5

set obs 30
gen Region = "African Region" in 1/5
replace Region = "Eastern Mediterranean Region" in 6/10
replace Region = "European Region" in 11/15
replace Region = "Region of the Americas" in 16/20
replace Region = "South-East Asia Region" in 21/25
replace Region = "Western Pacific Region" in 26/30

gen Eth = .
	foreach a of numlist 1 6 11 16 21 26 {
	replace Eth = 1 in `a'
local b = `a'+1	
	replace Eth = 2 in `b'
local c = `a'+2	
	replace Eth = 3 in `c'
local d = `a'+3	
	replace Eth = 4 in `d'
local e = `a'+4	
	replace Eth = 5 in `e'
	}

sort Eth Region 
svmat r2_fm_intext_ALL_eth
svmat calib_fm_intext_ALL_eth
svmat cil_fm_intext_ALL_eth

set obs 36
replace Region = "African Region" in 31
replace Region = "Eastern Mediterranean Region" in 32
replace Region = "European Region" in 33
replace Region = "Region of the Americas" in 34
replace Region = "South-East Asia Region" in 35
replace Region = "Western Pacific Region" in 36


label define Eth 1 "Whites" 2 "Blacks" 3 "" 4 "" 5 "Other"
label values Eth Eth


encode Region, gen(region)


gen Eth2 = "White" if Eth==1
replace Eth2 = "Blacks" if Eth==2
replace Eth2 = "South Asian" if Eth==3
replace Eth2 = "Other Asian" if Eth==4
replace Eth2 = "Other" if Eth==5


foreach x of numlist 1/6 {
meta set r2_fm_intext_ALL_eth1 r2_fm_intext_ALL_eth2 if region==`x', random(reml) studylabel(Eth2)
meta summarize if region==`x', random  predinterval(95)
replace r2_fm_intext_ALL_eth1 = round(r(theta)*100,0.01) if missing(Eth) & region==`x'
replace _meta_cil  = round(r(ci_lb)*100,0.01) if missing(Eth) & region==`x'
replace _meta_ciu = round(r(ci_ub)*100,0.01) if missing(Eth) & region==`x'
meta forestplot if region==`x', random esrefline xline(1, lwidth(thin) lcol(gs8)) xlab(0.5(0.1)1.0) title(Region `x') 

meta set calib_fm_intext_ALL_eth1 calib_fm_intext_ALL_eth2 calib_fm_intext_ALL_eth3  if region==`x', random(reml) studylabel(Eth2)
meta summarize if region==`x', random  predinterval(95)
replace calib_fm_intext_ALL_eth1 = r(theta) if missing(Eth) & region==`x'
replace calib_fm_intext_ALL_eth2  = r(ci_lb) if missing(Eth) & region==`x'
replace calib_fm_intext_ALL_eth3 = r(ci_ub) if missing(Eth) & region==`x'
meta forestplot  if region==`x', random  esrefline xline(1, lwidth(thin) lcol(gs8)) xlab(0.5(0.2)1.5) title(Region `x')

meta set cil_fm_intext_ALL_eth1 cil_fm_intext_ALL_eth2 cil_fm_intext_ALL_eth3  if region==`x', random(reml) studylabel(Eth2)
meta summarize  if region==`x', random  predinterval(95)
replace cil_fm_intext_ALL_eth1 = r(theta) if missing(Eth) & region==`x'
replace cil_fm_intext_ALL_eth2  = r(ci_lb) if missing(Eth) & region==`x'
replace cil_fm_intext_ALL_eth3 = r(ci_ub) if missing(Eth) & region==`x'
meta forestplot  if region==`x', random  esrefline xline(0, lwidth(thin) lcol(gs8)) xlab(-14(2)10, angle(vertical)) title(Region `x')
}
restore


