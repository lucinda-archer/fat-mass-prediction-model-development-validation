clear
use "DLW database v4.0 abbreviated 20240322 hudda v2_ADULTS.dta", clear
drop if DBID==10893 & IND_ID == 10151
destring BMI wt_av, replace
append using "DLW database v4.0 abbreviated 20240322 hudda v2_CHILD.dta"
keep if Age>=16
merge m:1 CountryISO using "Countries data"
replace Country="Scotland" if CountryISO=="SCT"
replace WHORegion = "WHO European Region"  if CountryISO=="SCT" | CountryISO=="SWE" | CountryISO=="SWE " 
drop if _merge==2
drop _merge

gen Region = "African" if WHORegion=="WHO African Region"
replace Region = "European" if WHORegion=="WHO European Region"
replace Region = "Americas" if WHORegion=="WHO Region of the Americas"
replace Region = "WesternPacific" if WHORegion=="WHO Western Pacific Region"
replace Region = "EastMed+SE.Asia" if WHORegion=="WHO Eastern Mediterranean Region" | WHORegion=="WHO South-East Asia Region"

*Dropping those with recorded disease
keep if Health == "H" | Health == "H " | Health == " H"| Health == " H " | Health=="NA" | Health ==""
*Dropping those who are confirmed pregnant
keep if rep_statF=="" | rep_statF=="NA" | rep_statF=="PM" | rep_statF=="PP" | rep_statF=="RANR" | rep_statF=="RANR " | rep_statF=="Reproductive age"
*Dropping professional/ semi-professional athletes
keep if Athlete==0

*Generating new /formatting existing variables
gen weight = wt_av 
gen htm = ht/100
gen male = 1 if Sex=="M" | Sex=="M "
replace male = 0 if Sex == "F"
tab male Sex
rename LBM DD_ffm_IAEACons
rename FTM DD_fm_IAEACons
destring FMkg FFMkg FM, replace
rename FFMkg Other_ffm
rename FMkg Other_fm
rename FM Other_fm_perc
rename FFMFMmethod Other_method
rename machine Other_machine
gen eth_group = 5 if eth1=="" | eth1=="H" | eth1=="NA"
replace eth_group = 1 if eth1=="C" | (eth1=="O" & eth2=="Italian") | (eth2=="Scottish, Irish")
replace eth_group = 2 if eth1== "A" | eth1=="AA" | (eth1=="O" & eth2=="Jamaican")
replace eth_group = 3 if (eth1=="AS" & missing(eth2)) | eth2 =="Indian"
replace eth_group = 4 if eth1=="AS" & missing(eth_group)
replace eth_group = 5 if missing(eth_group)
label define eth 1 "White" 2 "Blacks" 3 "South Asian" 4 "Other Asian" 5 "Other/Unknown"
lab values eth_group eth
gen BA = 1 if eth_group==2
replace BA = 0 if missing(BA)
gen SA = 1 if eth_group==3
replace SA = 0 if missing(SA)
gen AO = 1 if eth_group==4
replace AO = 0 if missing(AO)
gen Other = 1 if eth_group==5
replace Other = 0 if missing(Other)
encode CountryISO, gen(Country_coded)
encode WHORegion, gen(WHOR_coded)
destring  DD_ffm_IAEACons DD_fm_IAEACons , replace
gen lnDD_ffm_IAEACons = ln(DD_ffm_IAEACons)
gen lnDD_fm_IAEACons = ln(DD_fm_IAEACons)
gen ageg = 1 if Age>=16 & Age <25
replace ageg = 2 if Age >=25 & Age < 40
replace ageg = 3 if Age >=40 & Age < 65
replace ageg = 4 if Age >=65 & Age < 80
replace ageg = 5 if Age >=80 & !missing(Age)

*Checking for implausible Values of lnffm
summ lnDD_ffm_IAEACons,d
count if lnDD_ffm_IAEACons < (r(mean) - 10* r(sd))
qui summ lnDD_ffm_IAEACons,d
	qui summ lnDD_ffm_IAEACons,d
	drop if (lnDD_ffm_IAEACons < (r(mean) - 10* r(sd)))
count if lnDD_ffm_IAEACons > (r(mean) + 10* r(sd))

*Histograms in dev data by sex
hist DD_ffm_IAEACons, by(Sex) freq xtitle(Fat Free Mass (kg)) aspect(1)
hist lnDD_ffm_IAEACons, by(Sex) freq xtitle(ln(Fat Free Mass)) aspect(1)
hist DD_fm_IAEACons, by(Sex) freq xtitle(Fat Mass (kg)) aspect(1)

*Summary Stats of Age, ht, wt, fm, ffm, male, Ethnicity by WHO Regions and overall
foreach x of varlist Age htm weight DD_ffm_IAEACons DD_fm_IAEACons {
	foreach y of numlist 1/6 {
	qui summ `x' if WHOR_coded == `y',d

	mat `x'_R`y' = r(p50), r(p25), r(p75)
	}

qui summ `x' ,d
mat `x'_overall = r(p50), r(p25), r(p75)

mat `x' = `x'_R1 , `x'_R2 , `x'_R3 , `x'_R4 , `x'_R5 , `x'_R6 , `x'_overall
}

foreach y of numlist 1/6 {
count if male==1 & WHOR_coded == `y'
scalar male_n = r(N)
count if WHOR_coded == `y'
scalar n = r(N)
mat male_R`y' = male_n, male_n/n *100,.
}
count if male==1
scalar male_n = r(N)
count
scalar n = r(N)
mat male_overall = male_n, male_n/n *100,.
mat male = male_R1 , male_R2 , male_R3 , male_R4 , male_R5 , male_R6 , male_overall


foreach z of numlist 1/5 {
	foreach y of numlist 1/6 {
	count if eth_group==`z' & WHOR_coded == `y'
	scalar eth_n = r(N)
	count if WHOR_coded == `y'
	scalar n = r(N)
	mat eth`z'_R`y' = eth_n, eth_n/n *100,.
	}
	count if eth_group==`z'
	scalar eth = r(N)
	count
	scalar n = r(N)
	mat eth`z'_overall = eth, eth/n *100,.
	
	mat eth`z' = eth`z'_R1, eth`z'_R2, eth`z'_R3, eth`z'_R4, eth`z'_R5,eth`z'_R6, eth`z'_overall	
}
mat ethnic = eth1\eth2\eth3\eth4\eth5
mat Tab1 = Age \ htm \ weight \ DD_ffm_IAEACons \ DD_fm_IAEACons \ male \ ethnic


**Model Building
*Generate Dummy Variables / Interaction Terms
gen WE = 1 if eth_group==1
replace WE = 0 if eth_group!=1
gen htm_BA = htm*BA
gen htm_SA = htm*SA
gen htm_AO = htm*AO
gen htm_Oth = htm*Other
gen htm_WE = htm*WE
gen weight_BA = weight*BA
gen weight_SA = weight*SA
gen weight_AO = weight*AO
gen weight_Oth = weight*Other
gen weight_WE = weight*WE

gen wt_male = weight*male
gen ht_male = htm*male
gen Age_male = Age*male

gen wt_Age = weight*Age
gen ht_Age = htm*Age	

gen WHO_Afr = 1 if WHOR_c==1
replace WHO_Afr = 0 if WHOR_c!=1
gen WHO_EMed = 1 if WHOR_c==2
replace WHO_EMed = 0 if WHOR_c!=2
gen WHO_Eur = 1 if WHOR_c==3
replace WHO_Eur = 0 if WHOR_c!=3
gen WHO_Amer = 1 if WHOR_c==4
replace WHO_Amer = 0 if WHOR_c!=4
gen WHO_SEAsia = 1 if WHOR_c==5
replace WHO_SEAsia = 0 if WHOR_c!=5
gen WHO_WPac = 1 if WHOR_c==6
replace WHO_WPac = 0 if WHOR_c!=6
	
foreach x of varlist WHO_Afr WHO_EMed WHO_Eur WHO_Amer WHO_SEAsia WHO_WPac {
	gen wt_`x' = weight*`x'
	gen ht_`x' = htm*`x'
	gen Age_`x' = Age*`x'
}
gen ht_wt = htm*weight

capture drop pred_lnffm_sexcomb 
capture drop pred_ffm_sexcomb 
capture drop pred_fm_sexcomb 

xi: mfp, alpha(0.05) select(0.157) center(no) all: regress lnDD_ffm_IAEACons htm weight Age BA SA AO Other i.male (WHO_Afr WHO_EMed WHO_Eur WHO_SEAsia WHO_WPac) ht_Age	ht_male		wt_Age	wt_male
predict pred_lnffm_sexcomb, xb

capture gen pred_ffm_sexcomb = exp(pred_lnffm_sexcomb)
capture gen pred_fm_sexcomb = weight - exp(pred_lnffm_sexcomb)
scalar shrink = (-2*(e(ll_0)-e(ll)) - 10) / (-2*(e(ll_0)-e(ll)))
di shrink


**PERFORMANCE
regress DD_fm_IAEACons pred_fm
matrix temp = r(table)'
matrix slope_fm = temp[1,1], temp[1,5], temp[1,6]

constraint 1 pred_fm = 1
cnsreg DD_fm_IAEACons pred_fm, constraints(1)
matrix temp = r(table)'
matrix cil_fm = temp[2,1], temp[2,5], temp[2,6]

pmcalplot pred_fm DD_fm_IAEACons, cont nohist xtitle(Predicted FM (kg)) ytitle(Deuterium Dilution Assessed FM (kg)) legend(off) range(0 150)


*R2
corr DD_fm_IAEACons pred_fm 
scalar r2_DD_fm_IAEACons = r(rho)^2 
scalar r2_DD_fm_IAEACons_var = (4*r2_DD_fm_IAEACons*((1-r2_DD_fm_IAEACons)^2))/r(N)
matrix r2_fm = r2_DD_fm_IAEACons*100 , (r2_DD_fm_IAEACons  - 1.96*(sqrt(r2_DD_fm_IAEACons_var)))*100 , (r2_DD_fm_IAEACons  + 1.96*(sqrt(r2_DD_fm_IAEACons_var)))*100 


*RMSE
capture drop diff
gen diff = (DD_fm_IAEACons -  pred_fm)^2
summ diff
di sqrt((r(sum) / r(N)))
matrix rmse_fm =  sqrt((r(sum) / r(N)))

mat stats_fm = r2_fm , slope_fm , cil_fm , rmse_fm
mat colnames stats_fm = r2_fm  r2_fm_low r2_fm_up slope_fm  slope_fm_low slope_fm_up cil_fm  cil_fm_low cil_fm_up rmse_fm
