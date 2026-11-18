*load formatted data
use "${data}\formatted_NHANES", clear


********************************************
**# Define region 
/* Select as appropriate from below:
	1 = Africa
	2 = Eastern Mediterranean
	3 = European
	4 = Americas 
	5 = South-East Asia
	6 = Western Pacific
	
 e.g., KNHANES => South Korea => Western Pacific => region = 6 
*/

gen region = 6 // 

*label for regions, per WHO regions above
label define region_short  1 "Africa" ///
						   2 "Eastern Mediterranean" ///
						   3 "European" ///
						   4 "Americas" ///
						   5 "South-East Asia" ///
						   6 "Western Pacific"
label values region region_short

********************************************
**# Generate predicted values


*gen prediction variable
/*
gen pred_lnffm = 76.8548 + ///
    0.5392*heightm + ///
    (-0.0322)*(weight/100)^(-2) + ///
    0.4578*(weight/100) + ///
    (-0.0011)*age + ///
    0.0100*[eth_group==2] + ///      Black
    0.0113*[eth_group==3] + ///      South Asian
    (-0.0142)*[eth_group==4] + ///   Other Asian
    0.0184*[eth_group==5] + ///      Other/Unknown
    (-73.5434)*male + ///            Male (Female = reference)
    0.0008*[region==1] + ///         Africa
    (-0.0579)*[region==2] + ///      Eastern Mediterranean
    0.0140*[region==3] + ///         European
    (-0.0150)*[region==5] + ///      South-East Asia
    (-0.0304)*[region==6] + ///      Western Pacific
    0.0361*(((heightm*age)/100)^(-2)) + ///
    0.0262*((((heightm*age)/100)^(-2))*ln((heightm*age)/100)) + ///
    (-0.0771)*(heightm*male) + ///
    (-0.3164)*((((weight*male)/*+0.001816*/)/100)^(-0.5)) + ///
    (-0.0979)*((((weight*male)/*+0.001816*/)/100)^2)
	
*/
gen pred_lnffm = .

replace pred_lnffm = (0.4621*heightm - 322*weight^(-2) - 3.164*weight^(-0.5) + 0.004578*weight - 0.00000979*weight^2 - 0.0011*age + 0.000006*weight*age + 0.01*[eth_group==2] + 0.0113*[eth_group==3] - 0.0142*[eth_group==4] + 0.0184*[eth_group==5] + 361*(heightm*age)^(-2) + 0.0262*(((heightm*age)/100)^(-2))*ln((heightm*age)/100) + 3.3114) if male == 1

replace pred_lnffm = (0.5392*heightm - 322*weight^(-2) + 0.004578*weight - 0.0011*age + 0.000006*weight*age + 0.01*[eth_group==2] + 0.0113*[eth_group==3] - 0.0142*[eth_group==4] + 0.0184*[eth_group==5] + 361*(heightm*age)^(-2) + 0.0262*(((heightm*age)/100)^(-2))*ln((heightm*age)/100) + 2.6037) if male == 0
	
gen pred_ffm = exp(pred_lnffm)
gen pred_fm = weight - pred_ffm 

*save formatted data with predictions
save "${data}\formatted_NHANES", replace

*******************************
**# Update

*save formatted data with predictions
use "${data}\formatted_NHANES", clear

gen lnffm_dxa = ln(ffm_dxa)
gen upd_pred_lnffm = .
gen fmscale_upd_pred_fm = .


regress lnffm_dxa pred_lnffm if male == 0
	local female_slope = e(b)[1,1]
	local female_int = e(b)[1,2]
	replace upd_pred_lnffm = `female_slope'*pred_lnffm + `female_int' if male == 0
	
	
regress lnffm_dxa pred_lnffm if male == 1
	local male_slope = e(b)[1,1]
	local male_int = e(b)[1,2]
	replace upd_pred_lnffm = `male_slope'*pred_lnffm + `male_int' if male == 1
	

/*	**Decsion to update on the lnffm scale so no longer needed**

regress fm_dxa pred_fm if male == 0
	local female_slope = e(b)[1,1]
	local female_int = e(b)[1,2]
	replace fmscale_upd_pred_fm = `female_slope'*pred_fm + `female_int' if male == 0
	
regress fm_dxa pred_fm if male == 1
	local male_slope = e(b)[1,1]
	local male_int = e(b)[1,2]
	replace fmscale_upd_pred_fm = `male_slope'*pred_fm + `male_int' if male == 1
*/

gen upd_pred_fm = weight - exp(upd_pred_lnffm)
gen orig_pred_fm = pred_fm


tempname val_stats
tempfile val_stats_f
postfile `val_stats' str20 subgroup str20 group ///
					str20 model N ///
					cslope cslope_se cslope_lb cslope_ub ///
					citl citl_se citl_lb citl_ub ///
					rmse ///
					r2a r2a_se r2a_lb r2a_ub ///
					z_r2a z_r2a_se ///
					using `val_stats_f', replace
					

					
foreach model in orig upd /*fmscale_upd*/ {
	
	count
	local N = r(N)
	
	*Calibration slope
	regress fm_dxa `model'_pred_fm
					mat cslope_est = r(table)
					local cslope = cslope_est[1,1]
					local cslope_se = cslope_est[2,1]
					local cslope_lb = cslope_est[5,1]
					local cslope_ub = cslope_est[6,1]	

	*R-squared
	local r2 = e(r2)
	local n = e(N)
	local p = 19 // M[1,`model'+1]
	local r2a = `r2' // (((`n'-1)*`r2')-`p')/(`n'-`p'-1) - External validation

	local r2a_var = 4*(`r2a')*(1-`r2a')^2/`n'
	local r2a_lb = `r2a' - 1.96*sqrt(`r2a_var')
	local r2a_ub = `r2a' + 1.96*sqrt(`r2a_var')
	local r2a_se = sqrt(`r2a_var')/`n'

	local z_r2a = atanh(sqrt(`r2a'))
	local z_r2a_se = 1/(`n'-3)

		
	*Calculate calibration-in-the-large
	constraint 1 `model'_pred_fm=1
	cnsreg  fm_dxa `model'_pred_fm, c(1)
		mat citl_est = r(table)
		local citl = citl_est[1,2]
		local citl_se = citl_est[2,2]
		local citl_lb = citl_est[5,2]
		local citl_ub = citl_est[6,2]
		
	*Calculate RMSE
	gen sqerr_`model' = (fm_dxa - `model'_pred_fm)^2 //tempvar?
	sum sqerr_`model'
	local rmse = sqrt(r(mean))

	drop sqerr_`model'

				*Post results to file
				post `val_stats' ("all") ("all") ("`model'") (`N') ///
					(`cslope') (`cslope_se') (`cslope_lb') (`cslope_ub') ///
					(`citl') (`citl_se') (`citl_lb') (`citl_ub') ///
					(`rmse')  ///
					(`r2a') (`r2a_se') (`r2a_lb') (`r2a_ub') ///
					(`z_r2a') (`z_r2a_se') 



	foreach subgroup in male age_group eth_group {
		
		levelsof `subgroup', local(levels)
		
		foreach group of local levels {
						
			**CONTINUOUS**	
			di "`subgroup': `group'"
			
			count  if `subgroup' == `group' & `model'_pred_fm != .
			
				if r(N) < 10 {
												
				*Post results to file
					post `val_stats' ("`subgroup'") ("`group'") ("`model'") (0) ///
					(.) (.) (.) (.) ///
					(.) (.) (.) (.) ///
					(.)  ///
					(.) (.) (.) (.) ///
					(.) (.)
						
				}
				
				else if r(N) >= 0 {	
					
				local N = r(N)
				
					*first calculate the calibration slope
				regress fm_dxa `model'_pred_fm if `subgroup' == `group'
					mat cslope_est = r(table)
					local cslope = cslope_est[1,1]
					local cslope_se = cslope_est[2,1]
					local cslope_lb = cslope_est[5,1]
					local cslope_ub = cslope_est[6,1]
				
				*calculate R2				
				local r2 = e(r2)
				local n = e(N)
				local p = 19 // M[1,`model'+1]
				local r2a = `r2' // (((`n'-1)*`r2')-`p')/(`n'-`p'-1) - External validation
				
				local r2a_var = 4*(`r2a')*(1-`r2a')^2/`n'
				local r2a_lb = `r2a' - 1.96*sqrt(`r2a_var')
				local r2a_ub = `r2a' + 1.96*sqrt(`r2a_var')
				local r2a_se = sqrt(`r2a_var')/`n'
				
				local z_r2a = atanh(sqrt(`r2a'))
				local z_r2a_se = 1/(`n'-3)
				
					
				*Calculate calibration-in-the-large
				constraint 1 `model'_pred_fm=1
				cnsreg  fm_dxa `model'_pred_fm if `subgroup' == `group', c(1)
					mat citl_est = r(table)
					local citl = citl_est[1,2]
					local citl_se = citl_est[2,2]
					local citl_lb = citl_est[5,2]
					local citl_ub = citl_est[6,2]
					
				*Calculate RMSE
				gen sqerr_`model' = (fm_dxa - `model'_pred_fm)^2 //tempvar?
				sum sqerr_`model' if `subgroup' == `group'
				local rmse = sqrt(r(mean))
				
				drop sqerr_`model'
					
				*Post results to file
				post `val_stats' ("`subgroup'") ("`group'") ("`model'") (`N') ///
					(`cslope') (`cslope_se') (`cslope_lb') (`cslope_ub') ///
					(`citl') (`citl_se') (`citl_lb') (`citl_ub') ///
					(`rmse')  ///
					(`r2a') (`r2a_se') (`r2a_lb') (`r2a_ub') ///
					(`z_r2a') (`z_r2a_se') 
					
				}
					
			}
		}
	
}

postclose `val_stats'
use `val_stats_f', clear
save "${data}\Validation_stats_update", replace

use "${data}\Validation_stats_update", clear

*formatted string variables (3dp)
foreach s in N cslope citl r2a rmse {

    local fmt = cond("`s'" == "N", "%9.0f", "%9.3f")

    capture confirm variable `s'_lb
    if !_rc {
        gen str str_`s' = ///
            string(`s',    "`fmt'") + " (" + ///
            string(`s'_lb, "`fmt'") + " to " + ///
            string(`s'_ub, "`fmt'") + ")"
    }
    else {
        gen str str_`s' = string(`s', "`fmt'")
    }

    replace str_`s' = "" if missing(`s')
}

*display labels for subgroup
gen subgroup_lbl = ""
replace subgroup_lbl = "Overall"   if subgroup == "all"
replace subgroup_lbl = "Ethnicity" if subgroup == "eth_group"
replace subgroup_lbl = "Sex"       if subgroup == "male"
replace subgroup_lbl = "Age group" if subgroup == "age_group"

*display labels for group, conditional on subgroup
gen group_lbl = ""
replace group_lbl = "Overall"       if subgroup == "all"
replace group_lbl = "White"         if subgroup == "eth_group" & group == "1"
replace group_lbl = "Black"         if subgroup == "eth_group" & group == "2"
replace group_lbl = "South Asian"   if subgroup == "eth_group" & group == "3"
replace group_lbl = "Other Asian"   if subgroup == "eth_group" & group == "4"
replace group_lbl = "Other/Unknown" if subgroup == "eth_group" & group == "5"
replace group_lbl = "Female"        if subgroup == "male" & group == "0"
replace group_lbl = "Male"          if subgroup == "male" & group == "1"
replace group_lbl = "16-24"         if subgroup == "age_group" & group == "0"
replace group_lbl = "25-39"         if subgroup == "age_group" & group == "1"
replace group_lbl = "40-64"         if subgroup == "age_group" & group == "2"
replace group_lbl = "65-79"         if subgroup == "age_group" & group == "3"
replace group_lbl = "80+"           if subgroup == "age_group" & group == "4"


*keep tabel vars

foreach do in keep order {
	`do' subgroup_lbl group_lbl model ///
     str_N str_cslope str_citl str_r2a str_rmse
}



	
********************************************
**# Calibration plots

use "${data}\formatted_NHANES", clear

*overall
pmcalplot pred_fm fm_dxa, ///
			cont nostatistics nohist legend(off) ///
			range(0 100) xlabel(0(20)100) ylabel(0(20)100) ///
			xtitle(Predicted FM (kg)) ytitle(DXA assessed FM (kg)) ///
			title("") // UPDATE TO REGION/DATASET NAME, AS APPROPRIATE  
			
	graph export "${plots_file}\calplot_all.png", replace
	*graph export "${plots_file}\calplot_all", replace

*calplots over subgroups
	** exclude small subgroups, <10 observations
foreach subgroup in male age_group eth_group {
	
	levelsof `subgroup', local(levels)
	foreach group of local levels {
		
		count  if `subgroup' == `group' & pred_fm != .
				
		if r(N) < 10 {
				
			drop if `subgroup' == `group'
				
		}
	}	

	pmcalplot pred_fm fm_dxa, ///
			cont nostatistics nohist legend(off) ///
			range(0 100) xlabel(0(20)100) ylabel(0(20)100) ///
			xtitle(Predicted FM (kg)) ytitle(DXA assessed FM (kg)) ///
			by(`subgroup')   
			
	graph export "${plots_file}\calplot_`subgroup'.png", replace
	*graph export "${plots_file}\calplot_`subgroup'", replace
	
}

**update model cal plots
gen lnffm_dxa = ln(ffm_dxa)
gen upd_pred_lnffm = .

regress lnffm_dxa pred_lnffm if male == 0
	local female_slope = e(b)[1,1]
	local female_int = e(b)[1,2]
	replace upd_pred_lnffm = `female_slope'*pred_lnffm + `female_int' if male == 0
	
	
regress lnffm_dxa pred_lnffm if male == 1
	local male_slope = e(b)[1,1]
	local male_int = e(b)[1,2]
	replace upd_pred_lnffm = `male_slope'*pred_lnffm + `male_int' if male == 1
	
gen upd_pred_fm = weight - exp(upd_pred_lnffm)

*overall
pmcalplot upd_pred_fm fm_dxa, ///
			cont nostatistics nohist legend(off) ///
			range(0 100) xlabel(0(20)100) ylabel(0(20)100) ///
			xtitle(Predicted FM (kg)) ytitle(DXA assessed FM (kg)) ///
			title("") scatteropts(mcolor(navy)) // UPDATE TO REGION/DATASET NAME, AS APPROPRIATE  
			
	graph export "${plots_file}\upd_calplot_all.png", replace
	*graph export "${plots_file}\calplot_all", replace

*calplots over subgroups
	** exclude small subgroups, <10 observations
foreach subgroup in male age_group eth_group {
	
	levelsof `subgroup', local(levels)
	foreach group of local levels {
		
		count  if `subgroup' == `group' & pred_fm != .
				
		if r(N) < 10 {
				
			drop if `subgroup' == `group'
				
		}
	}	

	pmcalplot upd_pred_fm fm_dxa, ///
			cont nostatistics nohist legend(off) ///
			range(0 100) xlabel(0(20)100) ylabel(0(20)100) ///
			xtitle(Predicted FM (kg)) ytitle(DXA assessed FM (kg)) ///
			by(`subgroup') scatteropts(mcolor(navy))  
			
	graph export "${plots_file}\upd_calplot_`subgroup'.png", replace
	*graph export "${plots_file}\calplot_`subgroup'", replace
	
}

