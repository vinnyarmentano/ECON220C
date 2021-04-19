*************************************************
** Title: ECON220C Problem Set 1.5
** Author: Vincent Armentano
** Contact: varmenta@ucsd.edu
** Date Created: 4/12/2021
** Last Modified: 
*************************************************

** Purpose:
	{
	/*
	This Dofile contains the code for ECON220C Problem Set 1, Problems 5 & 6.
	*/
	}
*****

** Setup
	{
	** a. Setting Directory
		cd "C:\Users\Vincent\OneDrive - UC San Diego\B_CourseworkAdmin\Year_1\3_Spring2021\ECON220C\PSets\HW1"

	** b. User-Written Programs
		*ssc install blindschemes
		
	** c. RNG Seed
		set seed 462759
	}
*****

** Problem 5
	{
	** 1. Creating Dataset
	qui	{
		** A. Establishing i Structure
			clear
			set obs 500
			gen id = _n
		
		** B. Establishing t Structure
			expand 5
			sort id
			egen t = seq(), f(1) t(5)

		** C. Affirming Success
			isid id t
		}
	*****

	** 3. Parts A-C Simulations
	qui	{
		** A. Establishing Matrix to hold Simulation results
			mat five_a = J(1000,3,.)
			mat coln five_a = bhat sigma_hat sigma_tilde
		
		** B. Establishing sim loop
			qui forval sim = 1/1000 {
				if mod(`sim',100)==0{
					noi di "	Iteration `sim'..."
					}
				preserve
				** i. Creating RHS Variables
					gen x = rnormal(0,1)
					gen u = .
					d,s
					loc obs = `r(N)'
					qui forval i = 1/`obs' {
						loc sd = (abs(x[`i']))
						replace u = rnormal(0,`sd') in `i'
						}
				
				** ii. Estimating Mean Terms
					egen x_bar = mean(x), by(id)
					egen u_bar = mean(u), by(id)
				
				** iii. Creating Bhat Denominator Term
					{
					** a. Creating Squared difference term 
						gen dm_sqdiff_x = (x-x_bar)^2
						
					** b. Executing T & N Summations
						egen denom_sum = total(dm_sqdiff_x)
					}
				*****
				
				** iv. Creating Bhat Numerator Term
					{
					** a. Creating Differences Product term
						gen dm_product_xu = (x-x_bar)*(u-u_bar)
					
					** b. Executing T & N Summations
						egen numer_sum = total(dm_product_xu)
					}
				*****
				
				** v. Estimating Bhat using beta=1
					gen bhat = 1 + numer_sum / denom_sum
				
				** vi. Creating LHS Y given beta=1
					gen y = 1*x + u
					
				** vii. Creating Residual term now & mean across t
					gen u_hat = y - bhat*x
					egen u_hat_bar = mean(u_hat), by(id)
				
				** viii. Creating S_xx var 
					** (yes this is memory intensive but pedagogically nice)
						gen s_xx = denom_sum
				
				** ix. Creating Sigma hat
					{
					** a. Creating Inner product
						gen sighat_inner = (x-x_bar)*(u_hat-u_hat_bar)
					
					** b. Executing T summation & square 500 numbers made from 5 each
						egen sighat_tsum = total(sighat_inner), by(id)
						gen sighat_tsum_sq = sighat_tsum^2 if t==1
						** we need to be careful not to double count
					
					** c. Executing N summation, of just the 500 numbers above
						egen sighat_nsum = total(sighat_tsum_sq)
					
					** d. Multiply with 1/(s_xx^2) for variance
						gen var_hat = (1/s_xx^2)*sighat_nsum
					
					** e. Obtaining SD
						gen sd_hat = sqrt(var_hat)
					}
				*****
				
				** x. Creating Sigma Tilde
					{
					** a. Creating product of squared differences terms
						gen sqdiff_prod_xu_hat = ((x-x_bar)^2)*((u_hat-u_hat_bar)^2)
					
					** b. Executing T & N Summations
						egen sigtilde_sum = total(sqdiff_prod_xu_hat)
						
					** c. Multiply with 1/(s_xx^2) for variance
						gen var_tilde = (1/s_xx^2)*sigtilde_sum
					
					** d. Obtaining SD
						gen sd_tilde = sqrt(var_tilde)
					}
				*****
				
				** xi. Robustness checks
					{
					** a. We're going to compare to Stata's Estimates
						reg y x, r noconstant
						mat raw = r(table)
					
					** b. Comparing Stata's Beta to ours
						sum bhat
						loc beta = `r(mean)'
						if abs(`beta' - raw[1,1]) >= .05 {
							noi di "	- Robustnes Flag, Beta difference  `=abs(`beta'-raw[1,1])'"
							}
					
					** c. Sigma Hat
						sum sd_hat
						loc sd_hat = `r(mean)'
						if abs(`sd_hat' - raw[2,1]) >= .05 {
							noi di "	- Robustnes Flag, Sigma Hat difference  `=abs(`sd_hat'-raw[2,1])'"
							}
					
					** c. Sigma Tidle
					** No equivalent Stata command found
					/*	ren y x
						mat raw = r(table)
						sum sd_tilde
						loc sd_tilde = `r(mean)'
						if abs(`sd_tilde'-raw[1,2])>=.01 {
							noi di "	- Robustnes Flag, Sigma Hat difference  `=abs(`sd_tilde'-raw[1,2])'"
							}
					*/
					}
				*****
				
				** xii. Recording Results
					forval c = 1/3 {
						sum `=word("bhat sd_hat sd_tilde",`c')'
						mat five_a[`sim',`c'] = round(`r(mean)',.00001)
						}
				
				** xiii. Preparing for next iteration
					restore
			} // End Sim Loop
		
		** C. Establishing Matrix for this set of simulations
			{
			** i. Empty matrix
				mat five_t5 = J(8,1,.)
				
			** b. Adjusting row and column headers
				mat rown five_t5 =	std_bhat	///
									bias_hat	bias_tilde	///
									std_hat		std_tilde	///
									rmse_hat	rmse_tilde	///
									rmse_ratio
				mat coln five_t5 = T_5
			}
		*****
		
		** D. Working with Simulation Results
			{
			** i. Moving Results to Memory
				clear 
				svmat five_a, n(col)
			
			** ii. 5b, Bhat SD
				qui sum bhat
				noi di "Beta Hat has an SD of `=round(`r(sd)',.0001)'"
			
			** iii. Recording in matrix, as scalar and as variable
				mat five_t5[1,1] = round(`r(sd)',.0001)
				sca def bhat_sd = `r(sd)'
				egen sd_bhat = sd(bhat)
		
			** iv. Bias, E(estimator)-true
				qui foreach suff in hat tilde {
					sum sigma_`suff'
					loc bias = `r(mean)'-`=scalar(bhat_sd)'
					noi di "Sigma `suff' has bias of `=round(`bias',.0001)'"
					mat five_t5[`=cond("`suff'"=="hat",2,3)',1] = round(`bias',.0001)
					}
			
			** v. Std
				qui foreach suff in hat tilde {
					sum sigma_`suff'
					noi di "Sigma `suff' has an SD of `=round(`r(sd)',.0001)'"
					mat five_t5[`=cond("`suff'"=="hat",4,5)',1] = round(`r(sd)',.0001)
					}
			
			** vi. RMSE
				qui foreach suff in hat tilde {
					** a. Creating Sq difference
						gen sqdiff_`suff' = (sigma_`suff'-sd_bhat)^2
					
					** b. Creating Sum
						egen sum_`suff' = total(sqdiff_`suff')
					
					** c. Dividing by sample size and taking Sqrt
						gen rmse_`suff' = sqrt((sum_`suff'/1000))
						
					** d. Recording
						sum rmse_`suff'
						mat five_t5[`=cond("`suff'"=="hat",6,7)',1] = round(`r(mean)',.0001)
					}
				noi sum rmse_*
			
			** vii. RMSE Ratio
				mat five_t5[8,1] = round(five_t5[7,1]/five_t5[6,1],.0001)
				mat l five_t5
			}
		*****
		}
	*****

	** 4. Repeating For T=10 & T=20
		{
		** A. Starting T Loop
			qui foreach T in 10 20 {
				noi di "T=`T'"
				** i. Matrix to hold results
					mat results_t`T' = J(1000,3,.)
					mat coln results_t`T' = bhat sigma_hat sigma_tilde
				
				** ii. Creating Dataset
					{
					** a. Establishing i Structure
						clear
						set obs 500
						gen id = _n
					
					** b. Establishing t Structure
						expand `T'
						sort id
						egen t = seq(), f(1) t(`T')

					** c. Affirming Success
						isid id t
					}
				***** 
				
				** iii. Parts A-C Simulations
					{
					** a. Establishing sim loop
					qui forval sim = 1/1000 {
						if mod(`sim',100)==0{
							noi di "	Iteration `sim'..."
							}
						preserve
						** i. Creating RHS Variables
							gen x = rnormal(0,1)
							gen u = .
							d,s
							loc obs = `r(N)'
							qui forval i = 1/`obs' {
								loc sd = (abs(x[`i']))
								replace u = rnormal(0,`sd') in `i'
								}
						
						** ii. Estimating Mean Terms
							egen x_bar = mean(x), by(id)
							egen u_bar = mean(u), by(id)
						
						** iii. Creating Bhat Denominator Term
							{
							** a. Creating Squared difference term 
								gen dm_sqdiff_x = (x-x_bar)^2
								
							** b. Executing T & N Summations
								egen denom_sum = total(dm_sqdiff_x)
							}
						*****
						
						** iv. Creating Bhat Numerator Term
							{
							** a. Creating Differences Product term
								gen dm_product_xu = (x-x_bar)*(u-u_bar)
							
							** b. Executing T & N Summations
								egen numer_sum = total(dm_product_xu)
							}
						*****
						
						** v. Estimating Bhat using beta=1
							gen bhat = 1 + numer_sum / denom_sum
						
						** vi. Creating LHS Y given beta=1
							gen y = 1*x + u
							
						** vii. Creating Residual term now & mean across t
							gen u_hat = y - bhat*x
							egen u_hat_bar = mean(u_hat), by(id)
						
						** viii. Creating S_xx var 
							** (yes this is memory intensive but pedagogically nice)
								gen s_xx = denom_sum
						
						** ix. Creating Sigma hat
							{
							** a. Creating Inner product
								gen sighat_inner = (x-x_bar)*(u_hat-u_hat_bar)
							
							** b. Executing T summation & square 500 numbers made from 5 each
								egen sighat_tsum = total(sighat_inner), by(id)
								gen sighat_tsum_sq = sighat_tsum^2 if t==1
								** we need to be careful not to double count
							
							** c. Executing N summation, of just the 500 numbers above
								egen sighat_nsum = total(sighat_tsum_sq)
							
							** d. Multiply with 1/(s_xx^2) for variance
								gen var_hat = (1/s_xx^2)*sighat_nsum
							
							** e. Obtaining SD
								gen sd_hat = sqrt(var_hat)
							}
						*****
						
						** x. Creating Sigma Tilde
							{
							** a. Creating product of squared differences terms
								gen sqdiff_prod_xu_hat = ((x-x_bar)^2)*((u_hat-u_hat_bar)^2)
							
							** b. Executing T & N Summations
								egen sigtilde_sum = total(sqdiff_prod_xu_hat)
								
							** c. Multiply with 1/(s_xx^2) for variance
								gen var_tilde = (1/s_xx^2)*sigtilde_sum
							
							** d. Obtaining SD
								gen sd_tilde = sqrt(var_tilde)
							}
						*****
						
						** xi. Robustness checks
							{
							** a. We're going to compare to Stata's Estimates
								reg y x, r noconstant
								mat raw = r(table)
							
							** b. Comparing Stata's Beta to ours
								sum bhat
								loc beta = `r(mean)'
								if abs(`beta' - raw[1,1]) >= .05 {
									noi di "	- Robustnes Flag, Beta difference  `=abs(`beta'-raw[1,1])'"
									}
							
							** c. Sigma Hat
								sum sd_hat
								loc sd_hat = `r(mean)'
								if abs(`sd_hat' - raw[2,1]) >= .05 {
									noi di "	- Robustnes Flag, Sigma Hat difference  `=abs(`sd_hat'-raw[2,1])'"
									}
							
							** c. Sigma Tidle
							** No equivalent Stata command found
							/*	ren y x
								mat raw = r(table)
								sum sd_tilde
								loc sd_tilde = `r(mean)'
								if abs(`sd_tilde'-raw[1,2])>=.01 {
									noi di "	- Robustnes Flag, Sigma Hat difference  `=abs(`sd_tilde'-raw[1,2])'"
									}
							*/
							}
						*****
						
						** xii. Recording Results
							forval c = 1/3 {
								sum `=word("bhat sd_hat sd_tilde",`c')'
								mat results_t`T'[`sim',`c'] = round(`r(mean)',.00001)
								}
						
						** xiii. Preparing for next iteration
							restore
						} // End Sim Loop
		
					** b. Establishing Matrix for this set of simulations
						{
						** i. Empty matrix
							mat five_t`T' = J(8,1,.)
							
						** b. Adjusting row and column headers
							mat rown five_t`T' =	std_bhat	///
												bias_hat	bias_tilde	///
												std_hat		std_tilde	///
												rmse_hat	rmse_tilde	///
												rmse_ratio
							mat coln five_t`T' = T_`T'
						}
					*****
					
					** c. Working with Simulation Results
						{
						** i. Moving Results to Memory
							clear 
							svmat results_t`T', n(col)
						
						** ii. 5b, Bhat SD
							qui sum bhat
							noi di "Beta Hat has an SD of `=round(`r(sd)',.001)'"
						
						** iii. Recording in matrix, as scalar and as variable
							mat five_t`T'[1,1] = round(`r(sd)',.0001)
							sca def bhat_sd = `r(sd)'
							egen sd_bhat = sd(bhat)
						
						** iv. Bias, E(estimator)-true
							qui foreach suff in hat tilde {
								sum sigma_`suff'
								loc bias = `r(mean)'-`=scalar(bhat_sd)'
								noi di "Sigma `suff' has bias of `=round(`bias',.001)'"
								mat five_t`T'[`=cond("`suff'"=="hat",2,3)',1] = round(`bias',.0001)
					
								}
						
						** v. Std
							qui foreach suff in hat tilde {
								sum sigma_`suff'
								noi di "Sigma `suff' has an SD of `r(sd)'"
								mat five_t`T'[`=cond("`suff'"=="hat",4,5)',1] = round(`r(sd)',.0001)
								}
						
						** vi. RMSE
							foreach suff in hat tilde {
								** a. Creating Sq difference
									gen sqdiff_`suff' = (sigma_`suff'-sd_bhat)^2
								
								** b. Creating Sum
									egen sum_`suff' = total(sqdiff_`suff')
								
								** c. Dividing by sample size and taking Sqrt
									gen rmse_`suff' = sqrt((sum_`suff'/1000))
								
								** d. Recording
									sum rmse_`suff'
									mat five_t`T'[`=cond("`suff'"=="hat",6,7)',1] = round(`r(mean)',.0001)							
								}
							noi sum rmse_*
						
						** vii. RMSE Ratio
							mat five_t`T'[8,1] = round(five_t`T'[7,1]/five_t`T'[6,1],.0001)
							mat l five_t`T'
						}
					*****
					}
				*****
				} // End T Loop
		}
	*****

	** 5. Aggregating all the results into a Latex friendly table
		{
		** A. Combining Matrices from each set
			mat full = five_t5, five_t10, five_t20
		
		** B. Adjusting row and column headers
			mat rown full = sd_bhat	///
							bias_hat	bias_tilde	///
							std_hat		std_tilde	///
							rmse_hat	rmse_tilde	///
							rmse_ratio
			mat coln full = T_5 T_10 T_20
			
		** C. Displaying for Fun
			mat l full
			
		** D. Exporting to Latex
			esttab matrix(full, fmt(%15.4fc)) using "PS1_5_Results.tex", replace
		}
	*****
	}
*****

** Problem 6
	{
	 ** 1. Loading Data
		u "handguns.dta", clear
		isid year state
		
	** 2. QI Logged Dependent Variables
	qui	{
		** i. Creating Log Variables
			foreach var of varlist vio mur rob {
				gen ln_`var' = ln(`var')
				lab var ln_`var' "ln(`var')"
				}
				
		** ii. Executing Regressions, storing results
			eststo clear
			foreach var of varlist vio mur rob {
				eststo: reg ln_`var' i.shall, r
				}
		
		** iii. Exporting To Latex
			esttab using "PS1_6I_Results.tex",		///
				b(%15.3fc) se(%15.3fc) r2(%15.2fc)	///
				star compress label   				///
				noobs nogaps nonum nobaselevels 	///
				order(_cons 1.shall)				///
				coefl(	_cons 	"\hat{\beta}_0"		///
						1.shall	"\hat{\beta}_1(shall)") 	///
				replace
		}
	*****

	** 3. QII Including Control Variables
	qui	{
		** i. Executing Regressions, storing results
			eststo clear
			foreach var of varlist vio mur rob {
				eststo: reg ln_`var' i.shall incarc_rate density pop pm1029 avginc, r
				}
		
		** ii. Exporting To Latex
			esttab using "PS1_6II_Results.tex",		///
				b(%15.3fc) se(%15.3fc) r2(%15.2fc)	///
				star compress label   				///
				noobs nogaps nonum nobaselevels 	///
				order(_cons 1.shall)				///
				keep(_cons 1.shall)					///
				coefl(	_cons 	"\hat{\beta}_0"		///
						1.shall	"\hat{\beta}_1(shall)") 	///
				replace
		}
	*****
			
	** 4. QIV Taking Full Advantage of the Panel Setting
	qui	{
		** i. Creating New numeric State ID
			encode state, gen(num_stateid)
			
		** ii. Initializing Loop through Dependent Variables
			foreach var of varlist vio mur rob {
				noi di "`var'"
				** a. Setting up Results matrix for each variable now
					mat results_6IV_`var' 		= J(8,4,.)
					mat coln results_6IV_`var' 	= "Dependent=`var'&1" "2" "3" "Cluster Option"
					mat rown results_6IV_`var' 	= "\hat{\beta}_1(shall)" "(SE)" "State FE" "Year FE"	///
													"F-Test, State FE" "(PV)" "F-Test, Year FE" "(PV)"
				
				** i. Column 1 values
					reg ln_`var' i.shall incarc_rate density pop pm1029 avginc, r
					mat raw = r(table)
					mat results_6IV_`var'[1,1] = `=round(raw[1,"1.shall"],.001)'
					mat results_6IV_`var'[2,1] = `=round(raw[2,"1.shall"],.001)'
					mat results_6IV_`var'[3,1] = 0
					mat results_6IV_`var'[4,1] = 0
					
				** ii. Column 2
					reg ln_`var' i.shall						/// Treatment
							incarc_rate density pop pm1029 avginc	/// Controls
							i.num_stateid, r						// State FE
					mat raw = r(table)
					mat results_6IV_`var'[1,2] = `=round(raw[1,"1.shall"],.001)'
					mat results_6IV_`var'[2,2] = `=round(raw[2,"1.shall"],.001)'
					mat results_6IV_`var'[3,2] = 1
					mat results_6IV_`var'[4,2] = 0
					testparm i.num_stateid
					mat results_6IV_`var'[5,2] = round(r(F),.001)
					mat results_6IV_`var'[6,2] = round(r(p),.00001)
				
				** iii. Column 3
					reg ln_`var' i.shall							/// Treatment
							incarc_rate density pop pm1029 avginc	/// Controls
							i.num_stateid i.year, r						// State & Year FE
					mat raw = r(table)
					mat results_6IV_`var'[1,3] = `=round(raw[1,"1.shall"],.001)'
					mat results_6IV_`var'[2,3] = `=round(raw[2,"1.shall"],.001)'
					mat results_6IV_`var'[3,3] = 1
					mat results_6IV_`var'[4,3] = 1
					testparm i.num_stateid
					mat results_6IV_`var'[5,3] = round(r(F),.001)
					mat results_6IV_`var'[6,3] = round(r(p),.00001)
					testparm i.year
					mat results_6IV_`var'[7,3] = round(r(F),.001)
					mat results_6IV_`var'[8,3] = round(r(p),.00001)
					
				** iv. Column Cluster
					qui reg ln_`var' i.shall						/// Treatment
							incarc_rate density pop pm1029 avginc	/// Controls
							i.num_stateid i.year, cluster(num_stateid) r	// State & Year FE
					mat raw = r(table)
					mat results_6IV_`var'[1,4] = `=round(raw[1,"1.shall"],.001)'
					mat results_6IV_`var'[2,4] = `=round(raw[2,"1.shall"],.001)'
					mat results_6IV_`var'[3,4] = 1
					mat results_6IV_`var'[4,4] = 1
					
				** v. Displaying Results
					noi mat l results_6IV_`var'
				
				** vi. Exporting Matrix in Latex Friendly Format
					esttab matrix(results_6IV_`var') using "PS1_6IV_`var'_Results.tex", replace
				} // End Dep variable Loop
		}
	*****   
	}
*****

	*******
	* END *
	*******