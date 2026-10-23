*! version 1.1.0  13mar2019 Marginalized Zero-inflated Generalized Poisson LL functions
program define mzigp_llf
        version 6.0
	args todo b lnf g H sc1 sc2 sc3
quietly {
	tempname tau	
	tempvar eta zg 
	mleval `eta' = `b', eq(1)
	mleval `zg'  = `b', eq(2)
	mleval `tau' = `b', eq(3) scalar
	
	if `tau' < -18 { scalar `tau' = -18 }
	if `tau' >  18 { scalar `tau' =  18 }
	
	global SGLM_m = 1
	global MLINK = 2
		
	// We are using the GLM programs to handle the links, inverse links, and their derivatives. 
		// LOGIT :  glim_l02
		// PROBIT:  glim_l08
		// LOGLOG:  glim_l06
		// CLOGLOG: glim_l07
		

	/* Calculate the log-likelihood */

	tempname del gamma
	scalar `del' = (exp(2*`tau')-1) / (exp(2*`tau')+1)
	glim_l0$MLINK 1 `zg' `gamma'
	local y "$ML_y1"


	tempname maxdel mindel 
	tempvar  chk
	if 0 {
		qui gen double `chk' = (exp(`eta')-`y') 
		summ `chk' if $ML_samp, meanonly
		
		scalar `mindel' = max(-1, r(min))
	
		if `mindel' > -.001 {
			scalar `mindel' = -0.001
		}
	}
	else {
		qui gen double `chk' = exp(`eta')
		scalar `mindel' = max(-1, -r(max)/4)
	}
	scalar `mindel' = `mindel' + 0.001
	scalar `maxdel' =  0.999

	if `tau' > 300 | `del' > `maxdel' {
		scalar `del' = `maxdel'
	}
	if `del' < `mindel' {
		scalar `del' = `mindel'
	}

	tempname onemd
	scalar `onemd' = 1.0 - `del'

	tempvar mu den F G
	qui gen double `mu' = exp(`eta' - ln(1-`gamma'))	if $ML_samp
	qui gen double `den'= `onemd'*`mu'+`del'*`y' 		if $ML_samp
	qui gen double `F' = `onemd'*`mu'					if $ML_samp
	qui gen double `G' = `gamma'+(1-`gamma')*`F'		if $ML_samp
	
	#delimit ;
	mlsum `lnf' = cond(`y' == 0 ,
		ln(exp((-1+`del')*`mu')*(1-`gamma')+`gamma')
		,
		-(`y'*`del') - `F' - lngamma(`y'+1) + ln(1-`gamma') + ln(`F') + (-1+`y')*ln((`y'*`del')+`F')
		) ;
	#delimit cr 
	
	if `todo' == 0 { 
			exit 
		}
	
		forvalues k=1/3 {
			if "`sc`k''" == "" {
				tempvar sc`k'
				gen double `sc`k'' = .
			}
		}
		
		tempvar dmu1 v1 v2
		glim_l0$MLINK 2 `zg' `gamma' `dmu1'
		qui gen double `v1' = exp((2*(1-`gamma')*`mu')/((1+exp(2*`tau'))*(-1+`gamma')))	if $ML_samp
		qui gen double `v2' = ((1-`gamma')*`mu')/((1+exp(2*`tau'))*(-1+`gamma'))		if $ML_samp
		
	
		#delimit ;
		replace `sc1' = cond(`y' == 0,
				(2*`v1'*(1-`gamma')*`mu') / ((1+exp(2*`tau'))*(-`v1'+(-1+`v1')*`gamma')),
				(-4*(exp(2*(ln(1-`gamma')+ln(`mu'))))-`y'+exp(4*`tau')*`y'+(-1+exp(4*`tau'))*`y'*(`gamma'^2)+
				4*`y'*(1-`gamma')*`mu'-2*`y'*`gamma'*(-1+exp(4*`tau')+2*(1-`gamma')*`mu')) / 
				((1+exp(2*`tau'))*(-1+`gamma')*(`y'-exp(2*`tau')*`y'+(-1+exp(2*`tau'))*`y'*`gamma'-
				2*(1-`gamma')*`mu'))
			) if $ML_samp ;
			
		replace `sc2' = cond(`y' == 0,
				((1-`v1'+2*(exp((2*(1-`gamma')*`mu')/((1+exp(2*`tau'))*(-1+`gamma'))))*`v2')*`dmu1') / 
				(-`v1'*(-1+`gamma')+`gamma'),
				-((2*(1-`gamma')*`mu'*(1+(exp(2*`tau'))-2*`y'-(1+(exp(2*`tau'))-2*`y')*`gamma'+
				2*(1-`gamma')*`mu')*`dmu1') / ((1+exp(2*`tau'))*((-1+`gamma')^2)*(-`y'+
				(exp(2*`tau'))*`y'-(-1+exp(2*`tau'))*`y'*`gamma'+2*(1-`gamma')*`mu')))
			) if $ML_samp ;
			
		replace `sc3' = cond(`y' == 0,
				-((4*(exp(2*(`v2'+`tau')))*(1-`gamma')*`mu') / (((1+exp(2*`tau'))^2)*(-`v1'+(-1+`v1')*`gamma'))),
				-((2*(exp(2*`tau'))*(-4*exp(2*(ln(1-`gamma')+ln(`mu')))+`y'+2*exp(2*`tau')*`y'+
				exp(4*`tau')*`y'-4*(`y'^2)+(1+2*exp(2*`tau')+exp(4*`tau')-4*`y')*`y'*(`gamma'^2)+
				8*`y'*(1-`gamma')*`mu'+2*`y'*`gamma'*(-1-2*exp(2*`tau')-exp(4*`tau')+4*`y'-4*(1-`gamma')*`mu'))) / 
				(((1+exp(2*`tau'))^2)*(-1+`gamma')*(`y'-exp(2*`tau')*`y'+(-1+exp(2*`tau'))*`y'*`gamma'-2*(1-`gamma')*`mu')))
			) if $ML_samp ;
		#delimit cr

		tempname d1 d2 d3
		mlvecsum `lnf' `d1' = `sc1', eq(1)
		mlvecsum `lnf' `d2' = `sc2', eq(2)
		mlvecsum `lnf' `d3' = `sc3', eq(3)
		
		matrix `g' = (`d1', `d2', `d3')
		
			
		if (`todo'==1 | `lnf'>=.) {
			exit
		}
		
		tempvar dmu2
		glim_l0$MLINK 3 `zg' `gamma' `dmu2'

		tempvar t1 t2 t3 sc12 sc13 sc23
			
		#delimit ;
		gen double `t1' = cond(`y'==0, 
				(2*`v1'*(1-`gamma')*`mu'*((1+exp(2*`tau'))*`v1'+(1+exp(2*`tau'))*(-1+`v1')*(`gamma'^2)+
				`gamma'*((1+exp(2*`tau')-2*exp(2*(`v2'+`tau'))-2*`v1'-2*(1-`gamma')*`mu')))) / 
				(((1+exp(2*`tau'))^2)*(-1+`gamma')*((`v1'+`gamma'-`v1'*`gamma')^2))
				,
				(2*`v2'*(4*exp(2*(ln(1-`gamma')+ln(`mu')))+exp(4*`tau')*`y'-2*exp(2*`tau')*(`y'^2)
				+`y'*(-1+2*`y')+(-1+exp(2*`tau'))*(1+exp(2*`tau')-2*`y')*`y'*(`gamma'^2)-
				4*`y'*(1-`gamma')*`mu'+4*exp(2*`tau')*`y'*(1-`gamma')*`mu'-
				2*(-1+exp(2*`tau'))*`y'*`gamma'*(1+exp(2*`tau')-2*`y'+2*(1-`gamma')*`mu')))	/ 
				(`y'-exp(2*`tau')*`y'+(-1+exp(2*`tau'))*`y'*`gamma'-2*(1-`gamma')*`mu')^2
			) if $ML_samp ;
			
		gen double `sc12' = cond(`y'==0, 
				 (2*`v1'*(1-`gamma')*`mu'*(-(1+exp(2*`tau'))*(-1+`v1')-
				 (1+exp(2*`tau'))*(-1+`v1')*(`gamma'^2)+2*`gamma'*(-1-exp(2*`tau')+exp(2*(`v2'+`tau'))+
				 `v1'+(1-`gamma')*`mu'))*`dmu1') / 
				 (((1+exp(2*`tau'))^2)*((-1+`gamma')^2)*((`v1'+`gamma'-`v1'*`gamma')^2))
				,
				 -((2*(1-`gamma')*`mu'*(4*exp(2*(ln(1-`gamma')+ln(`mu')))+exp(4*`tau')*`y'-
				 2*exp(2*`tau')*(`y'^2)+`y'*(-1+2*`y')+(-1+exp(2*`tau'))*(1+exp(2*`tau')-2*`y')*`y'*(`gamma'^2)-
				 4*`y'*(1-`gamma')*`mu'+4*exp(2*`tau')*`y'*(1-`gamma')*`mu'-2*(-1+exp(2*`tau'))*`y'*`gamma'*(1+
				 exp(2*`tau')-2*`y'+2*(1-`gamma')*`mu'))*`dmu1') / 
				 ((1+exp(2*`tau'))*((-1+`gamma')^2)*((`y'-exp(2*`tau')*`y'+(-1+exp(2*`tau'))*`y'*`gamma'-
				 2*(1-`gamma')*`mu')^2)))
			) if $ML_samp ;
			
		gen double `sc13' = cond(`y'==0, 
				(4*exp(2*(`v2'+`tau'))*(1-`gamma')*`mu'*(-(1+exp(2*`tau'))*`v1'-
				(1+exp(2*`tau'))*(-1+`v1')*(`gamma'^2)+`gamma'*(-1-exp(2*`tau')+
				2*exp(2*(`v2'+`tau'))+2*`v1'+2*(1-`gamma')*`mu'))) /
				(((1+exp(2*`tau'))^3)*(-1+`gamma')*((`v1'+`gamma'-`v1'*`gamma')^2))
				,
				-((4*exp(2*`tau')*(1-`gamma')*`mu'*(4*exp(2*(ln(1-`gamma')+ln(`mu')))+
				`y'+exp(4*`tau')*`y'-2*exp(2*`tau')*`y'*(-1+2*`y')+(1+exp(4*`tau')+exp(2*`tau')*(2-4*`y'))*`y'*(`gamma'^2)-
				4*`y'*(1-`gamma')*`mu'+4*exp(2*`tau')*`y'*(1-`gamma')*`mu'+
				2*`y'*`gamma'*(-1-exp(4*`tau')+exp(2*`tau')*(-2+4*`y')+2*(1-`gamma')*`mu'-
				2*exp(2*`tau')*(1-`gamma')*`mu'))) / 
				(((1+exp(2*`tau'))^2)*(-1+`gamma')*((`y'-exp(2*`tau')*`y'+(-1+exp(2*`tau'))*`y'*`gamma'-
				2*(1-`gamma')*`mu')^2)))
			) if $ML_samp ;
			
		gen double `t2' = cond(`y'==0, 
				(-((-1+`v1'-2*exp((2*(1-`gamma')*`mu')/((1+exp(2*`tau'))*(-1+`gamma')))*`v2')^2)*(`dmu1'^2)+
				(-`v1'*(-1+`gamma')+`gamma')*(-((4*(exp(2*(`v2'+ln(1-`gamma')+ln(`mu'))))*(`dmu1'^2))/
				(((1+exp(2*`tau'))^2)*((-1+`gamma')^3)))+(1-`v1'+2*exp((2*(1-`gamma')*`mu')/
				((1+exp(2*`tau'))*(-1+`gamma')))*`v2')*`dmu2')) /
				((-`v1'*(-1+`gamma')+`gamma')^2)
				,
				(2*(1-`gamma')*`mu'*(2*(`dmu1'^2)-((2*(1+exp(2*`tau'))*(-1+`y')*(1-`gamma')*(-1+`gamma')*`mu'*(`dmu1'^2))/
				((`y'-exp(2*`tau')*`y'+(-1+exp(2*`tau'))*`y'*`gamma'-2*(1-`gamma')*`mu')^2))-((2*(1+exp(2*`tau'))*(-1+`y')*(-1+
				`gamma')*(`dmu1'^2))/(`y'-exp(2*`tau')*`y'+(-1+exp(2*`tau'))*`y'*`gamma'-2*(1-`gamma')*`mu'))-
				(-1+`gamma')*`dmu2'+ (((1+exp(2*`tau'))*(-1+`y')*((-1+`gamma')^2)*`dmu2')/(`y'-exp(2*`tau')*`y'+
				(-1+exp(2*`tau'))*`y'*`gamma'-2*(1-`gamma')*`mu')))) / ((1+exp(2*`tau'))*((-1+`gamma')^3))
			) if $ML_samp ;
			
		gen double `sc23' = cond(`y'==0, 
				(4*exp(2*(`v2'+`tau'))*(1-`gamma')*`mu'*((1+exp(2*`tau'))*(-1+`v1')+(1+exp(2*`tau'))*(-1+`v1')*(`gamma'^2)-
				2*`gamma'*(-1-exp(2*`tau')+exp(2*(`v2'+`tau'))+`v1'+(1-`gamma')*`mu'))*`dmu1') / 
				(((1+exp(2*`tau'))^3)*((-1+`gamma')^2)*((`v1'+`gamma'-`v1'*`gamma')^2))
				,
				(4*exp(2*`tau')*(1-`gamma')*`mu'*(4*exp(2*(ln(1-`gamma')+ln(`mu')))+`y'+exp(4*`tau')*`y'-
				2*exp(2*`tau')*`y'*(-1+2*`y')+(1+exp(4*`tau')+exp(2*`tau')*(2-4*`y'))*`y'*(`gamma'^2)-4*`y'*(1-`gamma')*`mu'+
				4*exp(2*`tau')*`y'*(1-`gamma')*`mu'+2*`y'*`gamma'*(-1-exp(4*`tau')+exp(2*`tau')*(-2+4*`y')+
				2*(1-`gamma')*`mu'-2*exp(2*`tau')*(1-`gamma')*`mu'))*`dmu1') / 
				(((1+exp(2*`tau'))^2)*((-1+`gamma')^2)*((`y'-exp(2*`tau')*`y'+(-1+exp(2*`tau'))*`y'*`gamma'-2*(1-`gamma')*`mu')^2))
			) if $ML_samp ;
			
		gen double `t3' = cond(`y'==0, 
				-((8*exp(2*(`v2'+`tau'))*(1-`gamma')*`mu'*(-(-1+exp(4*`tau'))*`v1'-(-1+exp(4*`tau'))*(-1+`v1')*(`gamma'^2)+
				`gamma'*(1-exp(4*`tau')+2*exp(2*`v2'+4*`tau')-2*`v1'+2*exp(2*`tau')*(1-`gamma')*`mu'))) / 
				(((1+exp(2*`tau'))^4)*(-1+`gamma')*(`v1'+`gamma'-`v1'*`gamma')^2))
				,
				(4*exp(2*`tau')*(-8*exp(3*(ln(1-`gamma')+ln(`mu')))+8*exp(2*`tau'+3*(ln(1-`gamma')+ln(`mu')))+
				20*exp(2*(ln(1-`gamma')+ln(`mu')))*`y'-20*exp(2*(`tau'+ln(1-`gamma')+ln(`mu')))*`y'+
				8*exp(4*`tau'+2*(ln(1-`gamma')+ln(`mu')))*`y'-exp(6*`tau')*(`y'^2)+
				(`y'^2)*(-1+4*`y')-exp(2*`tau')*(`y'^2)*(3+4*`y')+exp(4*`tau')*(`y'^2)*(-3+8*`y')-(`y'^2)*(-1-exp(6*`tau')+
				4*`y'-exp(2*`tau')*(3+4*`y')+exp(4*`tau')*(-3+8*`y'))*(`gamma'^3)+
				2*exp(6*`tau')*`y'*(1-`gamma')*`mu'-2*exp(4*`tau')*`y'*(-3+8*`y')*(1-`gamma')*`mu'-
				2*`y'*(-1+8*`y')*(1-`gamma')*`mu'+2*exp(2*`tau')*`y'*(3+8*`y')*(1-`gamma')*`mu'+
				`y'*(`gamma'^2)*(-3*exp(6*`tau')*`y'+3*`y'*(-1+4*`y')-3*exp(2*`tau')*`y'*(3+4*`y')+
				3*exp(4*`tau')*`y'*(-3+8*`y')+2*exp(6*`tau')*(1-`gamma')*`mu'+(2-16*`y')*(1-`gamma')*`mu'+
				exp(4*`tau')*(6-16*`y')*(1-`gamma')*`mu'+2*exp(2*`tau')*(3+8*`y')*(1-`gamma')*`mu')-
				`y'*`gamma'*(20*exp(2*(ln(1-`gamma')+ln(`mu')))-20*exp(2*(`tau'+ln(1-`gamma')+ln(`mu')))+
				8*exp(4*`tau'+2*(ln(1-`gamma')+ln(`mu')))-3*exp(6*`tau')*`y'+3*`y'*(-1+4*`y')-
				3*exp(2*`tau')*`y'*(3+4*`y')+3*exp(4*`tau')*`y'*(-3+8*`y')+4*exp(6*`tau')*(1-`gamma')*`mu'+
				(4-32*`y')*(1-`gamma')*`mu'-4*exp(4*`tau')*(-3+8*`y')*(1-`gamma')*`mu'+
				4*exp(2*`tau')*(3+8*`y')*(1-`gamma')*`mu'))) / 
				(((1+exp(2*`tau'))^3)*(-1+`gamma')*(`y'-exp(2*`tau')*`y'+(-1+exp(2*`tau'))*`y'*`gamma'-2*(1-`gamma')*`mu')^2)
			) if $ML_samp ;
	
		#delimit cr
		
		tempname d11 d12 d13 d22 d23 d33
		mlmatsum `lnf' `d11' = `t1'  , eq(1)
		mlmatsum `lnf' `d12' = `sc12', eq(1,2)
		mlmatsum `lnf' `d13' = `sc13', eq(1,3)
		mlmatsum `lnf' `d22' = `t2'  , eq(2)
		mlmatsum `lnf' `d23' = `sc23', eq(2,3)
		mlmatsum `lnf' `d33' = `t3'  , eq(3)
		
		matrix `H' = (`d11', `d12', `d13' \ `d12'', `d22', `d23' \ `d13'', `d23'', `d33')
	}
end
