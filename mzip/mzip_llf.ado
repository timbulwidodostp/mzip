*! version 1.1.0 13mar2019  Marginalized Zero-inflated Poisson LL functions
program define mzip_llf
	args todo b lnf g H sc1 sc2
	quietly {
		tempvar xb zg 
		mleval `xb' = `b', eq(1)
		mleval `zg' = `b', eq(2)

		// We are using the GLM programs to handle the links, inverse links, and their derivatives. 
		// LOGIT :  glim_l02
		// LOG   :  glim_l03
		// PROBIT:  glim_l08
		// LOGLOG:  glim_l06
		// CLOGLOG: glim_l07

		global SGLM_m = 1
		
		tempvar  mu gamma 
		glim_l0$MLINK 1 `zg' `gamma' 
		gen double `mu' = exp(`xb' - ln(1-`gamma'))

		local y "$ML_y1"
		
			#delimit ;
			mlsum `lnf' = cond(`y'==0,													
				/* y=0  */	ln(-exp(-`mu')*(-1+`gamma')+`gamma'), 						
				/* y!=0 */	-`mu'-lngamma(1+`y')+ln(1-`gamma')+`y'*ln(`mu')	
			) ;
			#delimit cr
				
		if (`todo'==0) { 
			exit
		}

		forvalues k=1/2 {
			if "`sc`k''" == "" {
				tempvar sc`k'
				gen double `sc`k'' = .
			}
		}
	
		tempvar dmudeta2 
		glim_l0$MLINK 2 `zg' `gamma' `dmudeta2'
		
			#delimit ;
			replace `sc1' = cond(`y'==0,
				((-1+`gamma')*`mu') / (1+(-1+exp(`mu'))*`gamma') , 
				`y'-`mu'
			) ;
			replace `sc2' = cond(`y'==0,
				((-1+exp(`mu')-`mu')*`dmudeta2') / (1+(-1+exp(`mu'))*`gamma') ,
				((1-`y'+`mu')*`dmudeta2') / (-1+`gamma')
			) ;
			#delimit cr
		
		tempname d1 d2 d3
		$ML_ec mlvecsum `lnf' `d1' = `sc1', eq(1)
		$ML_ec mlvecsum `lnf' `d2' = `sc2', eq(2)
				
		$ML_ec matrix `g' = (`d1', `d2')
		
		if (`todo'==1 | `lnf'>=.) {
			exit
		}
		
		tempvar sc11 sc12 sc22 dmu2deta2
		glim_l0$MLINK 3 `zg' `gamma' `dmu2deta2'
		
			#delimit ;
			gen double `sc11' = cond(`y'==0, 
				-((-1+`gamma')*(-1+`gamma'*(1+exp(`mu')*(-1+`mu')))*`mu' / 
				(1+(-1+exp(`mu'))*`gamma')^2 ) , 
				-`mu' 
				) ;
				
			gen double `sc12' = cond(`y'==0, 
				(-1+`gamma'+(exp(`mu')*(1+`gamma'*(-1+`mu'))))*`mu'*`dmudeta2' / 
				(1+(-1+exp(`mu'))*`gamma')^2 , 
				(`mu'*`dmudeta2') / (-1+`gamma')
				) ;
				
			gen double `sc22'   = cond(`y'==0, 
				(exp(-2*`mu')*(-((1-exp(`mu')+`mu')^2)*(`dmudeta2'^2)+ 
				(((1+(-1+exp(`mu'))*`gamma')*((-`mu'^2)*(`dmudeta2'^2)+
				(-1+`gamma')*(-1+exp(`mu')-`mu')*`dmu2deta2')) / (-1+`gamma')))) / 
				((-exp(-`mu')*(-1+`gamma')+`gamma')^2) , 
				((-1+`y'-2*`mu')*(`dmudeta2'^2)-(-1+`gamma')*(-1+`y'-`mu')*`dmu2deta2') / 
				((-1+`gamma')^2)
				) ;
			
			#delimit cr
		
		tempname d11 d12 d22 
		mlmatsum `lnf' `d11' = `sc11', eq(1)
		mlmatsum `lnf' `d12' = `sc12', eq(1,2)
		mlmatsum `lnf' `d22' = `sc22', eq(2)
		
		matrix `H' = (`d11', `d12' \ `d12'', `d22')
	}
end

