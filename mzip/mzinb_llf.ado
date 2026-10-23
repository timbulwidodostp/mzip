*! version 1.1.0  13mar2019 Marginalized Zero-inflated Negative Binomial LL functions
program define mzinb_llf
	args todo b lnf gvec H sc1 sc2 sc3
        
	quietly {
		local y "$ML_y1" 

		tempvar xb zg lnalpha gamma 
		mleval `xb'      = `b', eq(1)
		mleval `zg'      = `b', eq(2)
		mleval `lnalpha' = `b', eq(3) scalar
		
		glim_l0$MLINK 1 `zg' `gamma' 

		tempname alpha
		scalar `alpha' = max(min(exp(`lnalpha'),18),-18)
		scalar `lnalpha' = ln(`alpha')

		tempvar p mu m 
		gen double `m'  = 1/`alpha'          		if $ML_samp
		gen double `mu' = exp(`xb' - ln(1-`gamma'))	if $ML_samp
		gen double `p'  = 1/(1+`alpha'*`mu') 		if $ML_samp
		
		// We are using the GLM programs to handle the links, inverse links, and their derivatives. 
		// LOGIT :  glim_l02
		// PROBIT:  glim_l08
		// LOGLOG:  glim_l06
		// CLOGLOG: glim_l07
		
		mlsum `lnf' = cond(`y' == 0, 				///
				ln((`p'^`m')*(1-`gamma')+`gamma')	///
				, 					///
				- lngamma(`m') - lngamma(`y'+1) +       ///
				lngamma(`y'+`m') + `y'*ln(1-`p') +	///
				`m'*ln(`p') + ln(1-`gamma')		///
				) 
		
		if `todo' == 0 { 
			exit 
		}
	
		forvalues k=1/3 {
			if "`sc`k''" == "" {
				tempvar sc`k'
				gen double `sc`k'' = .
			}
		}
	
		tempvar dmu1 F G
		glim_l0$MLINK 2 `zg' `gamma' `dmu1'
		gen double `F' = (`p'^`m')
		gen double `G' = -`F'*(-1+`gamma')+`gamma'
		
		#delimit ;
		replace `sc1' = cond(`y' == 0,
				-(((`p'^(1+`m'))*(1-`gamma')*`mu') / (`G')),
				-((-`y'+`y'*`gamma'+(1-`gamma')*`mu') / (1-`gamma'+`alpha'*(1-`gamma')*`mu'))
			) if $ML_samp ;
			
		replace `sc2' = cond(`y' == 0,
				((1-`F'-(`p'^(`m'+1))*`mu')*`dmu1') / (`G')
				,
				((1-`y'+(-1+`y')*`gamma'+(1-`gamma')*`mu'+`alpha'*(1-`gamma')*`mu')*`dmu1') / 
				((-1+`gamma')*(1-`gamma'+`alpha'*(1-`gamma')*`mu'))
			) if $ML_samp ;
			
		replace `sc3' = cond(`y' == 0,
				(`F'*(1-`gamma')*(-(((1-`gamma')*`mu') / (1-`gamma'+`alpha'*(1-`gamma')*`mu')) - 
				`m'*ln(`p'))) / (`G')
				,
				((`y'-`y'*`gamma'-(1-`gamma')*`mu') / (1-`gamma'+`alpha'*(1-`gamma')*`mu')) + 
				`m'*digamma(`m')-`m'*digamma(`m'+`y')-`m'*ln(`p')
			) if $ML_samp ;
		#delimit cr

		tempname d1 d2 d3
		mlvecsum `lnf' `d1' = `sc1', eq(1)
		mlvecsum `lnf' `d2' = `sc2', eq(2)
		mlvecsum `lnf' `d3' = `sc3', eq(3)
		
		matrix `gvec' = (`d1', `d2', `d3')

		if (`todo'==1 | `lnf'>=.) {
			exit
		}
		
		tempvar dmu2
		glim_l0$MLINK 3 `zg' `gamma' `dmu2'

		tempvar t1 t2 t3 sc12 sc13 sc23
		
		#delimit ;
		gen double `t1' = cond(`y'==0, 
				-((`F'*(1-`gamma')*(-1+`gamma')*`mu'*(-`F'-(-1+`F')*(`gamma'^2)+
				`gamma'*(-1+2*`F'+(1-`gamma')*`mu'))) / (((`F'+`gamma'-`F'*`gamma')^2)*((1-`gamma'+
				(`alpha')*(1-`gamma')*`mu')^2)))
				,
				((1+(`alpha')*`y')*(1-`gamma')*(-1+`gamma')*`mu') / 
				((1-`gamma'+(`alpha')*(1-`gamma')*`mu')^2)
			) if $ML_samp ;
			
		gen double `sc12' = cond(`y'==0, 
				(`F'*(1-`gamma')*`mu'*(1-`F'-(-1+`F')*(`gamma'^2)+(`alpha')*(1-`gamma')*`mu'+
				`gamma'*(-2+2*`F'+(1-`gamma')*`mu'))*`dmu1') / (((`F'+`gamma'-`F'*`gamma')^2)*((1-`gamma'+
				(`alpha')*(1-`gamma')*`mu')^2))
				,
				-((((1+`alpha'*`y'))*(1-`gamma')*`mu'*`dmu1') / 
				(1-`gamma'+(`alpha')*(1-`gamma')*`mu')^2)
			) if $ML_samp ;
			
		gen double `sc13' = cond(`y'==0, 
				(`m'*`F'*(1-`gamma')*(-1+`gamma')*`mu'*(-`alpha'^2*`F'*(1-`gamma')*`mu'+
				(`gamma'^2)*ln(`p')+`gamma'*((`alpha')*(-1+(`alpha')*(-1+`F'))*(1-`gamma')*`mu'- 
				(1+(`alpha')*(1-`gamma')*`mu')*ln(`p')))) / (((`F'+`gamma'-`F'*`gamma')^2)*((1-`gamma'+
				(`alpha')*(1-`gamma')*`mu')^2))
				,
				((`alpha')*(1-`gamma')*`mu'*(-`y'+`y'*`gamma'+(1-`gamma')*`mu')) / 
				((1-`gamma'+(`alpha')*(1-`gamma')*`mu')^2)
			) if $ML_samp ;
			
		gen double `t2' = cond(`y'==0, 
					(
						-(-1+`F'+`p'^(1+`m')*`mu')^2 * `dmu1'^2 + 
						(-`F'*(-1+`gamma')+`gamma') * (
							exp(`lnalpha'+2*(ln(1-`gamma')+ln(`mu')))*(1+`m')*`p'^(2+`m')*`dmu1'^2 / (1-`gamma')^3 + 
							`dmu2' - `F'*`dmu2' - `p'^(1+`m')*`mu'*`dmu2'
						)
					) / (-`F'*(-1+`gamma')+`gamma')^2,

					(
						-(
							1 + exp(2*(`lnalpha'+ln(1-`gamma')+ln(`mu'))) + exp(`lnalpha'+2*(ln(1-`gamma')+ln(`mu'))) - 
							`y' - (-1+`y')*`gamma'^2 + 2*(1-`gamma')*`mu' + 2*`alpha'*(1-`gamma')*`mu' - 
							2*`gamma'*(1-`y'+(1-`gamma')*`mu' + `alpha'*(1-`gamma')*`mu')
						)*`dmu1'^2 + 
						(-1+`gamma')*(1-`gamma'+`alpha'*(1-`gamma')*`mu')*(
							1-`y'+(-1+`y')*`gamma'+(1-`gamma')*`mu'+`alpha'*(1-`gamma')*`mu'
						) * `dmu2'
					) / ((-1+`gamma')^2*(1-`gamma'+`alpha'*(1-`gamma')*`mu')^2) 
			) if $ML_samp ;
			
		gen double `sc23' = cond(`y'==0, 
				-(
					(
						`m'*`F'*(
							-`alpha'*(1-`gamma')*`mu'*(1+`alpha'*(1+`F')*(1-`gamma')*`mu') + 
							`gamma'^2*(-1+(1-`gamma')*`mu')*ln(`p') - 
							(1+`alpha'*(1-`gamma')*`mu')^2*ln(`p') + `gamma'*(
								`alpha'*(1-`gamma')*`mu'*(
									1-(1-`gamma')*`mu'+`alpha'*(-1+`F')*(1-`gamma')*`mu'
								) - (-2+(1-`gamma')*`mu')*(1+`alpha'*(1-`gamma')*`mu')*ln(`p')
							)
						) * `dmu1'
					) / 
					((`F'+`gamma'-`F'*`gamma')^2*(1-`gamma'+`alpha'*(1-`gamma')*`mu')^2)
				), 
				-(
					`alpha'*(1-`gamma')*`mu'*(-`y'+`y'*`gamma'+(1-`gamma')*`mu')*`dmu1' 
					
				) / ((-1+`gamma')*(1-`gamma'+`alpha'*(1-`gamma')*`mu')^2)
			) if $ML_samp ;
				
		gen double `t3' = cond(`y'==0, 
				(1/(-`F'*(-1+`gamma')+`gamma')^2)*`F'*(1-`gamma')*(-`F'*(1-`gamma')*(((((1-`gamma')*`mu')/
				(1-`gamma'+`alpha'*(1-`gamma')*`mu'))+`m'*ln(`p'))^2)+
				(-`F'*(-1+`gamma')+`gamma')*(((((1-`gamma')*`mu')/(1-`gamma'+`alpha'*(1-`gamma')*`mu'))+`m'*ln(`p'))^2)+
				(-`F'*(-1+`gamma')+`gamma')*((((1-`gamma')*`mu'*(1-`gamma'+2*`alpha'*(1-`gamma')*`mu')) / 
				((1-`gamma'+`alpha'*(1-`gamma')*`mu')^2))+`m'*ln(`p')))
				,	
				((exp(`lnalpha'+2*(ln(1-`gamma')+ln(`mu'))))/((1-`gamma'+`alpha'*(1-`gamma')*`mu')^2))-
				((`y'*(-1+`gamma')^2)/((1-`gamma'+`alpha'*(1-`gamma')*`mu')^2))+
				((2*`alpha'*`y'*(1-`gamma')*(-1+`gamma')*`mu')/((1-`gamma'+`alpha'*(1-`gamma')*`mu')^2))+
				((`y'-`y'*`gamma')/(1-`gamma'+`alpha'*(1-`gamma')*`mu'))+
				(((1-`gamma')*`mu')/(1-`gamma'+`alpha'*(1-`gamma')*`mu'))-
				`m'*digamma(`m')+`m'*digamma(`m'+`y')+`m'*ln(`p')-(exp(-2*`lnalpha'))*trigamma(`m')+
				(exp(-2*`lnalpha'))*trigamma(`m'+`y')
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
