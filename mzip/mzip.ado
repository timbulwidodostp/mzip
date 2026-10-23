*! version 1.1.0  13mar2019
program define mzip, eclass byable(onecall) prop(ml_score eform svyb svyj svyr)
	if _by() {
		local BY `"by `_byvars'`_byrc0':"'
	}
	`BY' _vce_parserun mzip, mark(INFlate OFFset CLuster) : `0'
	if "`s(exit)'" != "" {
		version 8: ereturn local cmdline `"mzip `0'"'
		exit
	}

	version 6.0, missing
	if replay() {
		if "`e(cmd)'" != "mzip" { error 301 }
		if _by() { 
			error 190 
		}
		Display `0'
		_prefix_footnote
		exit
	}
	if _caller() >= 11 {
		local vv : di "version " string(_caller()) ":"
	}
	`vv' `BY' Estimate `0'
	version 10: ereturn local cmdline `"mzip `0'"'
end

program define Estimate, eclass byable(recall)
	/* 	   
	   Model 0 is the target
	   We do not support constraints
	   
	   Regression coefficients and/or loglikelihoods from:
	   0) Zero-inf marginal Poisson model: FULL             B0_r, B0_0 with Off0_B, and G0_r, G0_0 with Off0_G
	   1) Zero-inf marginal Poisson model: CONSTANT-ONLY          B1_0 with Off0_B, and G1_r, G1_0 with Off0_G 
   
	   If Model 0 includes a constant 
	   (-nocons- option is not specified):
	       If either one of -nolr- or -from()- is not specifed: estimate Model 1 and get LL/pLL
		   
	       If there are no given starting values
		   (-from()- option is not specified):
               Model 1 provides candidate starting values for Model 0: {0, B1_0, G1_r, G1_0}
			   
		Select search(on)
	*/
	version 6.0, missing
	local cmd "`0'"
	local awopt = cond(_caller()<7,"aw","")
	#delimit ;
	syntax varlist(fv) [if] [in] [fweight `awopt' pweight iweight/] 
		, INFlate(string) [ 
				/* general options */
				noCONstant Level(cilevel) 
				
				/* optimization options */
				Robust CLuster(varname) SCore(string) EFORM 
				MLMETHOD(string) OFFset(varname numeric) noLRtest
				EXPosure(varname numeric) noLOg noDISPLAY
				FROM(string) CRITTYPE(passthru)
				ILink(string)
				
				/* additional optimization options */
				* 
		] ;
	#delimit cr
	
	local fvops1 = "`s(fvops)'" == "true" | _caller() >= 11
	if _by() { 
		_byoptnotallowed score() `"`score'"'
	}
	
	GetLink "`ilink'"
	if `r(rc)' {
		noi di as err "Illegal link specification in ilink(`ilink')"
		exit 198
	}
	local ilnk "`r(ilink)'"
	
	marksample touse
	markout `touse' `cluster', strok
	/* also see 2nd markout below for bionomial denominator */
	/* also see 3rd markout below for inflate varlist */

	tokenize `varlist'
	local dep "`1'"
	
	_fv_check_depvar `dep'
	mac shift
	local ind "`*'"

	if "`score'" != "" { 
		local n : word count `score'
		if `n'==1 & substr("`score'",-1,1)=="*" { 
			local score = substr("`score'",1,length("`score'")-1)
			local score `score'1 `score'2 
			local n 2
		}
		if `n' != 2 {
			di in red "score():  you must specify two new variables"
			exit 198
		}
	}

	if "`weight'" != "" {
		local wtexp `"[`weight'=`exp']"'
		local wtype "`weight'"
	}
	else    local exp 1

	_get_diopts diopts options, `options'
	local diopts `diopts' level(`level') `eform'
	mlopts mlopt, `options'
	local coll `s(collinear)'
	local cns `s(constraints)'
	local mlopt `mlopt' `crittype'
	
	if "`cluster'" != "" { local clopt "cluster(`cluster')" }
	_vce_parse, argopt(CLuster) opt(Robust oim opg) old: `wtexp', `vce' `clopt' `robust'
	local cluster `r(cluster)'
	local robust `r(robust)'

	if "`cluster'" != "" { local clopt "cluster(`cluster')" }

	local nc  "`constan'"

	PrseOff `inflate'
	local inflate "`s(inflate)'"
	local off2 "`s(off2)'"
	local nc2  "`s(nc2)'"
	if "`off2'" != "" {
		local off2a "offset(`off2')"
		local off2s "`off2'"
	}
	local fvops2 = "`s(fvops)'" == "true"
	
	if "`exposur'" != "" {
		tempvar off
		gen double `off' = ln(`exposur') if `touse'
		local offo "offset(`off')"
		local offstr "ln(`exposur')"
	}

	if "`offset'" != "" {
		if "`offstr'" != "" {
			di in red "cannot specify both exposure() and offset()"
			exit 198
		}
		tempvar off
		gen double `off' = `offset' if `touse'
		local offo "offset(`off')"
		local moff "-`off'"
		local poff "+`off'"
		local eoff "*exp(`off')"
		if "`offstr'" == "" { local offstr "`offset'" }
	}

	markout `touse' `inflate' `off2' `off' 
	
	if (`fvops1' | `fvops2') {
		local rmcoll "version 11: _rmcoll"
		local vv : di "version " string(max(11,_caller())) ":"
		local mm    e2
		local negh  ""
		local fvexp expand
	}
	else {
		local rmcoll _rmcoll
		local vv "version 8.1:"
		local mm d2
	}

	`rmcoll' `ind' if `touse' `wtexp' , `nc' `coll' `fvexp'
	local ind "`r(varlist)'"
	`rmcoll' `inflate' if `touse' `wtexp' , `nc2' `coll' `fvexp'
	local inflate "`r(varlist)'"

	// calculate correct df
	local df 0
	foreach var of local inflate {
		_ms_parse_parts `var'
		if !`r(omit)' {
			local ++df
		}
	}

	if "`nc2'" == ""     { local df = `df' + 1} 
	if "`mlmetho'" == "" { local mlmetho "`mm'" }
	if "`score'" != "" {
		local svar1 : word 1 of `score'
		local svar2 : word 2 of `score'
		confirm new variable `score'
		tempvar sc1 sc2 
		local scopt "score(`sc1' `sc2')"
	}

	qui count if `dep'<0 & `touse'
	if r(N) { 
		di in red "`dep'<0 in `r(N)' obs."
		exit 459
	}

	qui { 
		if "`wtype'" != "fweight" {
			count if `touse'
			local nobs = r(N)
			count if `touse' & `dep'==0
			local nzero = r(N)
		} 
		else {
			summ `dep' `wtexp' if `touse', meanonly
			local nobs = r(N)
			summ `dep' `wtexp' if `touse' & `dep'==0, meanonly
			local nzero = r(N)
		}
	}
	
	if `nzero'==`nobs' {
		di in red "`dep' never varies; always equal to zero"
		exit 459
	}
	if `nzero'==0 {
		di in bl _n "(note: " in ye "`dep'" in bl " never equal to zero; use glm instead)"
	}

	if "`log'`display'" == "" { local lll "noisily" }
	else                      { local lll "quietly" }
	if "`log'" == ""          { local log "noisily" }
	else                      { local log "quietly" }

	
	qui {
		local  lfprog  "mzip_llf"
		
		forvalues k=1/3 { local ll`k' = . }
			
		if "`nc'" == "" & "`lrtest'" == "" {
			/***********/
			/* Model 1 */
			/***********/
			`lll' di in gr _n "Fitting constant-only zero-inflated Poisson model:"
			
			global SGLM_m = 1
			#delimit ;
			`vv'
			`lll' ml model `mlmetho' `lfprog'
				(`dep': `dep' = , `offo')
				(inflate: `inflate', `nc2' `off2a')
				if `touse' `wtexp', wald(0) 
				collinear missing max nooutput nopreserve 
				`mlopt' search(on) `i1aopt' 
				title("Constant-only model:")
				nocnsnotes `negh' iter(8) nrtolerance(1e-1) ltolerance(1e-1) ;
			#delimit cr
			if e(converged)==1 {
				local ll1 = e(ll)
				local llarg "lf0(1 `e(ll)')"
			}
			if "`from'"=="" {
				tempname from1
				mat `from1' = get(_b)  // Keep for starting values in full model
				local initopt "init(`from1', skip)"
			}
		}
	}

	if "`from'" != "" {
		local initopt "init(`from', copy)"
	}
	
	/***********/
	/* Model 0 */
	/***********/
	global SGLM_m = 1
	`lll' di in gr _n "Fitting full marginalized zero-inflated Poisson model:"
	#delimit ;
	`vv'
	`log' ml model `mlmetho' `lfprog'
		(`dep': `dep' = `ind', `nc' `offo')
		(inflate: `inflate', `nc2' `off2a')
		if `touse' `wtexp',
		collinear missing max nooutput nopreserve `mlopt'
		title(Marginalized Zero-inflated Poisson regression) 
		`scopt' `robust' `clopt' 
		`initopt' search(on) `llarg' 
		`negh' ;
	#delimit cr

	if "`score'" != "" { 
		rename `sc1' `svar1' 
		rename `sc2' `svar2' 
		est local scorevars `svar1' `svar2'
	}
	
	// calculate correct df_m
	local coln : colfullnames e(b)
	local dfm 0
	foreach name of local coln {
		tokenize `name', parse(":")
		if "`1'" == "`dep'" {
			if "`3'" != "_cons" {
				_ms_parse_parts `3'
				// constraints are not marked as omitted
				if !r(omit) {
					local dfm = `dfm' + ///
					cond(_b[`3']==0 & _se[`3']==0,0,1)
				}
			}
		}
	}
	est local cmd
	est scalar N_zero = `nzero'
	est scalar df_m   = `dfm'
	est scalar df_c   = `df'
	if "`llb'" != "" { est scalar llb    = `llb' }

	est local narg    "`narg'"
	est local inflate "`ilnk'"
	est local offset1 "`offstr'"
	est local offset2 "`off2s'"
	est local predict "mzip_p"
	est local ilink   "`ilnk'"
	est local cmd     "mzip"
	
	if "`display'" == "" {
		Display, `diopts'
		_prefix_footnote
	}
end


program define Display
	if "`e(prefix)'" == "svy" {
		_prefix_display `0'
		exit
	}
	syntax [, Level(cilevel) EFORM *]

	_get_diopts diopts, `options'
	DispHd
	di
	version 9: ml di, level(`level') `eform' nohead nofootnote `diopts'

	local fmt "%8.2e"
	if ((e(chi2_c) > 0.005) & (e(chi2_c)<1e4)) { 
		local fmt "%8.2f" 
	}
end

program define DispHd
	if !missing(e(df_r)) {
		local model _col(51) as txt "F(" ///
			as res %3.0f e(df_m) as txt "," ///
			as res %6.0f e(df_r) as txt ")" _col(67) ///
			"=" _col(70) as res %9.2f e(F)
		local pvalue _col(51) as txt "Prob > F" _col(67) ///
			"=" _col(70) as res %9.4f Ftail(e(df_m),e(df_r),e(F))
	}
	else {
		if "`e(chi2type)'" == "" {
			local chitype Wald
		}
		else	local chitype `e(chi2type)'
		local model _col(51) as txt `"`chitype' chi2("' ///
			as res e(df_m) as txt ")" _col(67) ///
			"=" _col(70) as res %9.2f e(chi2)
		local pvalue _col(51) as txt "Prob > chi2" _col(67) ///
			"=" _col(70) as res %9.4f chiprob(e(df_m),e(chi2))
	}

	local crtype = upper(substr(`"`e(crittype)'"',1,1)) + /*
		*/ substr(`"`e(crittype)'"',2,.)

	local infl "" 
	local crlen = max(length("`infl'"),length(`"`crtype'"')) 
	di in gr _n "`e(title)'" _col(51) "Number of obs   =" ///
		in ye _col(70) %9.0g e(N) 
	di in gr _col(51) "Nonzero obs" _col(67) "=" ///
		in ye _col(70) %9.0g e(N) - e(N_zero) 
	di in gr "Inflation link : " as res "`e(ilink)'" as txt _col(51) "Zero obs" _col(67) "=" ///
		in ye _col(70) %9.0g e(N_zero) 
		
	di in gr  `model' 
	di in gr %-`crlen's "`crtype'" " = " in ye %9.0g e(ll) `pvalue' 
end


program define PrseOff, sclass
	local 0 : subinstr local 0 "_cons" "", word
	local 0 : subinstr local 0 "_cons," ",", word
	syntax [varlist(default=none fv)] [, OFFset(varname numeric) noCONstant]
	local fvops = "`s(fvops)'" == "true"
	sret clear
	sret local inflate "`varlist'"
	sret local off2    "`offset'"
	sret local nc2     "`constan'"
	sret local fvops = cond(`fvops', "true", "")
end

program define GetLink, rclass
	args ilink
	foreach mac in ilink {
		if "``mac''" == "" { local `mac' "logit" }
	}
	local ilinkrc = 1
	foreach mac in ilink {
		foreach lnk in log logit probit loglog cloglog {
			if "``mac''" == "`lnk'" { // substr("`lnk'", 1, max(3, length("``mac''"))) {
				local `mac' = "`lnk'"
				local `mac'rc = 0
			}
		}
	}
	if "`ilink'" == "log"      { global MLINK = 3 }
	if "`ilink'" == "logit"    { global MLINK = 2 }
	if "`ilink'" == "probit"   { global MLINK = 8 }
	if "`ilink'" == "loglog"   { global MLINK = 6 }
	if "`ilink'" == "cloglog"  { global MLINK = 7 }

	ret scalar rc = `ilinkrc'
	ret local ilink "`ilink'"
end

exit
