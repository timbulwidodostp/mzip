{smcl}
{* *! version 1.0.0  26jan2015}{...}
{viewerjumpto "Syntax" "mzigp##syntax"}{...}
{viewerjumpto "Description" "mzigp##description"}{...}
{viewerjumpto "Options" "mzigp##options"}{...}
{viewerjumpto "Examples" "mzigp##examples"}{...}
{viewerjumpto "Stored results" "mzigp##results"}{...}
{viewerjumpto "Author" "mzigp##author"}{...}
{viewerjumpto "Also see" "mzigp##alsosee"}{...}
{cmd:help mzigp}{right: ({browse "https://doi.org/10.1177/1536867X19874209":SJ19-3: st0563})}
{hline}

{title:Title}

{p2colset 5 14 16 2}{...}
{p2col :{cmd:mzigp} {hline 2}}Marginalized zero-inflated generalized Poisson regression{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 14 2}
{cmd:mzigp} {depvar} [{indepvars}] {ifin}
[{it:{help mzigp##weight:weight}}]{cmd:,}{break} 
   {opt inf:late}{cmd:(}{varlist}[{cmd:,} {opth off:set(varname)}]|{cmd:_cons)}
   [{it:options}]

{synoptset 28 tabbed}{...}
{synopthdr:options}
{synoptline}
{syntab:Model}
{p2coldent :* {opt inf:late()}}equation that determines whether the count is
zero{p_end}
{synopt :{opt nocons:tant}}suppress constant term{p_end}
{synopt :{opth exp:osure(varname:varname_e)}}include {opt ln(varname_e)} in model with
coefficient constrained to 1{p_end}
{synopt :{opth off:set(varname:varname_o)}}include {it:varname_o} in model with
coefficient constrained to 1{p_end}
{synopt :{cmdab:const:raints(}{it:{help estimation options##constraints():constraints}}{cmd:)}}apply specified linear constraints{p_end}
{synopt:{opt col:linear}}keep collinear variables{p_end}
{synopt :{opth il:ink(mzigp##linkname:linkname)}}link function; default is
{cmd:ilink(logit)}{p_end}

{syntab:SE/Robust}
{synopt :{opth vce(vcetype)}}{it:vcetype} may be {opt oim},
  {opt r:obust}, {opt cl:uster} {it:clustvar}, {opt opg}, {opt boot:strap},
  or {opt jack:knife}{p_end}

{syntab:Reporting}
{synopt :{opt l:evel(#)}}set confidence level; default is {cmd:level(95)}{p_end}
{synopt :{opt eform}}report exponentiated coefficients{p_end}
{synopt :{opt nocnsr:eport}}do not display constraints{p_end}
{synopt :{it:{help mzigp##display_options:display_options}}}control
INCLUDE help shortdes-displayoptall

{syntab:Maximization}
{synopt :{it:{help mzigp##maximize_options:maximize_options}}}control the maximization process; seldom used{p_end}
INCLUDE help shortdes-coeflegend
{synoptline}

{marker linkname}{...}
{synoptset 23}{...}
{synopthdr :linkname}
{synoptline}
{synopt :{opt l:ogit}}logit{p_end}
{synopt :{opt p:robit}}probit{p_end}
{synopt :{opt c:loglog}}cloglog{p_end}
{synopt :{opt logl:og}}loglog{p_end}
{synoptline}

{p2colreset}{...}
{p 4 6 2}
* {opt inf:late}{cmd:(}{it:varlist}[{cmd:,} {opt off:set(varname)}]|{cmd:_cons)}
is required.{p_end}
INCLUDE help fvvarlist2
{p 4 6 2}{cmd:bootstrap}, {cmd:by}, {cmd:fp}, {cmd:jackknife}, {cmd:rolling},
{cmd:statsby}, and {cmd:svy} are allowed; see {help prefix}.{p_end}
{p 4 6 2}Weights are not allowed with the {helpb bootstrap} prefix.{p_end}
{p 4 6 2}
{opt vce()} and weights are not allowed with the {helpb svy}
prefix.{p_end}
{marker weight}{...}
{p 4 6 2}{cmd:fweight}s, {cmd:iweight}s, and {cmd:pweight}s are allowed; see
{help weight}.{p_end}
{p 4 6 2}See {helpb mzigp postestimation} for features
available after estimation.{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:mzigp} estimates a marginalized zero-inflated generalized Poisson
regression of {depvar} on {indepvars}, where {it:depvar} is a nonnegative
count variable.


{marker options}{...}
{title:Options}

{dlgtab:Model}

{phang}
{cmd:inflate(}{varlist}[{cmd:,} {cmd:offset(}{varname}{cmd:)}]|{cmd:_cons)}
specifies the equation that determines whether the observed count is zero.
{opt inflate()} is required.  Conceptually, omitting {opt inflate()} would be
equivalent to fitting the model with {helpb poisson}.

{pmore}
{cmd:inflate(}{it:varlist}[{cmd:, offset(}{it:varname}{cmd:)}]{cmd:)}
specifies the variables in the equation.  You may optionally include an offset
for this {it:varlist}.

{pmore}
{cmd:inflate(_cons)} specifies that the equation determining whether the count
is zero contains only an intercept.  To run a zero-inflated model of {depvar}
with only an intercept in both equations, type
{bind:{cmd:mzigp} {it:depvar}{cmd:,} {cmd:inflate(_cons)}}.

{phang}
{opt noconstant}, {opth exposure:(varname:varname_e)}, {opt offset(varname_o)},
{opt constraints(constraints)}, {opt collinear}; see
{helpb estimation options:[R] estimation options}.

{phang}
{opth ilink:(mzigp##linkname:linkname)} specifies the link function.  The
default is {cmd:ilink(logit)}.

{dlgtab:SE/Robust}

INCLUDE help vce_asymptall

{dlgtab:Reporting}

{phang}
{opt level(#)}; see {helpb estimation options##level():[R] estimation options}.

{phang}
{opt eform} reports estimated coefficients transformed to incidence-rate
ratios.  Standard errors and confidence intervals are similarly transformed.
This option affects how results are displayed, not how they are estimated or
stored.  {opt eform} may be specified at estimation or when replaying
previously estimated results.

{phang}
{opt nocnsreport}; see
{helpb estimation options##nocnsreport:[R] estimation options}.

INCLUDE help displayopts_list

{marker maximize_options}{...}
{dlgtab:Maximization}

{phang}
{it:maximize_options}: {opt dif:ficult},
{opth tech:nique(maximize##algorithm_spec:algorithm_spec)}, 
{opt iter:ate(#)}, [{cmd:{ul:no}}]{cmd:{ul:lo}}{cmd:g}, {opt tr:ace}, 
{opt grad:ient}, {opt showstep}, {opt hess:ian}, {opt showtol:erance},
{opt tol:erance(#)}, {opt ltol:erance(#)},
{opt nrtol:erance(#)}, {opt nonrtol:erance}, and {opt from(init_specs)}; see
{helpb maximize:[R] maximize}.  These options are seldom used.

{pmore}
Setting the optimization type to {cmd:technique(bhhh)} resets the default
{it:vcetype} to {cmd:vce(opg)}.

{phang}
{opt coeflegend}; see
{helpb estimation options##coeflegend:[R] estimation options}.


{marker examples}{...}
{title:Example}

{pstd}Setup{p_end}
{phang2}{cmd:. webuse fish}{p_end}

{pstd}Fit marginalized zero-inflated generalized Poisson model{p_end}
{phang2}{cmd:. mzigp count persons livebait, inflate(child) ilink(loglog) nolog}{p_end}

{pstd}Replay results, displaying exponentiated coefficients{p_end}
{phang2}{cmd:. mzigp, eform}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:mzigp} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(N_zero)}}number of zero observations{p_end}
{synopt:{cmd:e(k)}}number of parameters{p_end}
{synopt:{cmd:e(k_eq)}}number of equations in {cmd:e(b)}{p_end}
{synopt:{cmd:e(k_eq_model)}}number of equations in overall model test{p_end}
{synopt:{cmd:e(k_dv)}}number of dependent variables{p_end}
{synopt:{cmd:e(df_m)}}model degrees of freedom{p_end}
{synopt:{cmd:e(ll)}}log likelihood{p_end}
{synopt:{cmd:e(ll_0)}}log likelihood, constant-only model{p_end}
{synopt:{cmd:e(ll_c)}}log likelihood, comparison model{p_end}
{synopt:{cmd:e(df_c)}}degrees of freedom for comparison test{p_end}
{synopt:{cmd:e(chi2)}}chi-squared{p_end}
{synopt:{cmd:e(p)}}significance of model test{p_end}
{synopt:{cmd:e(rank)}}rank of {cmd:e(V)}{p_end}
{synopt:{cmd:e(ic)}}number of iterations{p_end}
{synopt:{cmd:e(rc)}}return code{p_end}
{synopt:{cmd:e(converged)}}{cmd:1} if converged, {cmd:0} otherwise{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:mzigp}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(inflate)}}{cmd:logit} or {cmd:probit}{p_end}
{synopt:{cmd:e(ilink)}}link function{p_end}
{synopt:{cmd:e(title)}}title in estimation output{p_end}
{synopt:{cmd:e(chi2type)}}{cmd:Wald} or {cmd:LR}; type of model chi-squared
	test{p_end}
{synopt:{cmd:e(vce)}}{it:vcetype} specified in {cmd:vce()}{p_end}
{synopt:{cmd:e(opt)}}type of optimization{p_end}
{synopt:{cmd:e(which)}}{cmd:max} or {cmd:min}; whether optimizer is to perform
                         maximization or minimization{p_end}
{synopt:{cmd:e(ml_method)}}type of {cmd:ml} method{p_end}
{synopt:{cmd:e(user)}}name of likelihood-evaluator program{p_end}
{synopt:{cmd:e(technique)}}maximization technique{p_end}
{synopt:{cmd:e(properties)}}{cmd:b V}{p_end}
{synopt:{cmd:e(predict)}}program used to implement {cmd:predict}{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(ilog)}}iteration log (up to 20 iterations){p_end}
{synopt:{cmd:e(gradient)}}gradient vector{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix of the estimators{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}
{p2colreset}{...}


{marker author}{...}
{title:Author}

{pstd}James W. Hardin{p_end}
{pstd}Department of Epidemiology and Biostatistics{p_end}
{pstd}University of South Carolina{p_end}
{pstd}Columbia, SC{p_end}
{pstd}jhardin@sc.edu{p_end}


{marker alsosee}{...}
{title:Also see}

{p 4 14 2}
Article:  {it:Stata Journal}, volume 19, number 3: {browse "https://doi.org/10.1177/1536867X19874209":st0563}{p_end}

{p 7 14 2}
Help:  {helpb mzigp postestimation}, {helpb zigp} (if installed),
{manhelp svy_estimation SVY:svy estimation}{p_end}
