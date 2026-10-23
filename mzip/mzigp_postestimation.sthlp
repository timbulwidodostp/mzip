{smcl}
{* *! version 1.0.0  26jan2015}{...}
{viewerjumpto "Postestimation commands" "mzigp postestimation##description"}{...}
{viewerjumpto "Syntax for predict" "mzigp postestimation##syntax_predict"}{...}
{viewerjumpto "Options for predict" "mzigp postestimation##options_predict"}{...}
{viewerjumpto "Examples" "mzigp postestimation##examples"}{...}
{viewerjumpto "Author" "mzigp postestimation##author"}{...}
{viewerjumpto "Also see" "mzigp postestimation##alsosee"}{...}
{cmd:help mzigp postestimation}{right: ({browse "https://doi.org/10.1177/1536867X19874209":SJ19-3: st0563})}
{hline}

{title:Title}

{p2colset 5 29 31 2}{...}
{p2col :{cmd:mzigp postestimation} {hline 2}}Postestimation tools for mzigp{p_end}
{p2colreset}{...}


{marker description}{...}
{title:Postestimation commands}

{pstd}
The following postestimation commands are available after {cmd:mzigp}:

{synoptset 17 tabbed}{...}
{p2coldent :Command}Description{p_end}
{synoptline}
INCLUDE help post_contrast
INCLUDE help post_estatic
INCLUDE help post_estatsum
INCLUDE help post_estatvce
INCLUDE help post_svy_estat
INCLUDE help post_estimates
INCLUDE help post_forecast_star
INCLUDE help post_lincom
INCLUDE help post_lrtest_star
INCLUDE help post_margins
INCLUDE help post_marginsplot
{synopt :{helpb mzigp postestimation##predict:predict}}predictions, residuals,
influence statistics, and other diagnostic measures{p_end}
INCLUDE help post_predictnl
INCLUDE help post_pwcompare
INCLUDE help post_suest
INCLUDE help post_test
INCLUDE help post_testnl
{synoptline}
{p2colreset}{...}
{phang}
* {cmd:forecast} and {cmd:lrtest} are not appropriate with
{cmd:svy} estimation results.{p_end}


{marker syntax_predict}{...}
{marker predict}{...}
{title:Syntax for predict}

{p 8 16 2}
{cmd:predict} {dtype} {newvar} {ifin}
[{cmd:,} {it:statistic} {opt nooff:set}]

{p 8 16 2}
{cmd:predict} {dtype} {{it:stub}{cmd:*}|{it:{help newvar:newvar_reg}}
   {it:{help newvar:newvar_inflate}}}
   {ifin}{cmd:}

{marker statistic}{...}
{synoptset 17}{...}
{synopthdr:statistic}
{synoptline}
{synopt :{opt n}}number of events; the default{p_end}
{synopt :{opt ir}}incidence rate{p_end}
{synopt :{opt p:r}}probability of a degenerate zero{p_end}
{synopt :{opt xb}}linear prediction{p_end}
{synopt :{opt stdp}}standard error of the linear prediction{p_end}
{synoptline}
{p2colreset}{...}
INCLUDE help esample


{marker options_predict}{...}
{title:Options for predict}

{phang}
{opt n}, the default, calculates the predicted number of events, which is
exp(x_j b) if neither {opt offset()} nor {opt exposure()} was specified when
the model was fit; exp{(x_j b) + offset_j} if {opt offset()} was specified; or
{bind:{exp(x_j b) * exposure_j}} if {opt exposure()} was specified.

{phang}
{opt ir} calculates the incidence rate exp(xb), which is the predicted
number of events when exposure is 1.  This is equivalent to specifying both
the {opt n} and the {opt nooffset} options.

{phang}
{opt pr} calculates the probability Pr(y = 0), where this zero was obtained
from the degenerate distribution F(zg).  If {opt offset()} was specified
within the {opt inflate()} option, then F(zg + offset^g) is calculated.

{phang}
{opt xb} calculates the linear prediction, which is xb if neither
{opt offset()} nor {opt exposure()} was specified; {bind:xb + offset} if
{opt offset()} was specified; or {bind:xb + ln(exposure)} if
{opt exposure()} was specified; see {opt nooffset} below.

{phang}
{opt stdp} calculates the standard error of the linear prediction.

{phang}
{opt nooffset} is relevant only if you specified {opt offset()} or
{cmd:exposure()} when you fit the model.  It modifies the calculations made by
{cmd:predict} so that they ignore the offset or exposure variable; the linear
prediction is treated as xb rather than as {bind:xb + offset} or 
{bind:xb + ln(exposure)}.  Specifying {bind:{cmd:predict} ...{cmd:, nooffset}}
is equivalent to specifying {bind:{cmd:predict} ...{cmd:, ir}}.


{marker examples}{...}
{title:Examples}

{pstd}Setup{p_end}
{phang2}{cmd:. webuse fish}{p_end}
{phang2}{cmd:. mzigp count persons livebait, inflate(child) ilink(loglog) nolog}{p_end}

{pstd}Number of events{p_end}
{phang2}{cmd:. predict n}

{pstd}Linear prediction{p_end}
{phang2}{cmd:. predict lp, xb}

{pstd}Predictive margins for {cmd:livebait}{p_end}
{phang2}{cmd:. margins}

{pstd}Calculate the probability of a zero outcome based on the link function{p_end}
{phang2}{cmd:. predict p, pr}


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
Help:  {helpb mzigp} (if installed){p_end}
