Japan-hfmd-ev71-ca16
====================

Data and analyses for HFMD, EV71, and CA16 in Japan for 2000-2011




Summary of Results
==================


Stochastic results:
==============

The best models are:
1 species: SIR
2 species: SIQRI_HBRS_B and SIQRI_HBRS

The best models are indicated by the DIC plots (For DIC: the smaller
=> the better).

I have added the correlation plots and correlation matrix for all the
parameters. To make things simpler, the modal value of each parameters
(the "best" value) is reported in the diagonal of the correlation
matrix. So looking at the correlation matrix should give a complete
overview of the results. See for instance
corrmat-hfmd_siqri_hbrs_ev_ca.png

To make the matrix readable, I have shorten some parameter names.
rC: reporting rate for CA
rE: reporting rate for EV
r1: r0 for EV
r2: r0 for CA
v_1: duration of infection in days for EV
v_2: duration of infection in days for CA
q: duration of temporary full cross immunity in month
s: sigma: probability of susceptibility reduction in case of cross protection
e: intensity of seasonal forcing
d: dephasing of seasonal forcing
i1: iota_1: import parameter for EV
i2: iota_2: import parameter for CA
z: reinfection probability (SIRI like dynamics)
sto: intensity of noise on the transmission rates
Otherwise for the initial condition EV is always first so SR =>
susceptible for EV but recovered from CA.

As the SIQRI model converged better, it makes sense to look in detail
at the correlation structure of this model.
As an alternative to the matrix view, corr-siqri_hbrs_ev_ca.png gives
100 pairwise correlation plots.

As a side note, the correlation for the single species models are beautiful.

I have also added the ESS (ess.txt) that indicates how many
independent samples we do have in our posteriors.

The credible intervals are indicated by vertical dashed blue lines,
for instance in post_hfmd_siqri_hbrs_b_ev_ca_0.png. The red lines
represent the priors.

Deterministic results
===============

I added the same information for the deterministic model. The best models are:
- 1 species: SIR
- 2 species: SIQRI_HBRS_B

For the deterministic models, if we take the best parameter set of the
posteriors and plot a trajectory, the fit is always really good.