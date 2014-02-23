HFMD
====

If you end up in a position where you have to fit the HFMD data set
and can't run away, here are some scripts that you mind find useful.
Note that they all have to be run sequentially.

Usage
=====

##dependencies:

    npm install

##generate all the models

    node index.js

##building the models

    ./build.sh

(Edit the file to call a C compiler adapted to your system, you might also want to change the Makefile if you cannot write in ```/usr/local```)


Finding something like a mode
=============================

##LHS simplex

    ./lhs_simplex.sh

##store the results

    ./mle_lhs_simplex.sh

##plots

Note, if you run that on your local machine, you might have to regenerate binaries adapted to your machine

    ./post.sh

will do that for you.

    ./traj_lhs_simplex.sh

in graphs/

    res_lhs_simplex.r


Taking into account environmental stochasticity
===============================================

##ksimplex from the previous mode

    ./ksimplex.sh

##store the results

    ./mle_ksimplex.sh

##plots

    ./traj_ksimplex.sh

in graphs/

    res_ksimplex.r


Getting the posteriors
======================


##get covariance matrix starting from the previous results

    ./pmcmc_ode_cov.sh
    ./kmcmc_cov.sh


##plots

in graphs/

    post_kmcmc_cov.r
    post_pmcm_ode_cov.r


##store the results

    ./mle_pmcmc_ode_cov.sh
    ./mle_kmcmc_cov.sh


##iterate...

Sample using the previously found covariance matrix (you might have to
run that a couple of times, each time using a reduction before
restarting or simply chaining)

    ./pmcmc_ode.sh
    ./kmcmc.sh


Publishing
==========
