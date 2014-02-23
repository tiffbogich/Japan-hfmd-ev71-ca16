/**
 * contexts
 **/

var clone = require('clone');

var contexts = module.exports = {};

contexts.fit = {
  name: "Japan",
  description: "Hand, foot, and mouth disease in Japan reported from IASR. Virological and syndromic data",
  disease: ["HFMD"],

  population: [{id: "Japan__all"}],

  time_series: [{id: "EV71_JAP__IASR__inc",
                 population_id: ["Japan__all"]},
                {id: "CA16_JAP__IASR__inc",
                 population_id: ["Japan__all"]}],

  data: {
    comment: "HFMD data Japan (IASR) for EV71, CA16 and HFMD",
    t0: "1999-12-25",
    source: "../data/data.csv"
  },

  metadata: [
    {id: "prop",
     comment: "proportion of the population under surveillance",
     source: "../data/prop.csv"},
    {id: "ptest",
     comment: "proportion of hospitals testing patient with symptoms for EV71 and CA16",
     source: "../data/ptest.csv"},
    {id: "hfmd",
     comment: "syndromic HFMD incidence data used as a covariate",
     source: "../data/hfmd.csv"},
    {id: "N",
     comment: "population size",
     source: "../data/N.csv"},
    {id: "mu_b",
     comment: "birth rates",
     unit: "W",
     source: "../data/mu_b.csv"},
    {id: "mu_d",
     comment: "death rates",
     unit: "W",
     source: "../data/mu_d.csv"}
  ]

};


