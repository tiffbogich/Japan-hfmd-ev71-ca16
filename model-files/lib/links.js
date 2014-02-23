var clone = require('clone');

var links = module.exports = {
  sir: {
    name: "sir",
    description: "incidence at the recovery time",

    observed: [
      {id: "Inc_EV",
       definition: [{from: "IS", to: "RS"}, {from: "IR", to: "RR"}],
       time_series_id: ["EV71_JAP__IASR__inc"],
       observation_id: "common"},
      {id: "Inc_CA",
       definition: [{from: "SI", to: "SR"}, {from: "RI", to: "RR"}],
       time_series_id: ["CA16_JAP__IASR__inc"],
       observation_id: "common"}
    ],

    observation: [
      {
        id: 'common',
        parameter: [{id: "rep"}],        
        model: {
          distribution: "discretized_normal",
          comment: "We use the syndromic data (hfmd): Y ~ BINOMIAL(hfmd*ptest, x*prop*rep*ptest/(hfmd*ptest)",
          mean: "x*prop*ptest*rep",
          "var": "x*prop*ptest*rep*(1.0-GSL_MIN(x*prop*rep/hfmd, 1.0))"
        }      
      }
    ]
  }
};


//SIQR models
links.siqr = clone(links.sir);
links.siqr.name = 'siqr';
links.siqr.observed = [
  {id: "Inc_EV",
   definition: [{from: "IS", to: "QS"}, {from: "IR", to: "QR"}],
   time_series_id: ["EV71_JAP__IASR__inc"],
   observation_id: "common"},
  {id: "Inc_CA",
   definition: [{from: "SI", to: "SQ"}, {from: "RI", to: "RQ"}],
   time_series_id: ["CA16_JAP__IASR__inc"],
   observation_id: "common"}
];


