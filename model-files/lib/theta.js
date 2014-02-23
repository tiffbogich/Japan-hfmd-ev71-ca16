module.exports = {
  MM: {min: 0.0, guess: 0.01, max: 0.95, sd_transf: 0.2},

  SS: {min: 0.05, guess: 0.6,  max: 8.0, sd_transf: 0.2, transformation: "scale_pow10_neg"},
  RS: {min: 0.05, guess: 0.6, max: 8.0, sd_transf: 0.2, transformation: "scale_pow10_neg"},
  SR: {min: 0.05, guess: 0.6, max: 8.0, sd_transf: 0.2, transformation: "scale_pow10_neg"},

  IS: {min: 3.0, guess: 5.0, max: 9.0, sd_transf: 0.2, transformation: "scale_pow10_neg"},
  SI: {min: 3.0, guess: 5.0, max: 9.0, sd_transf: 0.2, transformation: "scale_pow10_neg"},

  ES: {min: 0.0, guess: 0.0, max: 0.0, sd_transf: 0.0},
  SE: {min: 0.0, guess: 0.0, max: 0.0, sd_transf: 0.0},
  ER: {min: 0.0, guess: 0.0, max: 0.0, sd_transf: 0.0},
  RE: {min: 0.0, guess: 0.0, max: 0.0, sd_transf: 0.0},
  IR: {min: 0.0, guess: 0.0, max: 0.0, sd_transf: 0.0},
  RI: {min: 0.0, guess: 0.0, max: 0.0, sd_transf: 0.0},
  SQ: {min: 0.0, guess: 0.0, max: 0.0, sd_transf: 0.0},
  QS: {min: 0.0, guess: 0.0, max: 0.0, sd_transf: 0.0},
  RQ: {min: 0.0, guess: 0.0, max: 0.0, sd_transf: 0.0},
  QR: {min: 0.0, guess: 0.0, max: 0.0, sd_transf: 0.0},

  m:      {min: 1.0,   guess: 3.0,   max: 6.0,  sd_transf: 0.04, unit: "M", type: "rate_as_duration"},
  r0_1:   {min: 2.0,   guess: 2.40,  max: 30.0, sd_transf: 0.04},
  r0_2:   {min: 2.0,   guess: 2.40,  max: 30.0, sd_transf: 0.04},
  v:      {min: 2.15,  guess: 2.46,  max: 2.81, sd_transf: 0.04, unit: "D", type: "rate_as_duration", prior:"normal"},
  l:      {min: 0.5,   guess: 3.0,   max: 7.0,  sd_transf: 0.04, unit: "D", type: "rate_as_duration"},
  e:      {min: 0.001 ,guess: 0.1,   max: 0.3,  sd_transf: 0.2,  transformation: "logit"},
  d:      {min: 0.6,   guess: 0.8,   max: 0.99, sd_transf: 0.04},
  sigma:  {min: 0.1,   guess: 0.6,   max: 0.98, sd_transf: 0.2,  transformation: "logit"},
  sto_1:  {min: 0.0,   guess: 0.01,  max: 0.8,  sd_transf: 0.02},
  sto_2:  {min: 0.0,   guess: 0.01,  max: 0.8,  sd_transf: 0.02},
  sto:    {min: 0.0,   guess: 0.05,  max: 0.8,  sd_transf: 0.02},
  iota_1: {min: -3.5,  guess: 0.0,   max: 4.0,  sd_transf: 0.04, transformation: "scale_pow10_bounded"},
  iota_2: {min: -3.5,  guess: 0.0,   max: 4.0,  sd_transf: 0.04, transformation: "scale_pow10_bounded"},

  g: {min: 0.0125, guess: 0.05, max: 1.0,  sd_transf: 0.04, unit: "Y"},
  q: {min: 0.05,   guess: 6.0,  max: 12.0, sd_transf: 0.04, unit: "M", type: "rate_as_duration"},
  z: {min: 0.0,    guess: 0.05, max: 0.95, sd_transf: 0.2,  transformation: "logit"},

  p: {min: 0.0, guess: 0.5, max: 1.0, sd_transf: 0.0, transformation: "logit"},

  rep: {min: 0.05,  guess: 0.175, max: 0.3, sd_transf: 0.2, transformation: "logit", partition_id: "variable_time_series", prior:"normal"}
};
