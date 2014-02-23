module.exports = function(model){

  var white_noise = [
    {
      "reaction": [],
      "sd": "sto"
    },
    {
      "reaction": [],
      "sd": "sto"
    }
  ];

  model.forEach(function(r){
    if (('tag' in r) && (r.tag.indexOf('transmission') !== -1)){

      if (r.rate.indexOf('r0_1') !== -1) {
        white_noise[0].reaction.push({from:r.from, to:r.to, rate:r.rate});
      } else if (r.rate.indexOf('r0_2') !== -1) {
        white_noise[1].reaction.push({from:r.from, to:r.to, rate:r.rate});
      }

    }

  });

  return white_noise;
}
