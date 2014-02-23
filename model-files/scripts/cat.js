/**
 * Concatenate all the pmcmc_ode runs stored on pubgrenfell
 */

var fs = require('fs')
  , path = require('path')
  , cat = require('plom-cat');

function sortFiles(a, b){
  var i = a.split('_');
  i = i[i.length-1];

  var j = b.split('_');
  j = j[j.length-1];

  return i - j;
}

var directories = fs.readdirSync('.')
  .filter(function(d){return d.indexOf('pmcmc_ode') !== -1})
  .sort(sortFiles);

var models = fs.readdirSync(directories[0])
  .filter(function(m){return m.indexOf('hfmd_') !== -1});

models.forEach(function(m){

  var bestFiles = fs.readdirSync(path.join(directories[0], m, 'model', 'results', 'pmcmc_ode'))
    .filter(function(b){return b.indexOf('best_') !== -1})
    .sort(sortFiles);

  bestFiles.forEach(function(b, h){

    var output = path.join(directories[0], m, 'trace_' + h +'.csv');
    var inputs = [];
    directories.forEach(function(d){
      inputs.push(path.join(d, m, 'model', 'results', 'pmcmc_ode', b));
    });

    cat(inputs, output, {thin:100});

  });

});
