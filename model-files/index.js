var util = require('util')
  , fs = require('fs')
  , clone = require('clone')
  , path =require('path')
  , _ =require('underscore')
  , contexts = require('./lib/contexts')
  , blocks = require('./lib/blocks')
  , add_noise = require('./lib/add_noise')
  , links = require('./lib/links')
  , theta = require('./lib/theta')
  , designs = require('./lib/designs')
  , Erlang = require('plom-erlangify');

/**
 * All the models
 **/

var models = {};

models.sir = {
  name: "SIR_HBRS",
  description: "SIR History based model with reduced susceptibility",
  state: ['SS', 'IS', 'SI', 'RS', 'SR', 'IR', 'RI', 'RR'],
  parameter: ['r0_1', 'r0_2', 'v', 'e', 'd', 'sigma', 'sto', 'iota_1', 'iota_2', 'mu_b', 'mu_d'],
  blocks: ['demography', 'infection', 'recovery'],
  link: links.sir,
  design: [designs.sir_simplex, designs.ksimplex, designs.kmcmc, designs.pmcmc_ode, designs.kmcmc_cov, designs.pmcmc_ode_cov]
};

models.sirs = {
  name: "SIRS_HBRS",
  description: "SIRS History based model with reduced susceptibility",
  state: models.sir.state,
  parameter: models.sir.parameter.concat('g'),
  blocks: models.sir.blocks.concat('waning_immunity'),
  link: links.sir,
  design: [designs.sir_simplex, designs.ksimplex, designs.kmcmc, designs.pmcmc_ode, designs.kmcmc_cov, designs.pmcmc_ode_cov]
};

models.siri = {
  name: "SIRI_HBRS",
  description: "SIRI History based model with reduced susceptibility",
  state: models.sir.state,
  parameter: models.sir.parameter.concat('z'),
  blocks: models.sir.blocks.concat('reinfection'),
  link: links.sir,
  design: [designs.siri_simplex, designs.ksimplex, designs.kmcmc, designs.pmcmc_ode, designs.kmcmc_cov, designs.pmcmc_ode_cov]
};

models.siqr = {
  name: "SIQR_HBRS",
  description: "SIQR History based model with reduced susceptibility!",
  state: models.sir.state.concat(["SQ", "QS", "QR", "RQ"]),
  parameter: models.sir.parameter.concat('q'),
  blocks: ['demography', 'infection', 'recovery_Q', 'waning_Q'],
  link: links.siqr,
  design: [designs.sir_simplex, designs.ksimplex, designs.kmcmc, designs.pmcmc_ode, designs.kmcmc_cov, designs.pmcmc_ode_cov]
};

models.siqrs = {
  name: "SIQRS_HBRS",
  description: "SIQRS History based model with reduced susceptibility",
  state: models.siqr.state,
  parameter: models.siqr.parameter.concat('g'),
  blocks: models.siqr.blocks.concat('waning_immunity', 'waning_immunity_Q'),
  link: links.siqr,
  design: [designs.sir_simplex, designs.ksimplex, designs.kmcmc, designs.pmcmc_ode, designs.kmcmc_cov, designs.pmcmc_ode_cov]
};

models.siqri = {
  name: "SIQRI_HBRS",
  description: "SIQRI History based model with reduced susceptibility",
  state: models.siqr.state,
  parameter: models.siqr.parameter.concat('z'),
  blocks: models.siqr.blocks.concat('reinfection'),
  link: links.siqr,
  design: [designs.siri_simplex, designs.ksimplex, designs.kmcmc, designs.pmcmc_ode, designs.kmcmc_cov, designs.pmcmc_ode_cov]
};

models.siqr_b = {
  name: "SIQR_HBRS_B",
  description: "SIQR (with boosting) History based model with reduced susceptibility!",
  state: models.siqr.state,
  parameter: models.siqr.parameter,
  blocks: models.siqr.blocks.concat('boosting_Q'),
  link: links.siqr,
  design: [designs.sir_simplex, designs.ksimplex, designs.kmcmc, designs.pmcmc_ode, designs.kmcmc_cov, designs.pmcmc_ode_cov]
};

//models.siqrs_b = {
//  name: "SIQRS_HBRS_B",
//  description: "SIQRS (with boosting) History based model with reduced susceptibility",
//  state: models.siqrs.state,
//  parameter: models.siqrs.parameter,
//  blocks: models.siqrs.blocks.concat('boosting_Q'),
//  link: links.siqr,
//  design: [designs.sir_simplex, designs.ksimplex, designs.kmcmc, designs.pmcmc_ode, designs.kmcmc_cov, designs.pmcmc_ode_cov]
//};

models.siqri_b = {
  name: "SIQRI_HBRS_B",
  description: "SIQRI (with boosting) History based model with reduced susceptibility",
  state: models.siqri.state,
  parameter: models.siqri.parameter,
  blocks: models.siqri.blocks.concat('boosting_Q_with_reinfection'),
  link: links.siqr,
  design: [designs.siri_simplex, designs.ksimplex, designs.kmcmc, designs.pmcmc_ode, designs.kmcmc_cov, designs.pmcmc_ode_cov]
};

models.sir_sbri = {
  name: "SIR_SBRI",
  description: "SIR Status based model with reduced infectivity",
  state: models.sir.state,
  parameter: models.sir.parameter,
  blocks: ['demography', 'infection_sbri', 'recovery'],
  link: links.sir,
  design: [designs.sir_simplex, designs.ksimplex, designs.kmcmc, designs.pmcmc_ode, designs.kmcmc_cov, designs.pmcmc_ode_cov]
};

models.sirs_sbri = {
  name: "SIRS_SBRI",
  description: "SIRS Status based model with reduced infectivity",
  state: models.sir.state,
  parameter: models.sir.parameter.concat('g'),
  blocks: models.sir_sbri.blocks.concat('waning_immunity'),
  link: links.sir,
  design: [designs.sir_simplex, designs.ksimplex, designs.kmcmc, designs.pmcmc_ode, designs.kmcmc_cov, designs.pmcmc_ode_cov]
};


////Models with an Exposed class
//function addExposed(blockList){
//
//  //replace key by values
//  var adapter = {
//    infection: ['infection_E', 'recovery_E'],
//    infection_sbri: ['infection_sbri_E', 'recovery_E'],
//    reinfection: ['reinfection_E'],
//    waning_immunity: ['waning_immunity', 'waning_immunity_E']
//  };
//
//  var a = clone(blockList);
//
//  a.forEach(function(block){
//    if (block in adapter){
//      a.splice(a.indexOf(block), 1, adapter[block]);
//    }
//  });
//
//  return _.flatten(a);
//}
//
//for (var m in models){
//  var e = m.replace(/i/,'ei');
//
//  models[e] = {
//    name: models[m].name.replace(/I/, 'EI'),
//    description: models[m].description.replace(/I/, 'EI'),
//    state: models[m].state.concat('ES', 'SE', 'ER', 'RE'),
//    parameter: models[m].parameter.concat('l'),
//    blocks: addExposed(models[m].blocks),
//    link: models[m].link,
//    design: models[m].design
//  };
//}
//
//

//add model with Erlang distributed duration of infections
//function addErlang(blockList, shape){
//
//  //replace key by values
//  var adapter = {
//    recovery: 'erlang_I_' + shape,
//    recovery_E: 'erlang_E_' + shape,
//    recovery_Q: 'erlang_I_' + shape
//  };
//
//  var erlang = [];
//  blockList.forEach(function(block){
//    if (block in adapter){
//      erlang.push(adapter[block]);
//    }
//  });
//
//  return erlang;
//}
//
//var nonErlang = Object.keys(models);
//
//[2,3].forEach(function(shape){
//
//  nonErlang.forEach(function(m){
//    var e = m + '_' + shape;
//
//    models[e] = clone(models[m]);
//    models[e].name +=  '_' + shape;
//    models[e].description +=  util.format(' gamma (%d)', shape);
//    models[e].erlang = addErlang(models[m].blocks, shape);
//
//  });
//
//});


////add models with maternal antibodies
//for (var m in models){
//  var tpl = models[m];
//
//  models[m + '_m'] = {
//    name: tpl.name + '_M',
//    descritpion: tpl.description + ' with maternal antibodies',
//    state: tpl.state.concat('MM'),
//    parameter: tpl.parameter.concat('m'),
//    blocks: tpl.blocks.slice(),
//    link: tpl.link,
//    design: tpl.design
//  };
//
//  //replace demography by demography_M
//  var index = tpl.blocks.indexOf('demography');
//  if (index !== -1) {
//    models[m + '_m']['blocks'][index] = 'demography_M';
//  }
//}


/**
 * Vaccination models
 */

var vmodels = {};

for(var m in models){

  vmodels[m+'_v'] = {
    name: models[m].name + '_V',
    description: 'Vaccination of EV71 for ' + models[m].description,
    state: models[m].state,
    parameter: models[m].parameter.concat('p'),
    blocks: models[m].blocks.slice(),
    link: links[models[m].link.name + '_simulation'],
    design: [designs.vaccination]
  };

  //replace demography by demography_vaccination
  vmodels[m+'_v'].blocks[models[m].blocks.indexOf('demography')] = 'demography_vaccination';
}


/**
 * Convert to verbose grammar
 **/

for(model in models){

  var m = models[model];

  var c = contexts.fit;
  c.type = "context";

  var p = {
    name: m.name,
    type: "process",
    state: m.state.map(function(x){
      if(x === 'RR'){
        return {id: x,tag: ["remainder"]};
      } else if (['IS', 'IR', 'SI', 'RI'].indexOf(x) !== -1){
        return {id: x, tag: ["infectious"]};
      } else {
        return {id: x};
      }      
    }),
    parameter: m.parameter.map(function(X){return {id: X};}),
    model: blocks.build(m.blocks, m.state)
  };

  p.white_noise = add_noise(p.model);

  var l = m.link;
  l.type = "link";

  var t = {
    name: m.name + " HFMD Japan",
    description: "...",
    type: "theta",
    parameter: {}
  };

  p.state.concat(p.parameter, l.observation[0].parameter).map(function(x){return x.id}).forEach(function(par){
    t.parameter[par] = theta[par];
  });


  //erlangify
  if('erlang' in m) {
    var erlang = new Erlang(blocks.build(m.erlang));
    p = erlang.ify(p);
    l = erlang.ify(l);
    t = erlang.ify(t);
  }

  var mydir = 'hfmd_' + m.name.toLowerCase();

  //write to JSON
  if(!fs.existsSync(mydir)) fs.mkdirSync(mydir);

  fs.writeFileSync(path.join(mydir, 'context.json'), JSON.stringify(c, null, 2));
  fs.writeFileSync(path.join(mydir, 'process.json'), JSON.stringify(p, null, 2));
  fs.writeFileSync(path.join(mydir, 'link.json'), JSON.stringify(l, null, 2));
  fs.writeFileSync(path.join(mydir, 'theta.json'), JSON.stringify(t, null, 2));

  m.design.forEach(function(d){
    fs.writeFileSync(path.join(mydir, 'design_' + d.name + '.json'), JSON.stringify(d, null, 2));
  });

  //fs.writeFileSync(path.join(mydir, 'design_invasion.json'), JSON.stringify(designs.invasion, null, 2));

  models[model] = {
    parent: m.parent,
    process: p,
    link: l,
    theta: t
  };
}
