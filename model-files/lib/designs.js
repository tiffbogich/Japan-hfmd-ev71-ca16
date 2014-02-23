var clone = require('clone');

var designs = module.exports = {

  sir_simplex: {
    name: "lhs_simplex",
    description: "with some magic to get the right scaling",

    id: "lhs",
    H: 500,
    seed: "it will converge",
    correlate: [
      {y: "r0_2:all",   x: "r0_1:all",   f: "identity", range: [-1, 1]},
      {y: "iota_2:all", x: "iota_1:all", f: "identity", range: [0, 0.25]}
    ],

    cmd: [
      {
        comment: "Get the initial conditions (no forcing). Note that even if we turn off the forcing, there could be an Hopf bifucation => we compute the average rep on 1000 week",
        pipe: "-D -I -S e:all:guess:0.0",
        algorithm: "simul ode -T 100000 -D 1000 --traj --freq W --quiet"
      },
      {
        comment: "First simplex + rescale reporting rate",
        pipe: "-D -X -r rep",
        algorithm: "simplex -M 10000 --no_trace --prior --quiet"
      },
      {
        comment: "We chain simplex",
        pipe: "-T -u 0.01",
        algorithm: "simplex -M 10000 --no_trace --prior --quiet",
        repeat: 19
      }
    ]
  },

  ksimplex: {
    name: "ksimplex",
    description: "chain 10 ksimplex",

    id: "replicate",
    H: 1,
    seed: "it will converge",

    cmd: [
      {
        comment: "First simplex fixing starting value for sto_1 and sto_2",
        pipe: "-S sto:all:guess:0.05 -u 0.01",
        algorithm: "ksimplex -M 10000 --prior --no_dem_sto --quiet",
        overwrite: false,
      },
      {
        comment: "We chain simplex",
        pipe: "-T -u 0.01",
        algorithm: "ksimplex -M 10000 --prior --no_dem_sto --quiet",
        repeat: 9,
        overwrite: false,
      }
    ]
  },

  pmcmc_ode_cov: {
    name: "pmcmc_ode_cov",
    description: "get a cov matrix pmcmc ode",

    id: "replicate",
    H: 5,
    seed: "it will converge",

    cmd: [
      {
        pipe: "-u 0.01",
        algorithm: "pmcmc ode -J 1 -M 100000 --full --acc --smooth -S 10000000 --quiet",
      },
      {
        pipe: "-C -T -u 0.01",
        algorithm: "pmcmc ode -J 1 -M 100000 --full --acc --smooth -S 10000000 --quiet",
        repeat: 4
      }
    ]
  },

  kmcmc_cov: {
    name: "kmcmc_cov",
    description: "get a cov matrix kmcmc full update",

    id: "replicate",
    H: 5,
    seed: "it will converge",

    cmd: [
      {
        pipe: "-u 0.01",
        algorithm: "kmcmc -M 100000 --full --no_dem_sto --acc --smooth -S 10000000 --quiet"
      },
      {
        pipe: "-C -T -u 0.01",
        algorithm: "kmcmc -M 100000 --full --no_dem_sto --acc --smooth -S 10000000 --quiet",
        repeat: 4
      }
    ]
  },

  pmcmc_ode: {
    name: "pmcmc_ode",
    description: "pmcmc ode",

    id: "replicate",
    H: 5,
    seed: "it will converge",

    cmd: [
      {
        pipe: "-u 0.01",
        algorithm: "pmcmc ode -J 1 -M 1000000 --full --acc --quiet"
      }
    ]
  },

  kmcmc: {
    name: "kmcmc",
    description: "kmcmc full update",

    id: "replicate",
    H: 5,
    seed: "it will converge",

    cmd: [
      {
        pipe: "-u 0.01",
        algorithm: "kmcmc -M 500000 --full --no_dem_sto --acc --quiet"
      }
    ]
  }

//  vaccination: {
//    name: "vaccination",
//    description: "vaccination",
//
//    id: "slice",
//    H: 40,
//    par: ['p:all'],
//
//    cmd: [
//      {
//        comment: "Get the initial conditions (no forcing)",
//        pipe: "-D -S e:all:guess:0.0",
//        algorithm: "simul ode -T 100000 -D 5000 --traj --freq W"
//      }
//    ]
//  },
//
//  invasion: {
//    name: "invasion",
//    description: "invasion",
//
//    id: "slice",
//    H: 20,
//    par: ['z:all'],
//
//    cmd: [
//      {
//        pipe: "-D -S e:all:guess:0.0,r0_1:all:guess:0.0,sigma:all:guess:0.0,iota_2:all:guess:0.1,RS:all:guess:0.0,IR:all:guess:0.0,IS:all:guess:0.0,RI:all:guess:0.0,SS:all:guess:0.5",
//        algorithm: "simul ode -T 100000 -D 1000 --traj --freq W"
//      }
//    ]
//  }

};

designs['siri_simplex'] = clone(designs['sir_simplex']);
designs.siri_simplex.correlate.push({y: "z:all", x: "r0_1:all", f: "inverse", range: [-0.05, 0.05]});
