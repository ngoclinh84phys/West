#!/bin/bash

${WGET} http://www.quantum-simulation.org/potentials/sg15_oncv/upf/H_ONCV_PBE-1.0.upf
${WGET} http://www.quantum-simulation.org/potentials/sg15_oncv/upf/Si_ONCV_PBE-1.1.upf

cat > pw.in << EOF
&control
calculation  = 'scf'
restart_mode = 'from_scratch'
pseudo_dir   = './'
outdir       = './'
prefix       = 'test'
wf_collect   = .TRUE.
/
&system
ibrav           = 1
celldm(1)       = 20
nat             = 5
ntyp            = 2
ecutwfc         = 25.0
nbnd            = 10
assume_isolated ='mp'
/
&electrons
diago_full_acc = .TRUE.
/
ATOMIC_SPECIES
Si 28.0855  Si_ONCV_PBE-1.1.upf
H  1.00794   H_ONCV_PBE-1.0.upf
ATOMIC_POSITIONS bohr
Si      10.000000   10.000000  10.000000
H       11.614581   11.614581  11.614581
H        8.385418    8.385418  11.614581
H        8.385418   11.614581   8.385418
H       11.614581    8.385418   8.385418
K_POINTS {gamma}
EOF


cat > wstat.in << EOF
{
  "input_west": {
    "qe_prefix": "test",
    "west_prefix": "test",
    "outdir": "./"
  },
  "wstat_control": {
    "wstat_calculation": "S",
    "n_pdep_eigen": 50
  }
}
EOF


cat > wfreq.in << EOF 
{
  "input_west": {
    "qe_prefix": "test",
    "west_prefix": "test",
    "outdir": "./"
  },
  "wstat_control": {
    "wstat_calculation": "S",
    "n_pdep_eigen": 50
  },
  "wfreq_control": {
    "wfreq_calculation": "XWGQ",
    "n_pdep_eigen_to_use": 50,
    "qp_bandrange": [1,5],
    "n_refreq": 300,
    "ecut_refreq": 2.0
  }
}
EOF


cat > westpp.in << EOF 
{
  "input_west": {
    "qe_prefix": "test",
    "west_prefix": "test",
    "outdir": "./"
  },
  "wstat_control": {
    "wstat_calculation": "S",
    "n_pdep_eigen": 50
  },
  "westpp_control" : {
     "westpp_calculation" : "E",
     "westpp_range" : [ 1, 2 ], 
     "westpp_format" : "C", 
     "westpp_sign" : true
  }
}
EOF
