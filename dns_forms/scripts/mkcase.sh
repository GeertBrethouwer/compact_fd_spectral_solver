#!/bin/bash
# mkcase.sh DIR ALPHA NL_FORM DY_OP SCAL_FORM SCAL_DY TMAX INIT_LINES [NX NY NZ]
d=$1; mkdir -p $d; cd $d
nx=${9:-32}; ny=${10:-65}; nz=${11:-32}
cat > input.nml <<NML
&grid
  nx = $nx, nz = $nz, ny = $ny, lx = 12.566370614359172, lz = 4.1887902047863905
  grid_type = 'tanh', grid_alpha = $2, fd_order = 6
/
&physics
  re = 2800., forcing = 'cfr', ubulk = 1.0, nl_form = '$3', dy_op = '$4', visc_theta = 0.0
/
&time
  nstep = 100000000, tmax = $7, cfl = 1.2, cfl_visc = 1.5
/
&init
  $8
/
&output
  iout_log = 50, iout_stats = 10, iout_save_stats = 100000000, iout_restart = 100000000, tstats_start = $(echo "$7*0.5" | bc)
/
&scalars
  nscal = 2, scal_form = '$5', scal_dy = '$6'
  pr(1) = 0.71, scal_theta(1) = 0.5, sbc_bot(1) = 'D', sval_bot(1) = 1.0, sbc_top(1) = 'D', sval_top(1) = -1.0, sinit(1) = 'linear'
  pr(2) = 1.e10, scal_theta(2) = 0.5, sbc_bot(2) = 'D', sval_bot(2) = 0.0, sbc_top(2) = 'D', sval_top(2) = 0.0, sinit(2) = 'noise'
/
&buoyancy
/
&parallel
  pdims = 0, 0
/
NML
