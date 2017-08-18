# simple inputs files for the unsplit CTU hydro scheme

[driver]
max_steps = 20000
tmax = 0.2

[compressible]
limiter = 1

[io]
basename = sod_x_
dt_out = 0.05

[mesh]
nx = 201
ny = 2

#----Shu-Osher---#
#xmin = -0.5
#xmax = 0.5
#ymax = .05

#----Sod tube---#
xmin = 0.0
xmax = 1.0
ymax = 0.05

xlboundary = outflow
xrboundary = outflow

[sod]
direction = x

#--N2--#
#dens_left = 192.85715
#dens_right = 25

#--CO2--#
#---validation--#
dens_left = 348.8
dens_right = 3.488

#--Super-crit--#
#dens_left = 845.0
#dens_right = 696.0

#--Near-crit---#
#dens_left = 610.0
#dens_right = 134.0

#--ideal--#
#dens_left = 3.857143
#dens_right = 0.125

#--CH4--#
#dens_left = 398.86
#dens_right = 394.08

#u_left = 842.459
#u_right = 0.0

u_left = 0.0
u_right = 0.0

#--N2--#
#p_left = 41.332E06
#p_right = 4.0E06

#--CO2--#
#--validation--#
p_left = 73.7E06
p_right = 0.737E06

#--Supercritical--#
#p_left = 29.264E06
#p_right = 14.632E06

#--near-crit--#
#p_left = 10.974E06
#p_right = 5.487E06

#--CH4--#
#p_left = 4.0E06
#p_right = 0.4E06
