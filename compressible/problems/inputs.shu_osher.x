# simple inputs files for the unsplit CTU hydro scheme

[driver]
max_steps = 20000
tmax = 0.2
cfl = 0.4

[compressible]
limiter = 1

[io]
basename = shu_x_
dt_out = 0.05

[mesh]
nx = 212
ny = 10
xmin = -5.0
xmax = 5.0
ymax = .05
xlboundary = outflow
xrboundary = outflow

[sod]
direction = x

#--N2--#
dens_left = 192.85
dens_right = 25

#--CO2--#
#dens_left = 348.8
#dens_right = 3.488

#--ideal--#
#dens_left = 1.0
#dens_right = 0.125

#--CH4--#
#dens_left = 307.33
#dens_right = 199.44

u_left = machref
u_right = 0.0

#--N2--#
p_left = 41.332E06
p_right = 4.0E06

#--CO2--#
#p_left = 73.7E06
#p_right = 0.737E06

#--ideal--#
#p_left = 1.0
#p_right = 0.1

#--CH4--#
#p_left = 14.7E06
#p_right = 4.7E06