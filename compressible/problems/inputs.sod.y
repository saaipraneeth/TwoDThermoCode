# simple inputs files for the unsplit CTU hydro scheme

[driver]
max_steps = 20000
tmax = 0.2

[io]
basename = sod_y_
dt_out = 0.05

[mesh]
nx = 10
ny = 201
xmax = .05
ymin = -0.5
ymax = 0.5
ylboundary = outflow
yrboundary = outflow

[sod]
direction = y

dens_left = 192.85715
dens_right = 0.125

u_left = 842.459
u_right = 0.0

p_left = 41.332E06
p_right = 4.0E06
