from __future__ import print_function

import importlib

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
from pdb import set_trace as keyboard
import compressible.BC as BC
import compressible.eos as eos
import compressible.derives as derives
import mesh.boundary as bnd
import mesh.patch as patch
from simulation_null import NullSimulation, grid_setup, bc_setup
import compressible.unsplit_fluxes as flx
import mesh.integration as integration

class Variables(object):
    """
    a container class for easy access to the different compressible
    variable by an integer key
    """
    def __init__(self, myd):
        self.nvar = len(myd.names)

        # conserved variables -- we set these when we initialize for
        # they match the CellCenterData2d object
        self.idens = myd.names.index("density")
        self.ixmom = myd.names.index("x-momentum")
        self.iymom = myd.names.index("y-momentum")
        self.iener = myd.names.index("energy")

        # if there are any additional variable, we treat them as
        # passively advected scalars
        self.naux = self.nvar - 4
        if self.naux > 0:
            self.irhox = 4
        else:
            self.irhox = -1

        # primitive variables
        self.nq = 4 + self.naux

        self.irho = 0
        self.iu = 1
        self.iv = 2
        self.ip = 3

        if self.naux > 0:
            self.ix = 4   # advected scalar
        else:
            self.ix = -1


def cons_to_prim(U, gamma, ivars, myg):
    """ convert an input vector of conserved variables to primitive variables """
    
    q = myg.scratch_array(nvar=ivars.nq)

    q[:,:,ivars.irho] = U[:,:,ivars.idens]
    q[:,:,ivars.iu] = U[:,:,ivars.ixmom]/U[:,:,ivars.idens]
    q[:,:,ivars.iv] = U[:,:,ivars.iymom]/U[:,:,ivars.idens]
    q[:,:,ivars.iv] = 0.0*q[:,:,ivars.iv]

    ## here iener will be in j/kg and it is rho*E
    e = (U[:,:,ivars.iener] - 0.5*q[:,:,ivars.irho]*(q[:,:,ivars.iu]**2 +  q[:,:,ivars.iv]**2))/q[:,:,ivars.irho] #why do we need rho here ?
    #keyboard()
    #q[:,:,ivars.iener] = e
    q[:,:,ivars.ip] = eos.pres(q[:,:,ivars.irho], e)
    if ivars.naux > 0:
        q[:,:,ivars.ix:ivars.ix+ivars.naux] = \
            U[:,:,ivars.irhox:ivars+naux]/q[:,:,ivars.irho]

    return q


def prim_to_cons(q, gamma, ivars, myg):
    """ convert an input vector of primitive variables to conserved variables """
    
    U = myg.scratch_array(nvar=ivars.nvar)

    U[:,:,ivars.idens] = q[:,:,ivars.irho] 
    U[:,:,ivars.ixmom] = q[:,:,ivars.iu]*U[:,:,ivars.idens]
    U[:,:,ivars.iymom] = q[:,:,ivars.iv]*U[:,:,ivars.idens]
    
    #rhoe = eos.rhoe(gamma, q[:,:,ivars.ip])
    rhoe = eos.rhoe(U[:,:,ivars.idens], q[:,:,ivars.ip])

    U[:,:,ivars.iener] = rhoe + 0.5*q[:,:,ivars.irho]*(q[:,:,ivars.iu]**2 + 
                                                       q[:,:,ivars.iv]**2)

    if ivars.naux > 0:
        U[:,:,ivars.irhox] = q[:,:,ivars.ix]*q[:,:,ivars.irho] 

    return U


class Simulation(NullSimulation):

    def initialize(self, extra_vars=None):
        """
        Initialize the grid and variables for compressible flow and set
        the initial conditions for the chosen problem.
        """
        my_grid = grid_setup(self.rp, ng=4)
        my_data = patch.CellCenterData2d(my_grid)

        # define solver specific boundary condition routines
        bnd.define_bc("hse", BC.user, is_solid=False)

        bc, bc_xodd, bc_yodd = bc_setup(self.rp)

        # are we dealing with solid boundaries? we'll use these for
        # the Riemann solver
        self.solid = bnd.bc_is_solid(self.rp)

        # density and energy
        my_data.register_var("density", bc)
        my_data.register_var("energy", bc)
        my_data.register_var("x-momentum", bc_xodd)
        my_data.register_var("y-momentum", bc_yodd)

        # any extras?
        if extra_vars is not None:
            for v in extra_vars:
                my_data.register_var(v, bc)

        # store the EOS gamma as an auxillary quantity so we can have a
        # self-contained object stored in output files to make plots.
        # store grav because we'll need that in some BCs
        my_data.set_aux("gamma", self.rp.get_param("eos.gamma"))
        my_data.set_aux("grav", self.rp.get_param("compressible.grav"))

        my_data.create()

        self.cc_data = my_data

        # some auxillary data that we'll need to fill GC in, but isn't
        # really part of the main solution
        aux_data = patch.CellCenterData2d(my_grid)
        aux_data.register_var("ymom_src", bc_yodd)
        aux_data.register_var("E_src", bc)
        aux_data.create()
        self.aux_data = aux_data

        self.ivars = Variables(my_data)

        # derived variables
        self.cc_data.add_derived(derives.derive_primitives)

        # initial conditions for the problem
        problem = importlib.import_module("{}.problems.{}".format(
            self.solver_name, self.problem_name))
        problem.init_data(self.cc_data, self.rp)

        if self.verbose > 0: print(my_data)


    def method_compute_timestep(self):
        """
        The timestep function computes the advective timestep (CFL)
        constraint.  The CFL constraint says that information cannot
        propagate further than one zone per timestep.

        We use the driver.cfl parameter to control what fraction of the
        CFL step we actually take.
        """

        cfl = self.rp.get_param("driver.cfl")

        # get the variables we need
        u, v, cs = self.cc_data.get_var(["velocity", "soundspeed"])

        # the timestep is min(dx/(|u| + cs), dy/(|v| + cs))
        xtmp = self.cc_data.grid.dx/(abs(u) + cs)
        ytmp = self.cc_data.grid.dy/(abs(v) + cs)

        self.dt = cfl*float(min(xtmp.min(), ytmp.min()))

    def evolve(self):
        """
        Evolve the equations of compressible hydrodynamics through a
        timestep dt.
        """

        tm_evolve = self.tc.timer("evolve")
        tm_evolve.begin()

        dens = self.cc_data.get_var("density")
        ymom = self.cc_data.get_var("y-momentum")
        ener = self.cc_data.get_var("energy")

        grav = self.rp.get_param("compressible.grav")

        # myg = self.cc_data.grid

        # ### ---- Normal time update ---- ####
        # Flux_x, Flux_y = flx.unsplit_fluxes(self.cc_data, self.aux_data, self.rp,
        #                                     self.ivars, self.solid, self.tc, self.dt)

        # old_dens = dens.copy()
        # old_ymom = ymom.copy()

        # # conservative update
        # dtdx = self.dt/myg.dx
        # dtdy = self.dt/myg.dy

        # for n in range(self.ivars.nvar):
        #     var = self.cc_data.get_var_by_index(n)

        #     var.v()[:,:] += \
        #         dtdx*(Flux_x.v(n=n) - Flux_x.ip(1, n=n)) + \
        #         dtdy*(Flux_y.v(n=n) - Flux_y.jp(1, n=n))

        # # gravitational source terms
        # ymom[:,:] += 0.5*self.dt*(dens[:,:] + old_dens[:,:])*grav
        # ener[:,:] += 0.5*self.dt*(ymom[:,:] + old_ymom[:,:])*grav

        myd = self.cc_data

        method = self.rp.get_param("compressible.temporal_method")

        rk = integration.RKIntegrator(myd.t, self.dt, method=method)
        rk.set_start(myd)

        for s in range(rk.nstages()):
            ytmp = rk.get_stage_start(s)
            ytmp.fill_BC_all()
            k = self.substep(ytmp)
            rk.store_increment(s, k)

        rk.compute_final_update()

        # increment the time
        self.cc_data.t += self.dt
        self.n += 1

        tm_evolve.end()


    def dovis(self):
        """
        Do runtime visualization.
        """

        plt.clf()
        plt.rc("font", size=10)

        # we do this even though ivars is in self, so this works when
        # we are plotting from a file
        ivars = Variables(self.cc_data)

        # access gamma from the cc_data object so we can use dovis
        # outside of a running simulation.
        gamma = self.cc_data.get_aux("gamma")

        q = cons_to_prim(self.cc_data.data, gamma, ivars, self.cc_data.grid)

        rho = q[:,:,ivars.irho]
        u = q[:,:,ivars.iu]
        v = q[:,:,ivars.iv]
        p = q[:,:,ivars.ip]
        #e = q[:,:,ivars.iener]
        e = eos.rhoe(rho, p)/rho
        #keyboard()

        magvel = np.sqrt(u**2 + v**2)

        myg = self.cc_data.grid

        # figure out the geometry
        L_x = self.cc_data.grid.xmax - self.cc_data.grid.xmin
        L_y = self.cc_data.grid.ymax - self.cc_data.grid.ymin

        f = plt.figure(1)

        cbar_title = False

        if L_x > 2*L_y:
            # we want 4 rows:
            axes = AxesGrid(f, 111,
                            nrows_ncols=(4, 1),
                            share_all=True,
                            cbar_mode="each",
                            cbar_location="top",
                            cbar_pad="10%",
                            cbar_size="25%",
                            axes_pad=(0.25, 0.65),
                            add_all=True, label_mode="L")
            cbar_title = True

        elif L_y > 2*L_x:
            # we want 4 columns:  rho  |U|  p  e
            axes = AxesGrid(f, 111,
                            nrows_ncols=(1, 4),
                            share_all=True,
                            cbar_mode="each",
                            cbar_location="right",
                            cbar_pad="10%",
                            cbar_size="25%",
                            axes_pad=(0.65, 0.25),
                            add_all=True, label_mode="L")

        else:
            # 2x2 grid of plots
            axes = AxesGrid(f, 111,
                            nrows_ncols=(2, 2),
                            share_all=True,
                            cbar_mode="each",
                            cbar_location="right",
                            cbar_pad="2%",
                            axes_pad=(0.65, 0.25),
                            add_all=True, label_mode="L")

        fields = [rho, magvel, p, e]
        field_names = [r"$\rho$", r"U", "p", "e"]

        for n in range(4):
            ax = axes[n]

            v = fields[n]

            img = ax.imshow(np.transpose(v.v()),
                            interpolation="nearest", origin="lower",
                            extent=[myg.xmin, myg.xmax, myg.ymin, myg.ymax],
                            cmap=self.cm)

            ax.set_xlabel("x")
            ax.set_ylabel("y")
            
            # needed for PDF rendering
            cb = axes.cbar_axes[n].colorbar(img)
            cb.solids.set_rasterized(True)
            cb.solids.set_edgecolor("face")

            if cbar_title:
                cb.ax.set_title(field_names[n])
            else:
                ax.set_title(field_names[n])

        if self.cc_data.t > 2.745E-04 :
            keyboard()
        plt.figtext(0.05, 0.0125, "t = {:10.5g}".format(self.cc_data.t))

        plt.pause(0.001)
        plt.draw()

    def substep(self, myd):
        """
        take a single substep in the RK timestepping starting with the 
        conservative state defined as part of myd
        """

        myg = myd.grid
        grav = self.rp.get_param("compressible.grav")

        # compute the source terms
        dens = myd.get_var("density")
        ymom = myd.get_var("y-momentum")

        ymom_src = myg.scratch_array()
        ymom_src.v()[:,:] = dens.v()[:,:]*grav

        E_src = myg.scratch_array()
        E_src.v()[:,:] = ymom.v()[:,:]*grav

        k = myg.scratch_array(nvar=self.ivars.nvar)

        flux_x, flux_y = flx.unsplit_fluxes(self.cc_data, self.aux_data, self.rp,
                                            self.ivars, self.solid, self.tc, self.dt)

        for n in range(self.ivars.nvar):
            k.v(n=n)[:,:] = \
               (flux_x.v(n=n) - flux_x.ip(1, n=n))/myg.dx + \
               (flux_y.v(n=n) - flux_y.jp(1, n=n))/myg.dy

        k.v(n=self.ivars.iymom)[:,:] += ymom_src.v()[:,:]
        k.v(n=self.ivars.iener)[:,:] += E_src.v()[:,:]

        return k