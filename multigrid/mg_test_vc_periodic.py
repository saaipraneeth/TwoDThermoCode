#!/usr/bin/env python3

"""
Test the variable-coefficient MG solver with periodic data.

Here we solve:

   div . ( alpha grad phi ) = f

with

   alpha = 2.0 + cos(2.0*pi*x)*cos(2.0*pi*y)

   f = -16.0*pi**2*(cos(2*pi*x)*cos(2*pi*y) + 1)*sin(2*pi*x)*sin(2*pi*y)

This has the exact solution:

   phi = sin(2.0*pi*x)*sin(2.0*pi*y)

on [0,1] x [0,1]

We use Dirichlet BCs on phi.  For alpha, we do not have to impose the
same BCs, since that may represent a different physical quantity.
Here we take alpha to have Neumann BCs.  (Dirichlet BCs for alpha will
force it to 0 on the boundary, which is not correct here)

"""

from __future__ import print_function

import os

import numpy as np
import matplotlib.pyplot as plt

import compare
import mesh.boundary as bnd
import mesh.patch as patch
import multigrid.variable_coeff_MG as MG
from util import msg, io

# the analytic solution
def true(x,y):
    return np.sin(2.0*np.pi*x)*np.sin(2.0*np.pi*y)

# the coefficients
def alpha(x,y):
    return 2.0 + np.cos(2.0*np.pi*x)*np.cos(2.0*np.pi*y)

# the righthand side
def f(x,y):
    return -16.0*np.pi**2*(np.cos(2*np.pi*x)*np.cos(2*np.pi*y) + 1)* \
        np.sin(2*np.pi*x)*np.sin(2*np.pi*y)



def test_vc_poisson_periodic(N, store_bench=False, comp_bench=False,
                             make_plot=False, verbose=1):
    """
    test the variable-coefficient MG solver.  The return value
    here is the error compared to the exact solution, UNLESS
    comp_bench=True, in which case the return value is the
    error compared to the stored benchmark
    """

    # test the multigrid solver
    nx = N
    ny = nx


    # create the coefficient variable
    g = patch.Grid2d(nx, ny, ng=1)
    d = patch.CellCenterData2d(g)
    bc_c = bnd.BC(xlb="periodic", xrb="periodic",
                  ylb="periodic", yrb="periodic")
    d.register_var("c", bc_c)
    d.create()

    c = d.get_var("c")
    c[:,:] = alpha(g.x2d, g.y2d)

    # check whether the RHS sums to zero (necessary for periodic data)
    rhs = f(g.x2d, g.y2d)
    print("rhs sum: {}".format(np.sum(rhs[g.ilo:g.ihi+1,g.jlo:g.jhi+1])))


    # create the multigrid object
    a = MG.VarCoeffCCMG2d(nx, ny,
                          xl_BC_type="periodic", yl_BC_type="periodic",
                          xr_BC_type="periodic", yr_BC_type="periodic",
                          coeffs=c, coeffs_bc=bc_c,
                          verbose=verbose, vis=0, true_function=true)


    # initialize the solution to 0
    a.init_zeros()

    # initialize the RHS using the function f
    rhs = f(a.x2d, a.y2d)
    a.init_RHS(rhs)

    # solve to a relative tolerance of 1.e-11
    a.solve(rtol=1.e-11)

    # alternately, we can just use smoothing by uncommenting the following
    #a.smooth(a.nlevels-1,10000)

    # get the solution
    v = a.get_solution()

    # get the true solution
    b = true(a.x2d,a.y2d)

    # compute the error from the analytic solution -- note that with
    # periodic BCs all around, there is nothing to normalize the
    # solution.  We subtract off the average of phi from the MG
    # solution (we do the same for the true solution to put them on
    # the same footing)
    e = v - np.sum(v.v())/(nx*ny) - (b - np.sum(b[a.ilo:a.ihi+1,a.jlo:a.jhi+1])/(nx*ny))

    enorm = e.norm()
    print(" L2 error from true solution = %g\n rel. err from previous cycle = %g\n num. cycles = %d" % \
          (enorm, a.relative_error, a.num_cycles))


    # plot the solution
    if make_plot:
        plt.clf()

        plt.figure(figsize=(10.0,4.0), dpi=100, facecolor='w')

        plt.subplot(121)

        plt.imshow(np.transpose(v.v()),
                   interpolation="nearest", origin="lower",
                   extent=[a.xmin, a.xmax, a.ymin, a.ymax])

        plt.xlabel("x")
        plt.ylabel("y")
        plt.title("nx = {}".format(nx))

        plt.colorbar()


        plt.subplot(122)

        plt.imshow(np.transpose(e.v()), 
                   interpolation="nearest", origin="lower",
                   extent=[a.xmin, a.xmax, a.ymin, a.ymax])

        plt.xlabel("x")
        plt.ylabel("y")
        plt.title("error")

        plt.colorbar()

        plt.tight_layout()

        plt.savefig("mg_vc_periodic_test.png")


    # store the output for later comparison
    bench = "mg_vc_poisson_periodic"
    bench_dir = os.environ["PYRO_HOME"] + "/multigrid/tests/"

    my_data = a.get_solution_object()

    if store_bench:
        my_data.write("{}/{}".format(bench_dir, bench))

    # do we do a comparison?
    if comp_bench:
        compare_file = "{}/{}".format(bench_dir, bench)
        msg.warning("comparing to {}".format(compare_file))
        bench = io.read(compare_file)

        result = compare.compare(my_data, bench)

        if result == 0:
            msg.success("results match benchmark\n")
        else:
            msg.warning("ERROR: {}\n".format(compare.errors[result]))

        return result


    # normal return -- error wrt true solution
    return enorm


if __name__ == "__main__":

    N = [16, 32, 64, 128, 256, 512]
    err = []

    plot = False
    store = False
    do_compare = False

    for nx in N:
        if nx == max(N):
            plot = True
            #store = True
            #do_compare = True
        enorm = test_vc_poisson_periodic(nx, make_plot=plot,
                                         store_bench=store, comp_bench=do_compare)
        err.append(enorm)


    # plot the convergence
    N = np.array(N, dtype=np.float64)
    err = np.array(err)

    plt.clf()
    plt.loglog(N, err, "x", color="r")
    plt.loglog(N, err[0]*(N[0]/N)**2, "--", color="k")

    plt.xlabel("N")
    plt.ylabel("error")

    f = plt.gcf()
    f.set_size_inches(7.0,6.0)

    plt.tight_layout()

    plt.savefig("mg_vc_periodic_converge.png")
