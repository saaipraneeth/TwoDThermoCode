#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import argparse
from util import io

# plot an output file using the solver's dovis script

def makeplot(plotfile, variable, outfile):

    sim = io.read(plotfile)
    myd = sim.cc_data
    myg = myd.grid

    if variable == "vort":
        vx = myd.get_var("x-velocity")
        vy = myd.get_var("y-velocity")

        v = myg.scratch_array()

        v[myg.ilo:myg.ihi+1,myg.jlo:myg.jhi+1] = \
             0.5*(vy[myg.ilo+1:myg.ihi+2,myg.jlo:myg.jhi+1] -
                  vy[myg.ilo-1:myg.ihi,myg.jlo:myg.jhi+1])/myg.dx - \
             0.5*(vx[myg.ilo:myg.ihi+1,myg.jlo+1:myg.jhi+2] -
                  vx[myg.ilo:myg.ihi+1,myg.jlo-1:myg.jhi])/myg.dy

    else:
        v = myd.get_var(variable)



    plt.figure(num=1, figsize=(6.0,6.0), dpi=100, facecolor='w')

    plt.imshow(np.transpose(v[myg.ilo:myg.ihi+1,myg.jlo:myg.jhi+1]),
               interpolation="nearest", origin="lower",
               extent=[myg.xmin, myg.xmax, myg.ymin, myg.ymax])

    plt.axis("off")
    plt.subplots_adjust(bottom=0.0, top=1.0, left=0.0, right=1.0)
    plt.savefig(outfile)


def get_args():

    parser = argparse.ArgumentParser()

    parser.add_argument("-o", type=str, default="plot.png",
                        metavar="plot.png", help="output file name")
    parser.add_argument("plotfile", type=str, nargs=1,
                        help="the plotfile you wish to plot")
    parser.add_argument("variable", type=str, nargs=1,
                        help="the name of the solver used to run the simulation")

    args = parser.parse_args()

    return args


if __name__== "__main__":

    args = get_args()

    makeplot(args.plotfile[0], args.variable[0], args.o)




