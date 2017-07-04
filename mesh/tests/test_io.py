import mesh.boundary as bnd
import mesh.patch as patch
import mesh.array_indexer as ai
import numpy as np
import util.io as io

from numpy.testing import assert_array_equal

def test_write_read():

    myg = patch.Grid2d(8, 6, ng = 2, xmax=1.0, ymax=1.0)
    myd = patch.CellCenterData2d(myg)

    bco = bnd.BC(xlb="outflow", xrb="outflow",
                 ylb="outflow", yrb="outflow")
    myd.register_var("a", bco)
    
    myd.create()

    a = myd.get_var("a")
    a.v()[:,:] = np.arange(48).reshape(8, 6)

    myd.write("io_test")


    # now read it in
    nd = io.read("io_test")
    ng = nd.grid

    anew = nd.get_var("a")

    assert_array_equal(anew.v(), a.v())

