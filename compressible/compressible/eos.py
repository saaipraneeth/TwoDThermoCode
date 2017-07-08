"""
This is a gamma-law equation of state: p = rho e (gamma - 1), where
gamma is the constant ratio of specific heats.
"""
from CoolProp.CoolProp import PropsSI
from CoolProp.CoolProp import PhaseSI
fluid = 'Oxygen'
from pdb import set_trace as keyboard
import numpy as np
#from preos import peng_robinson_fluid as PREOS
import thermodynamics_tools as tools

def pres(dens, eint):
    """
    Given the density and the specific internal energy, return the
    pressure

    Parameters
    ----------
    gamma : float
        The ratio of specific heats
    dens : float
        The density
    eint : float
        The specific internal energy

    Returns
    -------
    out : float
       The pressure

     """
    #p = dens*eint*(gamma - 1.0)
    # p = eint.copy()
    # if (eint.ndim == 2):
    #     for i in range(np.shape(eint)[0]):
    #         for j in range(np.shape(eint)[1]):
    #             if dens[i][j] < 0.1 :
    #                 p[i][j] = 0.0
    #                 continue
    #             p[i][j] = PropsSI('P', 'UMASS', eint[i][j],'DMASS', dens[i][j], fluid)
    #     p = np.array(p, order = 'F')
    #     return p
    # else:
    #     p = PropsSI('P', 'UMASS', eint,'DMASS', dens, fluid)
    #     p = np.array(p, order = 'F')
    #     return p
    MW = 32
    vol = tools.getvfromrho(dens, MW)
    T_in = PREOS.NewtonIterate_TemperaturefromEv(eint,vol, T_in = 300.0 ,eps=1E-6,omega=1.0)
    p = PREOS.getPressurefromVolumeTemperature(vol,T_in)
    return p


def dens(pres, eint):
    """
    Given the pressure and the specific internal energy, return
    the density

    Parameters
    ----------
    gamma : float
        The ratio of specific heats
    pres : float
        The pressure
    eint : float
        The specific internal energy

    Returns
    -------
    out : float
       The density

    """
    #dens = pres/(eint*(gamma - 1.0))
    keyboard()
    dens = PropsSI('DMASS', 'UMASS', eint,'P', pres, fluid)
    dens = np.array(dens, order = 'F')
    return dens


def rhoe(dens, pres):
    """
    Given the pressure, return (rho * e)

    Parameters
    ----------
    gamma : float
        The ratio of specific heats
    pres : float
        The pressure

    Returns
    -------
    out : float
       The internal energy density, rho e

    """
    #rhoe = pres/(gamma - 1.0)
    #eint = pres.copy()
    # if (eint.ndim == 2):
    #     for i in range(np.shape(pres)[0]):
    #         for j in range(np.shape(pres)[1]):
    #             if dens[i][j] < 0.1 :
    #                 eint[i][j] = 0.0
    #                 continue
    #             eint[i][j] = PropsSI('UMASS', 'P', pres[i][j],'DMASS', dens[i][j], fluid)
    #     rhoe = np.array(dens*eint, order = 'F')
    #     return rhoe
    # else:
    #     eint = PropsSI('P', 'UMASS', eint,'DMASS', dens, fluid)
    #     rhoe = np.array(dens*eint, order = 'F')
    #     return rhoe

    T_in = PREOS.NewtonIterate_TemperaturefromPrho(dens, pres, T_in = 300.0, eps= 1E-10, omega = 1.0)
    eint = PREOS.getEnergyfromVolumeTemperature(T_in)
    rhoe = dens*eint
    #rhoe = np.array(dens*eint, order = 'F')
    return rhoe

    # eint = PropsSI('UMASS', 'P', pres, 'DMASS', dens, fluid)
    # rhoe = np.array(dens*eint, order = 'F')
    #return rhoe
