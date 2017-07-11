"""
This is a gamma-law equation of state: p = rho e (gamma - 1), where
gamma is the constant ratio of specific heats.
"""
from CoolProp.CoolProp import PropsSI
from CoolProp.CoolProp import PhaseSI
fluid = 'Oxygen'
from pdb import set_trace as keyboard
import numpy as np
import preos
import thermodynamics_tools as tools
PREOS = preos.peng_robinson_fluid()
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
    MW = 32.0
    vol = tools.getVfromRho(dens, MW)
    T = np.ones(np.shape(vol))
    p = np.ones(np.shape(vol))
    if (vol.ndim == 2):
        for i in range(np.shape(vol)[0]):
            if dens[i].any() < 0.1:
                temp = dens[i]
                temp[temp < 0.1] = np.nan
                dens[i] = temp
                continue
            eint[i] = tools.convertMassToMolar(eint[i], MW)
            eint[i] = tools.convertMolarToMass(eint[i], MW)
            T[i] = PREOS.NewtonIterate_TemperaturefromEv(eint[i],vol[i], T_in = 300.0 ,eps=1E-6,omega=1.0)
            p[i] = PREOS.getPressurefromVolumeTemperature(vol[i],T[i])
        return p
    else:
        eint = tools.convertMassToMolar(eint, MW)
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
    MW = 32.0
    T = np.ones(np.shape(dens))
    eint = np.ones(np.shape(dens))
    if (dens.ndim == 2):
        for i in range(np.shape(dens)[0]):
            if dens[i].any() < 0.1:
                temp = dens[i]
                temp[temp < 0.1] = np.nan
                dens[i] = temp
                continue
            T[i] = PREOS.NewtonIterate_TemperaturefromPrho(dens[i], pres[i], T_in = 300.0, eps= 1E-10, omega = 1.0)
            eint[i] = PREOS.getEnergyfromVolumeTemperature(dens[i], T[i])
        return dens*eint
    elif (np.shape(dens)[0] == 1224 ):
        T = PropsSI('T', 'DMASS', dens, 'P', pres, 'Oxygen')
        dens = dens.reshape(68, 18)
        T = T.reshape(68,18)
        eint = np.ones(np.shape(dens))
        for i in range(68):
            eint[i] = PREOS.getEnergyfromVolumeTemperature(dens[i],T[i])
        
        dens = np.ndarray.flatten(dens)
        eint = np.ndarray.flatten(eint)
        eint = tools.convertMolarToMass(eint, MW) #The energy calculated using NASA coeff will be in J/Mol
        return dens*eint
    else:
        T = PREOS.NewtonIterate_TemperaturefromPrho(dens, pres, T_in = 300.0, eps= 1E-10, omega = 1.0)
        eint = PREOS.getEnergyfromVolumeTemperature(dens,T)
        eint = tools.convertMolarToMass(eint, MW)
        return dens*eint