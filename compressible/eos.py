"""
This is a gamma-law equation of state: p = rho e (gamma - 1), where
gamma is the constant ratio of specific heats.
"""
from CoolProp.CoolProp import PropsSI
from CoolProp.CoolProp import PhaseSI
fluid = 'Nitrogen'
from pdb import set_trace as keyboard
import numpy as np
import preos_cy as PREOS
#import preos
import thermodynamics_tools as tools

#PREOS = preos.peng_robinson_fluid()

eqofst = 'coolprop'

if eqofst == 'PREOS':

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


        T = np.ones(np.shape(dens))
        p = np.ones(np.shape(dens))
        if (dens.ndim == 2):
            for i in range(np.shape(dens)[0]):
                # if dens[i].any() < 0.1:
                #     temp = dens[i]
                #     temp[temp < 0.1] = np.nan
                #     dens[i] = temp
                #     continue

                eint[i] = tools.convertMassToMolar(eint[i], 28.0134)
                T[i] = PREOS.getTfromEandRho(eint[i],dens[i])
                p[i] = PREOS.getPfromTandRho(T[i], dens[i])
            keyboard()
            return p
        else:
            eint= tools.convertMassToMolar(eint, 28.0134)
            T_in = PREOS.getTfromEandRho(eint,dens)
            p = PREOS.getPfromTandRho(T, dens)
            keyboard()
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
        MW = 28.0134
        T = np.ones(np.shape(dens))
        eint = np.ones(np.shape(dens))
        if (dens.ndim == 2):
            for i in range(np.shape(dens)[0]):
                # if dens[i].any() < 0.1:
                #     temp = dens[i]
                #     temp[temp < 0.1] = np.nan
                #     dens[i] = temp
                #     continue
                print pres[i], dens[i]
                T[i] = PREOS.getTfromPandRho(pres[i], dens[i])
                #T[i] = PropsSI('T', 'P', pres[i], 'DMASS', dens[i], 'Nitrogen')
                eint[i] = PREOS.getEnergyfromTandRho(T[i], dens[i])
                eint[i] = tools.convertMolarToMass(eint[i], MW) #J/Kg
            keyboard()
            return dens*eint
        else:
            T = PREOS.getTfromPandRho(pres, dens)
            #T = PropsSI('T', 'P', pres, 'DMASS', dens, 'Nitrogen')
            eint = PREOS.getEnergyfromTandRho(T, dens)
            eint = tools.convertMolarToMass(eint, MW)
            keyboard()
        return dens*eint


        # T_in = PREOS.NewtonIterate_TemperaturefromPrho(dens, pres, T_in = 300.0, eps= 1E-10, omega = 1.0)
        # eint = PREOS.getEnergyfromVolumeTemperature(T_in)
        # rhoe = dens*eint
        #rhoe = np.array(dens*eint, order = 'F')
        # return rhoe

        # eint = PropsSI('UMASS', 'P', pres, 'DMASS', dens, fluid)
        # rhoe = np.array(dens*eint, order = 'F')
        #return rhoe

elif eqofst == 'coolprop':

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
        # p = dens*eint*(gamma - 1.0)
        p = np.zeros(np.shape(eint))
        if (eint.ndim == 2):
            for i in range(np.shape(eint)[0]):
                for j in range(np.shape(eint)[1]):
                    if dens[i][j] < 0.1 :
                        p[i][j] = 0.0
                        continue
                    p[i][j] = PropsSI('P', 'UMASS', eint[i][j],'DMASS', dens[i][j], fluid)
            return p
        else:
            p = PropsSI('P', 'UMASS', eint,'DMASS', dens, fluid)
            keyboard()
            #p = np.array(p, order = 'F')
            return p


        # T = np.ones(np.shape(dens))
        # p = np.ones(np.shape(dens))
        # if (dens.ndim == 2):
        #     for i in range(np.shape(dens)[0]):
        #         # if dens[i].any() < 0.1:
        #         #     temp = dens[i]
        #         #     temp[temp < 0.1] = np.nan
        #         #     dens[i] = temp
        #         #     continue

        #         eint[i] = tools.convertMassToMolar(eint[i], 28.0134)
        #         T[i] = PREOS.getTfromEandRho(eint[i],dens[i])
        #         p[i] = PREOS.getPfromTandRho(T[i], dens[i])
        #     keyboard()
        #     return p
        # else:
        #     eint= tools.convertMassToMolar(eint, 28.0134)
        #     T_in = PREOS.getTfromEandRho(eint,dens)
        #     p = PREOS.getPfromTandRho(T, dens)
        #     keyboard()
        #return p

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
        eint = np.zeros(np.shape(pres))
        if (eint.ndim == 2):
            for i in range(np.shape(pres)[0]):
                for j in range(np.shape(pres)[1]):
                    if dens[i][j] < 0.1 :
                        eint[i][j] = 0.0
                        continue
                    eint[i][j] = PropsSI('UMASS', 'P', pres[i][j],'DMASS', dens[i][j], fluid)
            rhoe = np.array(dens*eint, order = 'F')
            return rhoe
        else:
            eint = PropsSI('UMASS', 'P', pres,'DMASS', dens, fluid)
            rhoe = np.array(dens*eint, order = 'F')
            return rhoe
        # MW = 28.0134
        # T = np.ones(np.shape(dens))
        # eint = np.ones(np.shape(dens))
        # if (dens.ndim == 2):
        #     for i in range(np.shape(dens)[0]):
        #         # if dens[i].any() < 0.1:
        #         #     temp = dens[i]
        #         #     temp[temp < 0.1] = np.nan
        #         #     dens[i] = temp
        #         #     continue
        #         print pres[i], dens[i]
        #         T[i] = PREOS.getTfromPandRho(pres[i], dens[i])
        #         #T[i] = PropsSI('T', 'P', pres[i], 'DMASS', dens[i], 'Nitrogen')
        #         eint[i] = PREOS.getEnergyfromTandRho(T[i], dens[i])
        #         eint[i] = tools.convertMolarToMass(eint[i], MW) #J/Kg
        #     keyboard()
        #     return dens*eint
        # else:
        #     T = PREOS.getTfromPandRho(pres, dens)
        #     #T = PropsSI('T', 'P', pres, 'DMASS', dens, 'Nitrogen')
        #     eint = PREOS.getEnergyfromTandRho(T, dens)
        #     eint = tools.convertMolarToMass(eint, MW)
        #     keyboard()
        #return dens*eint


        # T_in = PREOS.NewtonIterate_TemperaturefromPrho(dens, pres, T_in = 300.0, eps= 1E-10, omega = 1.0)
        # eint = PREOS.getEnergyfromVolumeTemperature(T_in)
        # rhoe = dens*eint
        #rhoe = np.array(dens*eint, order = 'F')
        # return rhoe

        # eint = PropsSI('UMASS', 'P', pres, 'DMASS', dens, fluid)
        # rhoe = np.array(dens*eint, order = 'F')
        #return rhoe