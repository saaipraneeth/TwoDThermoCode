'''
 class containing the perfect gas thermodynamics for single species
 For efficiency, we only compute all mixture properties in the consToPrim function
'''
import numpy as np
import thermodynamics_tools as tools
import pdb
from pdb import set_trace as keyboard
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI

class peng_robinson_fluid():
    '''
      Peng_Robinson Equation of State

                   RT              A
            P =  ______  _   _____________

                 v - B      v^2 + 2*B*v - B^2
      where, 
                   0.457236*R^2*T_c*2                 0.07780*R*T_c
            A = ( ___________________ )*alpha  , B = _______________ ,
                          p_c                             p_c

            alpha = (1 + k*(1 - T_r^0.5))**2 ,  k = 0.37464 + 1.54226*w - 0.26992*w**2  , T_r = T/T_c

            T_c, p_c, accentric_factor (w), k1, k2, k3 (PRSV modification..)

    '''
    def __init__(self):
      #self.setup    = setup
      self.Rcst     = 8.3144621	                #J/(mol*K)
      self.Nxc      = 18		# Loading number of grid points (with ghost cells pre-accounted for) from grid object
      self.gamma    = np.ones(self.Nxc)		# INCORRECT # what about CH3OH and C12H26 ?????

      self.mu       = np.zeros(self.Nxc)
      self.c_p      = np.zeros(self.Nxc)
      self.c_v      = np.zeros(self.Nxc) 
      self.G        = np.zeros(self.Nxc)
      self.MW       = np.zeros(self.Nxc)
      self.A        = np.zeros(self.Nxc)
      self.B        = np.zeros(self.Nxc)
      self.Told     = np.ones(self.Nxc)*100.0    # allows a speed up of the newton solver
      self.Pold     = np.ones(self.Nxc)*10000.0

      self.viscous   = True #setup.settings['viscous']
      self.Prandtl  = 0.7 #setup.settings['Prandtl']
      self.species  = 'O2' #setup.settings['species']      
      self.N_species= 1 #setup.settings['N_species']
      self.visc_Tref= 273.15 #setup.settings['viscosity_Sutherland_Tref']
      self.visc_muref= 1.716e-05 #setup.settings['viscosity_Sutherland_muref']
      self.visc_Sref= 110.4 #setup.settings['viscosity_Sutherland_Sref']


      self.MW         = 0.0
      self.RoM        = 0.0
      self.Pr         = 0.0
      self.Tcrit      = 0.0
      self.Pcrit      = 0.0
      self.rhocrit    = 0.0
      self.vcrit      = 0.0
      self.Zcrit      = 0.0
      self.omega      = 0.0
      self.coef_Cp    = 0.0

      self.Tbounds    = [None]       # temperature bounds for NASA polynomials
      self.NASAcoeff_lowT  = [None]       # NASA coefficients for low temperature
      self.NASAcoeff_highT = [None]  # NASA coefficients for high temperature

      # get critical properties of single species
      #critProp_db   = setup.settings['path_crit_prop_database']
      critProp_db   = '/home/administrator/pyro2/compressible/databases/thermodynamic_CritPropDatabase.dat'
      MW,Pcrit,Tcrit, rhocrit,omega = tools.readCritPropDatabase(critProp_db,self.species)
      self.MW       = MW
      self.RoM      = self.Rcst/self.MW
      self.Pcrit    = Pcrit
      self.Tcrit    = Tcrit
      self.rhocrit  = rhocrit
      self.vcrit    = self.MW/(self.rhocrit * 1000.0)
      self.Zcrit    = Pcrit*self.vcrit/(self.Rcst*Tcrit)
      self.omega    = omega

     # get NASA polynomials
      #NASA_db             = setup.settings['path_NASA_poly_database']
      NASA_db             = '/home/administrator/pyro2/compressible/databases/thermodynamic_NASAdatabase.xml'
      outArray            = tools.parse_database_xml(NASA_db,self.species)
      self.NASAbounds     = np.array(outArray[0])
      self.NASAcoeff      = [None]*(self.Nxc)
      self.NASAcoeff_lowT = np.array([np.array(outArray[1])]*self.Nxc)  # create an array of coefficients! (efficiency reasons)
      self.NASAcoeff_highT= np.array([np.array(outArray[2])]*self.Nxc)

      
     # compute here fro single species
      self.computeB()
      self.computeC()

    def getDensity(self,T, P):
      self.computeA(T)
      return self.getDensityfromPressureTemperature(P,T)


    def computeA(self, T_in):
      '''
       Computes the A coefficient for the Peng-Robinson equation (temperature dependant)
       @param T     Temperature [K] (array) 
       @param R     Universal gas constant [J/(mol*K)] 
       @param Pcrit Critical pressure [Pa]
       @param Tcrit Critical temperature [K]
       @param A     Peng Robinson constant [J/(mol*Pa) = m^3/mol]
      '''
      R     = self.Rcst 
      Tcrit = self.Tcrit
      Pcrit = self.Pcrit     
      omega = self.omega
      ones  = np.ones(self.Nxc)
      kappa = 0.37464 + 1.54226*omega-0.26992*omega**2
      #k_0   = 0.378893 + 1.4897153*omega - 0.17131848*omega**2 + 0.0196554*omega**3
      #kappa = k_0 + [k1 + k2(k3 - T_r)(1-T_r**0.5)](1 + T_r**0.5)(0.7 - T_r)
      alpha = (ones+kappa*(ones-(T_in/Tcrit)**0.5))**2
      self.A = (0.457236 * alpha * (R*Tcrit)**2/Pcrit)
      return self.A


    def computeB(self):
      ''' Computes the B coefficient for the Peng-Robinson equation
       @param B     Peng Robinson constant [m^3/mol]     '''
      R  = self.Rcst
      Tcrit = self.Tcrit
      Pcrit = self.Pcrit
      ones = np.ones(self.Nxc)
      self.B = 0.077796*R*Tcrit/Pcrit*ones
      return

    def computeC(self):
      ''' Computes the C coefficient for the Peng-Robinson equation
       @param C     Peng Robinson constant [-]     '''
      self.C = 0.37464 + 1.54226 * self.omega - 0.26992 *self.omega**2
      return

    def computeG(self,T_in):   #Miller 2001 (appendix)
      '''  Computes the G coefficient for the Peng-Robinson equation  '''
      self.G = self.C*((T_in/self.Tcrit)**0.5)/(1.0+self.C*(1.0-(T_in/self.Tcrit)**0.5))
      return self.G


    def computeK1(self,v_in):
      ''' Computes the K1 coefficient for the Peng-Robinson equation (Miller 2001 eq 2.13)  '''
      B = self.B
      ones = np.ones(self.Nxc)
      self.K1 = 1.0/(2.0*np.sqrt(2.0)*B) * np.log((v_in+(1.0-np.sqrt(2))*B)/(v_in+(1.0+np.sqrt(2))*B*ones))
      return self.K1


    def setNASAcoeff(self,T_in):
      ''' Sets the correct NASA coefficients depending on the temperature '''
      if isinstance(T_in, float):
        T_in = np.ones(self.Nxc)*T_in
      bounds = self.NASAbounds
      lowT   = np.where(T_in <= bounds[1])[0]
      highT  = np.where(T_in  > bounds[1])[0]
      self.NASAcoeff      = [None]*(self.Nxc)
      # Not optimal programming... (needs to be fixed)
      for index in lowT:    
        self.NASAcoeff[index] = self.NASAcoeff_lowT[index]
      for index in highT:    
        self.NASAcoeff[index] = self.NASAcoeff_highT[index]
      #convert to an array of (Nx,7)
      self.NASAcoeff = np.asarray(self.NASAcoeff)



    def computed2AdT2(self, T_in):  # According to Miller2001   
      R   = self.Rcst
      C   = self.C
      self.d2AdT2 =  (0.457236*R**2/(2.0*T_in))*C*(1.0+C)*self.Tcrit/self.Pcrit*(self.Tcrit/T_in)**0.5



    def getAveragedGasConstant(self,Rcst,MW):
      '''  Return the universal constant over molar mass [J/(g*K)]  '''
      return Rcst/MW

    def getdAdT(self,T_in,A,G):          # According to Miller2001
      return(-A*G/T_in)

    def getdPdT(self,v_in,T_in):			#Coussement 3.70
      R     = self.Rcst
      B     = self.B
      dAmdT = self.getdAdT(T_in,self.A,self.G) 	
      return R/(v_in-B) - dAmdT/(v_in**2 + 2*v_in*B - B**2)  #Miller 2001 2.11

    def getdPdv(self,v_in,T_in): 		#Coussement 3.71
      R   = self.Rcst
      A   = self.A
      B   = self.B
      coef = R*T_in/(v_in-B)**2
      para = (R*T_in*(v_in+B)*(v_in/(v_in-B)+ B/(v_in+B))**2)
      return - coef * (1.0 - 2.0*A/para)	#Miller 2001 

	

    def getCv(self,v_in,T_in):
      '''         Returns the molar value of cv [J/(mol K)]    '''
      Cv0     = tools.getCv_ideal(self.NASAcoeff,T_in,self.Rcst)
      d2AmdT2 = self.d2AdT2
      deltaCv = -self.K1 * T_in * d2AmdT2	
      return Cv0 + deltaCv



    def getCp(self,v_in,T_in):
      '''     Returns the molar value of cp [J/(mol K)]    '''
      Cp0     = tools.getCp_ideal(self.NASAcoeff,T_in,self.Rcst)
      d2AdT2 = self.d2AdT2
      dPdv    = self.getdPdv(v_in,T_in)
      dPdT    = self.getdPdT(v_in,T_in)
      if isinstance(T_in, float):
        ones    = np.ones(1)
      else:
        ones    = np.ones(len(T_in))
      deltaCp = - self.K1 * T_in * d2AdT2 - T_in*(dPdT)**2/dPdv - self.Rcst*ones
      return Cp0 + deltaCp

    
    def getPressurefromVolumeTemperature(self,v_in,T_in):
      '''   get pressure given molar volume, T  '''
      R = self.Rcst
      B = self.B
      A = self.A
      ones = np.ones(np.size(T_in))
      return (R*T_in/(v_in-B*ones))-A/((v_in**2 + 2*v_in*B - ones*B**2))



    def getDensityfromPressureTemperature(self,P_in,T_in):
      ''' Computes the density from the pressure and temperature.
         The values of A and B must be updated before!!!'''
      R     = self.Rcst
      A     = self.A
      B     = self.B
      Nxc   = np.size(T_in)
      ones  = np.ones(Nxc)
      rho   = np.ones(Nxc)
      coef_a = ones
      coef_b = B*ones - R * T_in/P_in
      coef_c = - 2.0* B * R * T_in/P_in  +  A/P_in - 3*B**2*ones
      coef_d = B**2*R*T_in/P_in + B**3*ones - A*B/P_in
      # for ii in range(Nxc):
      #   roots = tools.cubicSolver(coef_a[ii],coef_b[ii],coef_c[ii],coef_d[ii])
      #   v     = roots[0]
      #   rho[ii] = tools.getRhofromV(v,self.MW)   
      rho = PropsSI('DMASS', 'P', P_in, 'T', T_in, 'Oxygen')
      return rho

    def getTemperature(self,rho,p):
      return  self.NewtonIterate_TemperaturefromPrho(rho,p)
    

    def updateGamma(self, v_in, T_in):
	Cp = tools.convertMolarToMass(self.getCp(v_in, T_in), self.MW)
	Cv = tools.convertMolarToMass(self.getCv(v_in, T_in), self.MW)
	self.gamma = Cp/Cv

    def getCompressibilityFactor(self, p, rho, T):
	v_in = tools.getVfromRho(rho, self.MW)
	self.setRealFluidThermodynamics(v_in, T)
	self.updateGamma(v_in, T)
	return (p/(rho*self.RoM*1000*T))

    def getIsothermalComp(self, v_in, T_in):      #Miller 2001 2.16
      dPdv    = self.getdPdv(v_in,T_in)
      K_T = -1.0/(v_in*(dPdv))
      #print "Isothermal comp ", K_T[0]
      return K_T

    def getExpansivity(self, v_in, T_in):         #Miller2001 2.16
      dPdv    = self.getdPdv(v_in,T_in)
      dPdT    = self.getdPdT(v_in,T_in)
      alpha_v = -(dPdT / (v_in*dPdv)) 
      return alpha_v 

    def getIsentropicComp(self, v_in, T_in): #Miller 2001 2.15
      K_T = self.getIsothermalComp(v_in,T_in)
      alpha_v = self.getExpansivity(v_in,T_in)
      Cpt=self.getCp(v_in, T_in)
      Cp = tools.convertMolarToMass(Cpt, self.MW)/(self.MW*1000)
      K_S = K_T - (v_in*T_in*(alpha_v**2)/Cp)
      return K_S


    def getSpeedOfSound(self, P_in, T_in): #Miller 2001 2.14

      if isinstance(T_in, float):
        self.computed2AdT2(T_in)
        rho = self.getDensityfromPressureTemperature(P_in,T_in)[0]  
        v_in = tools.getVfromRho(rho, self.MW)
        self.setRealFluidThermodynamics(v_in, T_in)
        self.computeK1(v_in)
        K_S = self.getIsentropicComp(v_in, T_in)
        return np.sqrt(1.0/(rho*K_S))
      rho = self.getDensityfromPressureTemperature(P_in,T_in)[0]     
      v_in = tools.getVfromRho(rho, self.MW)
      K_S = self.getIsentropicComp(v_in, T_in)
      return np.sqrt(1.0/(rho*K_S))

    def getPrimitiveFromConservative(self,rho,rhou,rhoE,rhoScalar):
      '''       get primitives from conservatives variables     '''
      u           = rhou/rho
      scalar      = rhoScalar/rho
      Eint_mass   = rhoE/rho - 0.5*u*u     
      Eint_mol    = tools.convertMassToMolar(Eint_mass,self.MW)
  
      self.v      = tools.getVfromRho(rho,self.MW)

      T           = self.NewtonIterate_TemperaturefromEv(Eint_mol,self.v, self.Told)
      P           = self.getPressurefromVolumeTemperature(self.v,T)

      self.Told   = T
 
      return u,P,T,scalar
     
             
    def getConservativeFromPrimitive(self,u,P,T,scalar):
      '''    get conservatives from primitives  variables  '''
      self.computeA(T)
      rho     = self.getDensityfromPressureTemperature(P,T)
      v       = tools.getVfromRho(rho,self.MW)
      self.setRealFluidThermodynamics(v,T)
      Eint_mol = self.getEnergyfromVolumeTemperature(v,T)
      Eint_mass= tools.convertMolarToMass(Eint_mol,self.MW)		#

      rhou = rho*u
      rhoE = rho*(Eint_mass+0.5*u*u)			# Total energy
      rhoScalar=rho*scalar
      return rho,rhou,rhoE,rhoScalar


    def getEnthalpyFromTemperature(self,T_in,rho_in,p_in):      
      v_in       = tools.getVfromRho(rho_in,self.MW)
      return self.getEnthalpyfromVolumeTemperature(v_in,T_in,p_in)


    def getTransportProperties(self,T_in):
      '''    Currently constant!!! (naive)      '''
      if self.viscous:
        mu = tools.viscosity_Sutherland(T_in,self.visc_muref,self.visc_Tref,self.visc_Sref)
        Cp_ideal = tools.getCp_ideal(self.NASAcoeff,T_in,self.Rcst/self.MW*1000)
      else:
        mu = np.zeros(self.Nxc)
        Cp_ideal = np.zeros(self.Nxc)
      return mu, Cp_ideal*mu/self.Prandtl

    def NewtonIterate_TemperaturefromEv(self,e_target,v_target, T_in,eps=1E-6,omega=1.0):
      #    Newton Iteration to find temperature. Needs to pass energy in molar form!!!  
      T_n = T_in * np.ones(v_target.size)
      self.setRealFluidThermodynamics(v_target,T_n)
      e_n = self.getEnergyfromVolumeTemperature(v_target,T_n)
      cv_n = self.getCv(v_target,T_n)
      diff = (e_n - e_target)
      prevdiff = diff.copy()
      itera=0 
      #Newton solver
      range_convergence_complete = False
      pointwise_converged=np.zeros(diff.shape)
      flagged = False
      # while not(range_convergence_complete):    
      #   diff0 = max(abs(diff))            # old max diff
      #   prevdiff = diff.copy()            # old diff array
      #   T_correction = omega * diff * e_target/cv_n;      # corrections
      #   T_correction[np.where(pointwise_converged==True)] = 0.    # filtering corrections
      #   T_n -= T_correction           # MODIFY
      #   keyboard()
      #   self.setRealFluidThermodynamics(v_target,T_n)     # set new state
      #   e_n = self.getEnergyfromVolumeTemperature(v_target,T_n)   # evaluate new energy
      #   cv_n = self.getCv(v_target,T_n)         # new Cv
      #   diff = (e_n - e_target)/e_target            # new difference
      #   diff1 = max(abs(diff))            # new max diff

      #   if (diff1 > diff0):
      #     if (pointwise_converged[0]) and (pointwise_converged[-1]):  # if ghost cells have converged (for continuity)
      #       break;
      #     elif itera > 10:    # Taking too long AND ends not converged yet
      #       print "Taking too long"
      #       flagged = True
      #       break;
      #     # elif np.mean(prevdiff) < np.mean(diff):     # average diverging 
      #     #   flagged = True
      #    #    print "Average diverging failed"
      #     #   break;

      #   pointwise_converged = (abs(diff)<eps)
      #   range_convergence_complete = not(False in pointwise_converged)
      #   pointwise_converged = np.array(pointwise_converged)
      #   itera+=1

      while (max(abs(diff) > 1.0E-4)):
        dpdt = self.getdPdT(v_target, T_n)
        T_n = T_n - (diff/dpdt)
        print T_n
        self.setRealFluidThermodynamics(v_target, T_n)
        e_n = self.getEnergyfromVolumeTemperature(v_target, T_n)
        diff = e_n - e_target
      keyboard()
      return T_n

    def NewtonIterate_TemperaturefromPrho(self,rho_target,P_target, T_in=300.0,eps=1E-10,omega=1.0):
      T_n = T_in * np.ones(rho_target.shape)  
      v_target= tools.getVfromRho(rho_target,self.MW)
      self.setRealFluidThermodynamics(v_target,T_n)
      P_n  = self.getPressurefromVolumeTemperature(v_target,T_n)
      # if P_target.any() == 0:
      #   keyboard()
      diff = (P_n - P_target) / P_target
      prevdiff = diff
      diff0 = max(abs(diff))
      #print "diff , diff0", diff, diff0
      itera=0
      #Newton solver
      range_convergence_complete = False
      pointwise_converged = np.zeros(diff.shape)
      flagged = False
      while not range_convergence_complete:
        diff0 = abs(max(diff))
        prevdiff = diff
        dPdT = self.getdPdT(v_target,T_n)
        T_correction = omega * diff * P_target/dPdT
        T_correction[np.where(pointwise_converged == True)] = 0.
        T_n -= T_correction
        self.setRealFluidThermodynamics(v_target,T_n)        
        P_n  = self.getPressurefromVolumeTemperature(v_target,T_n)
        diff = (P_n - P_target) / P_target
        diff1 = max(abs(diff))

        if (diff1 > diff0):
          if (pointwise_converged[0]) and (pointwise_converged[-1]):  # if ghost cells have converged (for continuity)
            break;
          elif itera > 10:    # Taking too long AND ends not converged yet
            print "Taking too long"
            flagged = True
            break;
          elif np.mean(prevdiff) < np.mean(diff):     # average diverging 
            print "The average diverging is not optimistic"
            flagged = True
            break;

        pointwise_converged = abs(diff) < eps  #Still converging means FALSE, if converged then True
        range_convergence_complete = not(False in pointwise_converged) #Still converging means False, if converged True
        pointwise_converged = np.array(pointwise_converged)
        itera+=1

      return T_n

    def setRealFluidThermodynamics(self,v_in, T_in):
      self.computeA(T_in)
      self.computeG(T_in)
      self.computeK1(v_in)
      self.computed2AdT2(T_in)
      self.setNASAcoeff(T_in)
      
      
  
    def getEnergyfromVolumeTemperature(self,v_in,T_in):
      self.setRealFluidThermodynamics(v_in,T_in)
      e0   = tools.getE_ideal(self.NASAcoeff,T_in,self.Rcst)
      dAdT = self.getdAdT(T_in,self.A,self.G)

      int_e = e0 + self.K1* ( self.A - T_in*dAdT )
      return e0 + self.K1* ( self.A - T_in*dAdT )

    def getEnthalpyfromVolumeTemperature(self,v_in,T_in,p_in):    
      return self.getEnergyfromVolumeTemperature(v_in,T_in) + p_in*v_in

    def getCriticalInformation(self):
	critdict = {}
	critdict["EOS"] = self.__class__.__name__		# peng robinson
	critdict["species"] = self.species
	critdict["MW"]     = self.species
	critdict["Tcrit"] = self.Tcrit
	critdict["Pcrit"] = self.Pcrit
	critdict["rhocrit"] = self.rhocrit
	critdict["vcrit"] = self.vcrit
	return critdict
    def getPressure(self,rho,T_in):
      '''
       get pressure given rho, T
      '''
      v_in = tools.getVfromRho(rho, self.MW)
      return self.getPressurefromVolumeTemperature(v_in,T_in)
