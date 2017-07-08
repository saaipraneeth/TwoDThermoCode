'''
 generic tools used for the thermodynamic computations
'''
from xml.etree import ElementTree
import numpy as np
import pdb
from mpi4py import MPI
from pdb import set_trace as keyboard

def getMolarFraction(MW_mix,MW,scalar):
    N_scalars,Nxc = np.shape(scalar)
    X = np.zeros(np.shape(scalar))
    if (N_scalars==1):
      return np.ones(N_scalars)
    else:
      for ii in range(self.NSPECIES):
        X[:,ii] = MW_mix/MW[ii]*scalar[:][ii]
      return X

def convertMassToMolar(phi,MW):
    #given the MW molar weight [g/mol], compute the mass 
    return phi/(MW*1E3)

def convertMolarToMass(phi,MW):
    #given mol/m^3 to
    return phi*MW*1E3

def getMolecularWeightFluid(MW,scalar):
    '''given the partial densities, compute the MW of single fluid
      @param MW molar weight [g/mol] '''
    N_scalars,Nxc = np.shape(scalar)
    MW  = np.zeros(Nxc)
    MW  = np.ones(Nxc)*MW
    return MW

def getMolecularWeightMixture(MW,scalar):
    '''given the partial densities, compute the MW of mixture'''
    N_scalars,Nxc = np.shape(scalar)
    MW_mix = np.zeros(Nxc)
    if (N_scalars==1):
      MW_mix = np.ones(Nxc)*MW
    else:
      for ii in range(N_scalars):
        MW_mix += scalar[:][ii]*MW[ii]
    return MW_mix


def getRhofromV(v,MW):
    '''computes density from molar volume
    @param MW molar weight [kg/mol]
     @param v molar volume [g/mol]
     @param rho density [kg/m^3] '''
    return (MW/v)/1E3


def getVfromRho(rho,MW):
    '''computes  molar volume from density'''
    return MW/(rho*1E3)


#------------------------------------------------------------------
# Finds the Cp and Cv of the ideal gas. Output in J/(kg K)
#------------------------------------------------------------------
def getCp_ideal(coef,T,RoM):  
    return RoM*(coef[:,0]  + coef[:,1]*T + coef[:,2]*T**2  + coef[:,3]*T**3  + coef[:,4]*T**4 )

def getCv_ideal(coef,T,RoM):
    Cp0 = getCp_ideal(coef,T,RoM)
    Cv0 = Cp0 - RoM
    return(Cv0)

def getH_ideal(coef,T,RoM):

#    mpi_rank = MPI.COMM_WORLD.Get_rank()
#    print "\n ###### ############## ############### \n"
#    print "\n I am mpi_rank = " + repr(mpi_rank)
#    print "\n ######  coef.shape" + repr(coef.shape)

    #if (True in np.isnan(T)):
	#os._exit("NAN found in temperature in processor " + repr(mpi_rank))

    H_ideal = T*RoM*(coef[:,0] + coef[:,1]*T/2.0 + coef[:,2]*T**2/3.0 + coef[:,3]*T**3/4.0 + coef[:,4] * T**4/5.0+ coef[:,5] / T)
    return H_ideal

def getE_ideal(coef,T,RoM):
    return getH_ideal(coef,T,RoM) - RoM*T



def cubicSolver(a,b,c,d):	
	# We solve for the equation of the form: 
	# ax^3 + bx^2 + cx + d =0
	# Followed the steps in http://www.proofwiki.org/wiki/Cardan%27s_Formula
  Q = (3.0*a*c-b**2)/(9.0*a**2)
  R = (9.0*a*b*c - 27.0*a**2*d - 2.0*b**3)/(54.0*a**3.0)
  roots=np.zeros([3])
  D=Q**3.0+R**2.0	
  if D<0.0:
    #print "three real roots"		
    theta=np.arccos(R/np.sqrt(-Q**3))
    roots[0] = 2.0*np.sqrt(-Q)*np.cos(theta/3)-b/(3*a)
    roots[1] = 2.0*np.sqrt(-Q)*np.cos(theta/3 + 2.0*np.pi/3.0)-b/(3*a)
    roots[2] = 2.0*np.sqrt(-Q)*np.cos(theta/3 + 4.0*np.pi/3.0)-b/(3*a)	
  elif D>0.0:
    #print "one real and two complex conjugates", a,b,c,d
    St = R + (Q**3+R**2)**0.5
    Tt = R - (Q**3+R**2)**0.5
    S= abs(St)**(1.0/3.0)
    T= abs(Tt)**(1.0/3.0)
    if St<0.0: S=-S
    if Tt<0.0: T=-T
		
    roots[0]=   S + T    - b/(3.0*a)
	#roots[1]= -(S + T)/2.0 - b/(3.0*a) +1j*np.sqrt(3)/2.0*(S-T)
	#roots[2]= -(S + T)/2.0 - b/(3.0*a) -1j*np.sqrt(3)/2.0*(S-T)
	
  else:
    #print "all real but two are equal"
    if R>=0.0:	
      R3 = math.pow(R, 1.0/3.0)
    elif R<0.0:
      R3 = - math.pow(abs(R), 1.0/3.0)
      roots[0] = 2.0*R3 - b/(3.0*a)
      roots[1] =  -R3 - b/(3.0*a)
      roots[2] =   x2
   #print "the Q and R are:", Q,R
  return(roots)



def viscosity_Sutherland(T_in,muref,Tref,Sref):
  return muref*(T_in/Tref)**(3./2.)*((Tref+Sref)/(T_in+Sref))



def readCritPropDatabase(path_CritProp,specie):
  ''' Reads the critical property database'''
  species =[]
  Pcrit=[] 
  Tcrit=[]
  rhocrit=[]
  omega =[]
  Gasphase=[]
  dipole=[]
  MW=[]
  f = open(path_CritProp, 'r')
  lines = f.readlines()
  for line in lines[6:46]:     #I am lazy...I should add the correct MW
    splitline=line.split()
    species.append(splitline[0])
    Pcrit.append(float(splitline[2]))
    Tcrit.append(float(splitline[3]))
    rhocrit.append(float(splitline[4]))
    omega.append(float(splitline[5]))
    dipole.append(float(splitline[6]))
    MW.append(float(splitline[7]))  
  index=species.index(specie)
  return MW[index],Pcrit[index],Tcrit[index], rhocrit[index],omega[index] 


def parse_database_xml(path_NASA,myspecie):
  ''' Reads the NASA polynomial database and outputs the temperature bounds and the coefficiants'''


  if myspecie.lower() == 'c12h26':
    # Manual read for dodecane	>>> Not available in xml database
    tempBounds = np.array([200.,1000.,1000.,3500.])
    lowT = np.array([3.70187925E+01, 5.54721488E-02, -1.92079548E-05, 3.08175574E-09, -1.84800617E-13, -5.26984458E+04, -1.61453501E+02])
    highT = np.array([2.13264480E+01, -3.86394002E-02, 3.99476113E-04, -5.06681097E-07, 2.00697878E-10, -4.22475053E+04, -4.85848300E+01])

  elif myspecie.lower() == 'r134a':
    # Manual read for R-134a (C2H2F4)	>>> Not available in xml database
    tempBounds = np.array([200.,1000.,1000.,6000.])
    lowT = np.array([2.29239681, 0.030310848, -5.33714E-6, -2.19457E-8, 1.2997E-11, -111790.431, 16.2830568])
    highT = np.array([12.5551115, 0.008401861, -3.12077E-6, 5.12285E-10, -3.1011E-14, -114846.319, -38.0374329])


  else:
	  #ElementTree to open the database
	  document = ElementTree.parse(path_NASA)

	  #Find all the species elements
	  species=document.findall(".//species")

	  #Loop to ID the desired specie
	  selected=-1
	  for specie in species:
	    if specie.attrib.get('name')==myspecie:
	      selected=specie


	  if selected==-1: 
	    quit("Species () was not found in the database")

	  # the database structure tells us that 'thermo' is the third
	  thermo=selected[2]

	  # Temperature bounds
	  tempBounds=np.zeros(4)
	  tempBounds[0]=thermo[0].attrib.get('Tmin')
	  tempBounds[1]=thermo[0].attrib.get('Tmax')
	  tempBounds[2]=thermo[1].attrib.get('Tmin')
	  tempBounds[3]=thermo[1].attrib.get('Tmax')

	  lowT=[float(x) for x in thermo[0][0].text.replace("\n", "").split(",")]
	  highT=[float(x) for x in thermo[1][0].text.replace("\n", "").split(",")]
	  lowT = np.array(lowT)
	  highT = np.array(highT)

  return   tempBounds,lowT,highT



