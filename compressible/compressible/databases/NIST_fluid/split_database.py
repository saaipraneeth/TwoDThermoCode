import numpy as np
import pickle
from pdb import set_trace as keyboard

filename = 'NIST_tables_2000.pkl'

print "Loading pickle file... " + filename
reader = open(filename, 'rb')
datadict = pickle.load(reader)
reader.close()

species = datadict['species']

ext = '.pkl'

for specie in species:
	minidict = datadict[specie]
	
	minidict['info']['species_name'] = specie	# Saving identifying bit of information to info dict within species dict

	writefile = specie + ext
	print "\n Writing to pickle file... " + writefile
	writer = open(writefile, 'wb')
	pickle.dump(minidict, writer)
	writer.close()
	#keyboard()

print "\n\n#### End of Program Reached ####"

