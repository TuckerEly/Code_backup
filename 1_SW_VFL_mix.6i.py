########################################################################################################################################
################################    1_SW_VFL_mix.6i.py    To prepare 6i VFL-SW mixes.    ###############################################
####################################      Tucker Ely, 27 DEC 2016   -   By Dillons Bedside    ##########################################
########################################################################################################################################	

import os, re, sys, shutil
import numpy as np
import pandas as pd

from itertools import *
from collections import *
import fileinput


import networkx as nx
import matplotlib.pyplot as plt
from compiler.ast import flatten
from scipy.integrate import ode
from time import gmtime, strftime


from PIL import *


"""	EDIT BELOW	"""

directory = '/Users/tuckerely/Google Drive/Ship/C-DEBI-DCO/EQ36_Test/2_6o_to_vfl.3i'
template = '/Users/tuckerely/Google Drive/Ship/C-DEBI-DCO/EQ36_Test/3_mix/mix_template.txt'

SW_3p = '/Users/tuckerely/Google Drive/Ship/C-DEBI-DCO/EQ36_Test/3_mix/tdf45.3p'

"""	EDIT ABOVE	"""

#################################################################
########################    Main    #############################
def main():

	"""Read in all 3p files"""
	file_name, file_list = read_3p()
	### Eventual index of files for which Zi did not reach maximum (failed)
	fail_list = [] 														#	for failed files

	### SW name

	print SW_3p[([pos for pos, char in enumerate(SW_3p) if char == '/'][-1]+1):-3]
	SW_name = SW_3p[([pos for pos, char in enumerate(SW_3p) if char == '/'][-1]+1):-3]
	
	print file_name

	
	"""Process all 6o files"""
	### Cycle through 6o files
	w = 0 																#	.6o file index
	while w < len(file_name): 											#	cycle through all 6o files, producing one 3i for each if zi is reached in 6o
		### Load 6o file as o
		o = open(file_name[w], "r") 									#	Read current file
		o_lines = o.readlines() 										#	stored 6o file lines
		o.close()

		print file_name[w]
		print 'Reaction of  ' + file_name[w][:-3] + '  Vent fluid back into seawater  ' + SW_name
		### open template
		### change rxn notes
		### load univeral tdf SW as special reactant
		### load risge specific 3p VFL into bottom of 6i
		







		w += 1  									#	onto the next 6o file

#################################################################
#####################    Functions     ##########################

def read_3p():
	### this function can find all .6o files in all downstream folders, hence the additio of file_list
	file_name = [] 						#	file names
	file_list = []						# 	file names with paths
	for root, dirs, files in os.walk(directory):# this works. 
		for file in files:
			if file.endswith(".3p"):
				file_name.append(file)
				file_list.append(os.path.join(root, file)) 					
	
	return file_name, file_list





def write_3i(w, file_name, new_basis_value, t, fO2, pH):

	new_file = file_name[w][:-3] + '_mix.6i' 					#	new file title
	shutil.copyfile(template, new_file) 						#	copy 3i template for editing
	title ='EQ3NR input file name= ' 							#	line 1 old header
	newtitle = title + file_name[w][:-3] 						#	line 1 new header
	




#################################################################
#################################################################


if __name__=='__main__':
    main()






