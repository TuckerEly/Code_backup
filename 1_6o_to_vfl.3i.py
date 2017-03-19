########################################################################################################################################
###########    1_6o_to_VFL.3i.py    To generate solidless data0.tdf vent fluids from dat0.tde   .6o fies   #############################
####################################      Tucker Ely, 25 DEC 2016   -   By Dillons Bedside    ##########################################
########################################################################################################################################
###  	(1) Read output of equalibrated 6o.
###		(2) If basis not involved in redox:
### 			(2a) Determine total molal Values for all basis species needed to recreate the fluid.
###				(2b) For species involved in redox: Determine totals for each oxidation state and enter
### 				enter as new basis in e-limited basis set
###
###
###   	NOTE: 	The following searches a 6o file set with the following print opions in the header. It may not work with other
### 			print configurations as needed informaiton may be missing.
###
###   				 iopr1-10=     0    0    1    1    1    1    0    1    1    0
###					iopr11-20=     1    0    0    0    0    0    0    0    0    0
###
###		Step 1: fill in basis_search [ ] with elements used in redox reactions
###		Step 2: Check 3o template for accurate values other than fep and t
###		Step 3: run code
###					1 			

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
template = '/Users/tuckerely/Google Drive/Ship/C-DEBI-DCO/EQ36_Test/2_6o_to_vfl.3i/3i_vfl_template.txt'

tdf_file = "Data0_7.tdf.R1" 											#	e_restricted data0 file

run_ID = 1

### All elements involved in redox reactions
basis_search = ["Fe","Fe3i"] 											#	 Reg edit syntax required


"""	EDIT ABOVE	"""

#################################################################
########################    Main    #############################
def main():

	"""Read in all 6o files"""
	file_name, file_list = read_6o()
	### Eventual index of files for which Zi did not reach maximum (failed)
	fail_list = [] 														#	for failed files



	"""Load list of all basis species and species with multiple oxidation states"""
	### This only happens once, for the whole set, an dis based on data0.tdf
	### Load species list from e_inhibited data0 file	
	tdf = open(directory + '/' + tdf_file, "r")
	tdf_lines = tdf.readlines()
	tdf.close()
	### Determine what elements are asssociated with which basis species from data0
	basis_list_tot = generate_basis_species_list(tdf_lines) 			#	needed to determine the basis species called in 3i file for each element grabbed from the 6o
	### Determine what aqueous species are built for each basis species
	basis_composite_list = [] 											#	this is the list of lists (aqueous species linked a specific basis species) Same index as the basis_search list.
	for item in basis_search:
		species_list = generate_redox_species_list(tdf_lines, item) 	#	find all non-solids which contain the basis element (item)
		basis_composite_list.append(species_list)   					#	Add item species list to 




	
	"""Process all 6o files"""
	### Cycle through 6o files
	w = 0 																#	.6o file index
	while w < len(file_name): 											#	cycle through all 6o files, producing one 3i for each if zi is reached in 6o
		### Load 6o file as o
		o = open(file_name[w], "r") 									#	Read current file
		o_lines = o.readlines() 										#	stored 6o file lines
		o.close()


		###Grab needed info from 6o
		fail_trigger, t, pH, fO2, h2o, six_O_basis, basis_composite_list = grab_6o_piecies(o_lines, file_name[w], fail_list, basis_composite_list)
		

		### abort if fail trigger tripped by zi not reaching zi_max
		if fail_trigger == 1:
			w += 1
			continue


		### Prepare the new basis species list for .3i file from combination of .6o and the basis_search items listed by user
		new_basis_value = prep_species_list(basis_search, basis_list_tot, basis_composite_list, six_O_basis)


		### Write to new .3i file
		write_3i(w, file_name, new_basis_value, t, fO2, pH)
		

		w += 1  									#	onto the next 6o file

#################################################################
#####################    Functions     ##########################

def read_6o():
	### this function can find all .6o files in all downstream folders, hence the additio of file_list
	file_name = [] 						#	file names
	file_list = []						# 	file names with paths
	for root, dirs, files in os.walk(directory):# this works. 
		for file in files:
			if file.endswith(".6o"):
				file_name.append(file)
				file_list.append(os.path.join(root, file)) 					
	
	return file_name, file_list

def generate_basis_species_list(tdf_lines):

	basis_list_tot = [] 												#	basis species in data0, [species, element component, sto]

	### build list of all basis species in data0, with associated element
	x = 0 																#	start at top of data0
	while not re.findall('^basis species', tdf_lines[x]): 				#	Move to top of basis species
		x += 1
	x += 1 																#	move to first casis (h2o) 
	while x < len(tdf_lines):
		if re.findall('\+-----', tdf_lines[x]): 						#	This should find the fist as well, as the x indexed begins on top of the +---- above h2o
			x += 1
			a = str(tdf_lines[x])
			b = a.split()
			spe = b[0] 													#	basis species name of this block
			### Exit upon hitting the auxiliary set
			if re.findall('^auxiliary', spe):
				return basis_list_tot
				break
			### deal with H2O, H+, and O2(g)
			elif re.findall('^[hH]2[oO]$', spe) or re.findall('^[hH]\+$', spe) or re.findall('^[oO]2\([gG]\)$', spe):							#	if it is h2o
				ele = 'NONE' 											#	default to an element of none
				sto = 'NONE'
				basis_list_tot.append([spe, ele, sto])
				x += 1 													# 	back to while loop
			### All other basis species
			else:														#	seach for elements other than H or O in remainging basis species
				x += 1	
				while not re.findall('\+-----', tdf_lines[x]): 			#	stop if next basis is reached
					if re.findall('chemical elements', tdf_lines[x]): 	#	data of interest is ononeline bloew this line
						x += 1
						a = str(tdf_lines[x])
						b = a.split()
						c = [item for item in b[1::2] if item != 'H' and item != 'h' and item != 'O' and item != 'o'] 	#	for element that is not O or H, should only be 1
						d = re.search(c[0], a)
						e = d.start()
						sto = a[(e - 7) : (e - 1)] 						#	grab element stochiometry (dont knwo if i will end oup needing this)
						basis_list_tot.append([spe, c[0], sto]) 			#	append this basis species info to list
						x += 1
					else:
						x += 1
		else:
			x += 1

	return basis_list_tot

def generate_redox_species_list(tdf_lines, item):
	###	Search the lines of data0.tdf.Rxxx
	ss = ' ' + item 									#	search string
	species_list = []

	### This returns a list of all [sto, species] pairs that contina the element (item) being called form the basis_search
	x = 0
	while x < len(tdf_lines):			
		
		### Break out of loop if soilds are reached
		if re.findall('^solids', tdf_lines[x]):
			return species_list


		### find iten (basis species) in basis, aux, and aqueous species
		if re.findall(ss + '[ \n]', tdf_lines[x]):	 				#	must catch elements ending in a space, or a newline hence [ \n]
			### stochiometry = previous 6 spaces (how though, reg expr just found the  line, not the space)
			s = str(tdf_lines[x])
			a = re.search(ss, s)
			b = a.start()
			sto = s[(b - 6) : (b - 1)]
			### now that I have grabbed stoiciometry, walk back to beginning of block to find the name of the species
			while not re.findall('\+-----', tdf_lines[x]):
				x -= 1
			x += 1
			a = str(tdf_lines[x])
			b = a.split() 											#	need to split to get away from newline.
			spe = b[0]												#	species name
			### Move to next block to aviod hitting the same item twice
			while not re.findall('\+-----', tdf_lines[x]):
				x += 1
			species_list.append([float(sto),spe]) 					#	Add Stoichiometry + Species Names
			

		else:
			x += 1

	return species_list

def grab_6o_piecies(o_lines, file_name, fail_list, basis_composite_list):
	
	t = 'NaN'
	pH = 'NaN'
	fO2 = 'NaN'
	h2o = 'NaN'
	six_O_basis = []


	### Zi and Zi Max determination
	### search from beginning untill zimax is found, then stop	
	x = 0
	while x < len(o_lines):			
		if re.findall('^     zistrt=', o_lines[x]):	
			a = str(o_lines[x])
			b = a.split()
			zi_max = float(b[3])
			break
		else:
			x += 1
	### then search from the bottom up unill 
	### check zi completion
	x = len(o_lines) -1
	while x > 0: 							#	count down from end of file
		if re.findall('^ {21}reaction progress {8}=', o_lines[x]):
			a = str(o_lines[x])
			b = a.split()		
			zi_reached = float(b[3])	
			break 							#	this should leave x where i want it, on the line number of the last reaction step
		else:
			x -= 1		
	### Exit if zi_max is not reached
	fail_trigger = 0
	if zi_reached != zi_max:
		print 'Zi_max not reached on ' + file_name +  ': zi_reached = ' + str(zi_reached) 
		fail_list.extend(file_name)
		fail_trigger = 1		
		return fail_trigger, t, pH, fO2, h2o, six_O_basis, basis_composite_list		#	return sends me back out of this function
	### Zi max has been reached, continue grabbing needed info from the file.






	### Grab needed components from last step (x should already be at the correct index)
	while not re.findall('^ {3}species {16}moles {8}grams', o_lines[x]): 				#	Now count back to the end of the file
		### Grab t
		if re.findall('^ {21}temperature {4}=', o_lines[x]):
			a = str(o_lines[x])
			b = a.split()		
			t = "%.4f" % float(b[2]) 							#	needs to be formatted for 4 decimal places.
			x += 1
		### Grab pH
		elif re.findall('^     modified nbs ph scale ', o_lines[x]):
			a = str(o_lines[x])
			b = a.split()		
			pH = -float(b[4])
			x += 1
		### Grab log fO2
		elif re.findall('^          log oxygen fugacity =', o_lines[x]):
			a = str(o_lines[x])
			b = a.split()		
			fO2 = float(b[4])
			print fO2
			x += 1
		### Grab log aH2O
		elif re.findall('^               log activity of water =', o_lines[x]):
			a = str(o_lines[x])
			b = a.split()		
			h2o = "%.5f" % float(b[5]) 							#	needs to be formatted for 5 decimal places.
			x += 1
		### Grab aq element totals
		elif re.findall('^      --- element totals for the aqueous phase ---', o_lines[x]):	
			x += 4
			while not re.findall('^\n', o_lines[x]):
				a = str(o_lines[x])
				b = a.split()
				six_O_basis.append([str(b[0]),str(b[2])]) 			#	[name, molal conc.]
				x += 1		
		else:
			x += 1

	


	### Should now be at the top of aq species.
	x += 2
	top = x 										#	mark top of aq species
	ba_el = 0										#	index within basis_composite_list
	### make basis composite list re compatible
	for item in basis_composite_list:
		for species in item:
			if '+' in species[1]:
				species[1] = species[1].replace('+', '\+') 
	while ba_el < len(basis_composite_list): 								#	iterate through basis elements in basis_search / basis_composite_list
		aq = 0 																# 	index within basis_composite_list[ba_el]
		while aq < len(basis_composite_list[ba_el]):						#	iterate trough all aq species of a given basis element in basis_search
			if re.findall('^   ' + str(basis_composite_list[ba_el][aq][1]) + ' ', o_lines[x]): 	#	the trailing space should keep partial hits from happening.
				a = str(o_lines[x])
				b = a.split()
				basis_composite_list[ba_el][aq].append(b[3])
				x = top 													# 	send index back to top of aq species to look for next aq
				aq += 1
			else:
				x += 1
		ba_el += 1

	return fail_trigger, t, pH, fO2, h2o, six_O_basis, basis_composite_list

def prep_species_list(basis_search, basis_list_tot, basis_composite_list, six_O_basis):

	### build new_basis_value for new basis sp. in basis_search from 6o files aq section
	new_basis_value = []											#	list of eventual  basis_sepcies names and values: same index as basis_search
	ba = 0 															#	index for basis_search
	while ba < len(basis_search):
		summ = 0.0 													#	float for summing together the basis composite species
		s = 0 														#	index for composite species within a given basis element
		while s < len(basis_composite_list[ba]):
			summ = summ + float(basis_composite_list[ba][s][2])		#	sum togeher of the of the species which contain the basis_search element
			s += 1
		summ = str(summ)
		if 'e' in summ:												#	make e in 1.000000e-01 uppercase
			summ = summ.replace('e', 'E') 
		new_basis_value.append(summ) 								#	add tot he basis_elemet the total molol conc. of its summed composite species
		ba += 1


	### Conbine basis_search elements with their derived values
	s = 0
	while s < len(basis_search):
		new_basis_value[s] = [basis_search[s], new_basis_value[s]]
		s += 1 


	### build new elements list from basis_search elemetns, and remaining new_basis_values form 6o (Excluding H and O). 
	s = 0
	while s < len(six_O_basis):
		if six_O_basis[s][0] not in basis_search and six_O_basis[s][0] != 'O' and six_O_basis[s][0] != 'H':  		#	start with basis_searhc elements, then add all item in 6o element, if they are not already present in the basis search
			new_basis_value.append(six_O_basis[s])
			s += 1
		else:
			s += 1
	
	### replace the elemtns of new_bais_values with the correct associated basis_species, for input in .3i file
	r = 0
	while r < len(new_basis_value):
		s = 0
		while s < len(basis_list_tot):
			if new_basis_value[r][0] == basis_list_tot[s][1]:
				new_basis_value[r][0] = basis_list_tot[s][0]
				s += 1
			else:
				s += 1
		r += 1
	return new_basis_value

def write_3i(w, file_name, new_basis_value, t, fO2, pH):

	new_file = file_name[w][:-3] + '_vfl.3i' 					#	new file title
	shutil.copyfile(template, new_file) 						#	copy 3i template for editing
	title ='EQ3NR input file name= ' 							#	line 1 old header
	newtitle = title + file_name[w][:-3] 						#	line 1 new header
	for line in fileinput.input(new_file, inplace=1):			#	lines to replace
		line = re.sub(title, newtitle, line.rstrip())
		line = re.sub('     tempc=', '     tempc=' + (13 - len(t))*' ' + str(t), line.rstrip())   								#	print t:  13 spaces total to get the end of the 4 decimals places alotted to t
		line = re.sub('       fep=               ', '       fep=' + (15 - len(str(fO2)))*' ' + str(fO2), line.rstrip()) 		# 	print fep: 15 spaces
		print(line)												#	write into file
	fileinput.close()											#	


	### reopen and write from the end of the file
	with open(new_file, 'a') as build:
		line_1 = 'data file master species= '
		line_2 = '   switch with species='
		line_3 = '   jflag= 0     csp= ' 						#	0 for total_molality

		for item in new_basis_value:
			species_add = line_1 + item[0] + '\n' + line_2 + '\n' + line_3 + str(item[1]) + '\n'
			build.write(species_add) 

		pH_add = line_1 + 'h+' + '\n' + line_2 + '\n' + '   jflag= 16    csp= ' + str(pH) + '\n'		#	Add pH
		build.write(pH_add)
		end = 'endit.'
		build.write(end) 								# 	Add ending and close
	build.close()


#################################################################
#################################################################


if __name__=='__main__':
    main()






