########################################################################################################################################
###################################################     Search_6o.py    ################################################################
##################    to grab 6o output phases as needed    Tucker Ely, 17 March 2017        ###########################################
########################################################################################################################################
###		This code is modified from the origional  6o.any_species3.py. It has been updated in its organization
###		And has been modified to handle version 8 outputs.
### 
###
###
###
###




from collections import defaultdict

from string import digits
from operator import itemgetter
from collections import OrderedDict
 


import os, re, sys, shutil, fileinput, getopt, argparse, csv, copy
import numpy as np
import pandas as pd
from itertools import *
from collections import *
from subprocess import * 
from time import *





""" EDIT BELOW """

pwd = os.getcwd()

rock_stats = '1_Dalton_Segment_catalog_total.csv'			#	file continaing all of the rock stats per segment


aq_grab = 0 						# 1 = yes   0 = no
aq_sp = ['SO4-2','HCO3-','CO2,AQ', 'Na\+','Mg\+2','Fe\+2','Fe\+3','Ca\+2', 'METHANE,AQ'] 	# List desired
# 'SO4-2','CO,aq','CO2,aq','HCO3-','CO3-2','O2,aq','HS-','H2S,aq','H2,aq','METHANE,AQ','FORMATE,AQ','ACETATE,AQ','MTHEANOL,AQ','ACETIC-ACID,AQ','FORMALDEHYDE,AQ','ACETALDEHYDE,AQ','ETHANE,AQ','SiO2,aq','HPO4-2','Fe\+2','Fe\+3','Mn\+2','MnSO4,aq'


s_grab = 0 



""" EDIT ABOVE """
#fitting coefficients, eq# denotes order

# 2-350K and 500bar
neutral_ph_eq0 = 7.34334 
neutral_ph_eq1 = - 0.0193815
neutral_ph_eq2 = 0.0000796789 
neutral_ph_eq3 = - 1.75343E-7
neutral_ph_eq4 = 1.83503E-10

fmq_eq0 = -87.9119 
fmq_eq1 = 0.335374
fmq_eq2 = - 0.000959826
fmq_eq3 = 1.84084E-6
fmq_eq4 = - 1.57724E-9


#################################################################
########################    Main    #############################

def main(argv):
	project_folder_name, six_o_type, retrieve_type = process_arguments(argv, pwd) 		#	process the arguements given.
	project_path = pwd + '/' + project_folder_name
	six_o_folder = project_path + '/' + six_o_type + '_6o'
	file_name, file_list = read_inputs('6o', six_o_folder)				#	grab 6o's to be analyzed
	
	visual_folder = project_path + '/' + six_o_type + '_search' 		#	mass output folder
	mk_check_del_directory(visual_folder) 								#	make that folder


	os.chdir(six_o_folder) 												#	move into the folder with the 6o files
	
	### build lists for final
	if retrieve_type == 'final':
		# t_global = [], ph_global = [], modified_ph_global = [], neutral_ph_global = [], fo2_global = [], modified_fo2_global = [], aq_array_global = []
		global_add_list = []


	w = 0
	while w < len(file_name): 											#	cycle through all files
		print 'Processing ', file_name[w]
		lines = grab_lines(file_name[w]) 								#	6o file lines

		### Grab whole 6o path
		if retrieve_type == 'path':
			print 'path retrieval not yet supported. Use "final".'
			sys.exit()
			# t = grab(lines) 										
			# modified_fo2, fmq, fo2 = grab_fo2(lines)
			# ph, modified_ph, neutral_ph, eh, pe = grab_ph_eh_pe(lines)



		### Grab jsut the vent fluid at teh final step. If this is chosen, the oupt will be an amended version of 
		### the rock stats data.
		if retrieve_type == 'final':
			t, ph, modified_ph, neutral_ph, fo2, modified_fo2, aq_array = grab_6o_piecies(lines)
			### Append to global
			global_add_list.append( [t, ph, modified_ph, fo2, modified_fo2] + aq_array)

		w += 1



	if retrieve_type == 'final':
		print_final_data(rock_stats, global_add_list, file_name, project_path, project_folder_name)

#################################################################
#####################   Primary Functions     ###################

def process_arguments(argv, pwd):
	### Grab Arguments and chack for their accuracy.
	
	try: 														#	test for the correct arguemnts in command line
		opts, args = getopt.getopt(argv,"h",["help"]) 								
	except getopt.GetoptError: 									#	if an unrecognized flag is passed, an exception is raised. right now only -h is recognized
		print 'Unrecognized argement.'
		help() 													#	Show user what is expected to input
		sys.exit(2) 											#	2 here is the exit code and means (command line syntax errors)


	### Check for exceptions and errors in arguments.
	for opt, arg in opts: 										#	If user called help, then give it to them
		if opt in ("-h", "--help"):
			help()
			sys.exit()
	if len(args)<3:												#	if fewer than 5 arguments are given, FIX THAT SHIT !
		help() 													#	FIX THAT SHIT ! !  Notify user of correct arguments


	### Argument 1: project folder in pwd
	if not os.path.isdir(pwd + '/' + args[0]): 					#	Check for existance of folder containing the inputs
		print 'Project directry not present in current folder.'
		sys.exit()
	project_folder_name = args[0] 								#	set project folder name


	### Argument 2: dslop file to be used
	if args[1] != 'mix' and args[1] != 'rxn':					
		print 'Options are either rxn, or mix, currently'
		sys.exit() 	
	six_o_type = args[1] 	

	### Argument 3: dslop file to be used
	if args[2] != 'path' and args[2] != 'final':					
		print 'Options are either path, or final, currently'
		sys.exit() 	
	retrieve_type = args[2] 										


	return project_folder_name, six_o_type, retrieve_type	#	Grab all that shit !

def grab_t(lines):
	t = []
	x = 0
	while x < len(lines): 		# builds t and zi in master file,   per .6o
		if re.findall('^ Temperature=', lines[x]):
			a = str(lines[x])
			b = a.split()
			t.append(float(b[1]))
			x += 1
		else:
			x += 1      
	return t

def grab_fo2(lines):
	fo2 = []
	x = 0
	while x < len(lines):
		if re.findall('^ {16}Log oxygen fugacity=', lines[x]): 
			a = str(lines[x])
			b = a.split()
			fo2.append(float(b[3])) 	# 	log value 	
			x += 1
		else:
			x += 1  

	fmq = [(fmq_eq0 + fmq_eq1*i + fmq_eq2*i**2 + fmq_eq3*i**3 + fmq_eq4*i**4) for i in  t]
	modified_fo2 = [0]*len(fmq)
	x = 0
	while x < len(fmq):
		modified_fo2[x] = fo2[x] - fmq[x] 			#	This will mean that - values are more acidic thn neutral, and + more basic
		x += 1 

	return modified_fo2, fmq, fo2

def grab_ph_eh_pe(lines):
	ph = []
	eh = []
	pe = []
	x = 0
	while x < len(lines): #	Modified nbs ph scale
		if re.findall('^ NBS pH scale', lines[x]): 
			a = str(lines[x])
			b = a.split()
			ph.append(float(b[3])) 
			eh.append(float(b[4]))
			pe.append(float(b[5]))
			x += 1
		else:
			x += 1
	
	# incorporat neutral ph change. modified then becomes the pH as if 7 were always neutral
	neutral_ph = [(neutral_ph_eq0 + neutral_ph_eq1*i + neutral_ph_eq2*i**2 + neutral_ph_eq3*i**3 + neutral_ph_eq4*i**4) for i in  t]
	modified_ph = [0]*len(neutral_ph)
	x = 0
	while x < len(neutral_ph):
		modified_ph[x] = ph[x] - neutral_ph[x] 			#	This will mean that - values are more acidic thn neutral, and + more basic
		x += 1 



	return ph, modified_ph, neutral_ph, eh, pe

def grab_6o_piecies(lines):
	
	### defaults to be replaced
	t = 'NaN'
	pH = 'NaN'
	fo2 = 'NaN'
	aq_array = [0]*len(aq_sp)

	### search from the bottom up unill last xi is encountered
	x = len(lines) -1
	while x > 0: 							#	count down from end of file
		if re.findall('^ {20}Xi=', lines[x]):
			a = str(lines[x])
			b = a.split()		
			xi_reached = float(b[1])	
			break 							#	this should leave x where i want it, on the line number of the last reaction step
		else:
			x -= 1		

	### Grab needed components from last step (x should already be at the correct index)
	while not re.findall('^ {16}--- Distribution of Aqueous Solute Species ---', lines[x]): 				#	Now count back to the end of the file
		### Grab t
		if re.findall('^ Temperature=', lines[x]):
			a = str(lines[x])
			b = a.split()		
			t = "%.4f" % float(b[1]) 							#	needs to be formatted for 4 decimal places.
			x += 1
		### Grab pH
		elif re.findall('^ NBS pH scale', lines[x]):
			a = str(lines[x])
			b = a.split()		
			ph = float(b[3])
			x += 1
		### Grab log fO2
		elif re.findall('^ {16}Log oxygen fugacity=', lines[x]):
			a = str(lines[x])
			b = a.split()		
			fo2 = float(b[3])		
			x += 1
		### Grab aq element totals
		else:
			x += 1



	### Grab aq species of interest, as identified in aq_sp
	### Should be at teh top of the aq species given the previous while loops end point
	x += 4
	top = x 										#	mark top of aq species
	sp_n = 0										#	index within aq_sp
	### Find the end of the aq sections
	while not re.findall('      --- Major Species by Contribution to Aqueous Mass Balances ---', lines[x]):
		x += 1
	end_aq_line = x

	### search for aq species.
	x = top																	#	put on top of the aq stack
	while sp_n < len(aq_sp): 												#	iterate through basis elements in basis_search / aq_sp
		x = top 															#	put back on top of the aq stack
		while x < end_aq_line: 												#	stop if species is not found id aq block
			if re.findall('^ ' + aq_sp[sp_n] + ' ', lines[x]): 				#	the trailing space should keep partial hits from happening.
				a = str(lines[x])
				b = a.split()
				aq_array[sp_n] = b[1]										#	1 grabs molality
				x = end_aq_line 											#	if species found, send to end of section
			else:
				x += 1
		sp_n += 1




	### Neutrals comparison
	
	t = float(t)
	neutral_ph = neutral_ph_eq0 + neutral_ph_eq1*t + neutral_ph_eq2*t**2 + neutral_ph_eq3*t**3 + neutral_ph_eq4*t**4
	modified_ph = ph - neutral_ph
	fmq = fmq_eq0 + fmq_eq1*t + fmq_eq2*t**2 + fmq_eq3*t**3 + fmq_eq4*t**4
	modified_fo2 = fo2 - fmq	
	t = str(t)

	return t, ph, modified_ph, neutral_ph, fo2, modified_fo2, aq_array

def print_final_data(rock_stats, global_add_list, file_name, project_path, project_folder_name):

	os.chdir(pwd) 									#	need to be here to load the rock db
	seg_names = []									#	list of segment names as they are ordered in the rock stat csv
	seg_names_local = []							#	lsit of segment from this run with information being added to the rock stat file

	for file in file_name:
		seg_names_local.append(file[:(file.find('_'))])					#	Grab the segment name which precedes the  _  in the six_o_name

	lines_as_lists = [] 												#	to add rock stats csv data, list of line
	with open(rock_stats, 'rU') as fin:
		fin.seek(0)
		reader = csv.reader(fin)
		for row in reader:
			lines_as_lists.append(row)
	pre_addition_row_len = len(lines_as_lists[0])						#	what is the row length of the origional rock array
	
	old_col_names = lines_as_lists[0][1:]								#	grab old col names	
	lines_as_lists = lines_as_lists[1:]									#	remove that first head row (list)
	seg_add_info = []													#	list of new row info to be added
	y = 0 																#	y = the index for a given segment in all fo the info list to add (t, ph, etc.)
	while y < len(seg_names_local):
		seg_add_info.append(global_add_list[y])
		y += 1

	s = 0
	while s < len(seg_names_local):										#	seg_name_local and global_add_list have same index
		for item in lines_as_lists:
			if item[0] == seg_names_local[s]:
				item = item.extend(global_add_list[s])					#	if line in lines_as_lists[item[0]] == segname[s], then copy the global add info the the end of that line
		s += 1
	

	### bring list sof list to = dimensions, by adding 0's to missing data. this will aloows for datafram conversion
	blank_add = [0]*len(global_add_list[0])								
	for item in lines_as_lists:
		if len(item) == pre_addition_row_len:
			item = item.extend(blank_add)

	### remove the first position (seg name) for all lists within lines_as_lists. This will instad be added using row names in the dataframe
 	row_names = [line[0] for line in lines_as_lists]					#	use index 0 for pandas printing
 	lines_as_lists = [line[1:] for line in lines_as_lists] 				#	now get rid of it, as it will be in the index, not the 'array'

	aq_sp_re = copy.deepcopy(aq_sp) 	 								#	deepcopy is needed here to keep alterations downstream from affect both lists. As 						#	copy an re compatible version
	aq_sp_re = [item.replace('\+', '+') for item in aq_sp_re] 			#	make basis composite list re compatible.
	new_col_names = ['T','ph','mod_ph','fo2','mod_fo2'] + aq_sp_re
	col_names = old_col_names + new_col_names

	os.chdir(project_path)
	### Print new dataframe to project folder (not pwd)
	out_name =  project_folder_name + '_rock_stats.csv'
	df = pd.DataFrame(lines_as_lists, index=row_names, columns=col_names) # INDEX = row titles
	df.to_csv(out_name)


#################################################################
#####################    Helper Functions     ###################

def help():
	### Displayed to user of any of the arguments are enetered wrong
	print '2 arguements needed:' 
	print '		(1) Project folder name. Must be in local directory'
	print '		(2) rxn or mix'
	print '		(3) whole "path", or "final" step'

def read_inputs(file_type, location):
	### this function can find all .6o files in all downstream folders, hence the additio of file_list
	file_name = [] 						#	file names
	file_list = []						# 	file names with paths
	for root, dirs, files in os.walk(location):# this works. 
		for file in files:
			if file.endswith(file_type):
				file_name.append(file)
				file_list.append(os.path.join(root, file)) 					
	
	return file_name, file_list

def mk_check_del_directory(path):
	###  This code checks for the dir being created, and if it is already present, deletes it, before recreating it
	if not os.path.exists(path): 							#	Check if the dir is alrady pessent
		os.makedirs(path) 									#	Build desired output directory
	else:
		shutil.rmtree(path) 								#	Remove directory and contents if it is already present.
		os.makedirs(path)

def mk_check_del_file(path):
	###  This code checks for the file being created/moved already exists at the destination. And if so, delets it.
	if os.path.isfile(path): 								#	Check if the file is alrady pessent
		os.remove(path) 									#	Delete file

def grab_lines(file):
	f = open(file, "r")
	lines = f.readlines()
	f.close()
	return lines


#################################################################
#####################     Call Main    ##########################


if __name__=='__main__':
    main(sys.argv[1:]) 						

#################################################################
##################   Troubleshooting Notes    ###################






