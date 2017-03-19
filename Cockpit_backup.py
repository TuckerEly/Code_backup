########################################################################################################################################
################      Cockpit.py:  To dynamically call all required subroutines of the geochemical modeling process.     ###############
##################################################      Tucker Ely, 27 Jan 2017       ##################################################
########################################################################################################################################
###  	Curently not included:
###				slop  	--> 	CPRONS 		--> dslop
###
###		Current Limitations:
### 			build_6i_files function currently assumes that the generated 6i files only contain 1 solid reactant
###					comprised of oxides.
###
### 	To include:
###				data0 	--> 	DBCreate 	--> data0.new
### 			data0 	-->		EQPT 		-->	data1  +  data1f
### 			
###				name.3i (en masse) + data1	--> EQ3 -->	name.3o + name.3p (en masse)
### 			
### 				build six_i template(three_p)
###
###				6i_rxn_template(solids db) + sw.3p	--> 1_6i_build.py --> rxn.6i (en masse)   (redone in version 8)
### 				#####  Temporary work around  ###### Issues with where to get * Reaction info for version8 6i's
###					
###					- run 1_6i_build.p on version 7 formate.
###						(1) convert to version 8
###						(2) remove pickup from bottom
###					 	(3) re attach new version 8 SW pickup
### 				
###					- pickup_syntax_fix_1(project_path)
### 					it was discovered that the 3p files used at the base of my 6i's had syntax errors
###						I have built this temporary patch, untill I can determine why eq3 generted these  
###						faulty pickup files.		
### 		
###				rxn.6i (en masse) + data1  -->  EQ6  -->  rxn.6o (en masse) + rxn.6p (en masse) 
###
### 			rxn.6o (enmasse) + data0   -->  6o_to_VFL.3i.py  -->  rock_vfl.3p (en masse)
###
###				rock_vfl.3p (en masse) + tde_sw.3p	--> 1_SW_VFL_mix.6i.py --> mix.6i (en masse)   (redone in version 8)			
###
### 			Include vis function, 
###
### 			Include output to  mathematica for geochemical modeling
### 			INclude output to mathematica for Metabolism modeling
###
###		Use:  	(1) This code currently assumes tde. for SW.3i and initial rxn.6i sets, and .tdf for mixing.
###				(2) Code currently runs both tde and tfd through same CON file. It doesnt need to be this way 
###				, as they only need to match up at their highest temperatures.
###
###


import os, re, sys, shutil, fileinput, getopt, argparse, csv, copy
import numpy as np
import pandas as pd
from itertools import *
from collections import *
from subprocess import * 
from time import *



"""	EDIT BELOW	"""

pwd = os.getcwd()
slop_db_p = pwd + '/1_SLOP/working_db/' 				#	Storage for slop and dslop files
data0_storage_p = pwd + '/3_data0/working_db/' 			#	For long term stroage of data0s
dbc_exe_p = pwd 	#+ '/4_DBCreate_v2/' 				#	DBCreate location for running
eq_home = '/Users/tuckerely/EQ3_6v8.0a' 				#	EQ3.6 software home


six_i_v8_template_empty = '6i_v8_template.6i' 			#	template to be used in conjunction with sw pickup on db.
six_i_v7_template_empty = '6i_v7_template.6i' 	
six_i_build_code = '1_6i_build.v2.py' 					#	note correct version number

meta_rxn_list = '1_rxn_list' 							#	The list of metabolic reactions
e_list_file = '1_rxn_list_e-transfered'					#	the list of e-transfered associated with the reaction sin the aboe rxn list


### Thermo constants for  A = RTln(K/Q)
R = 8.31446 / 4.184 				#	with calory conversion (cal K-1 mol-1)
euler = 2.7182818284 
kelvin = 273.15 					#	at 0 C
cal_convert = 0.23900574



"""	EDIT ABOVE	"""


eq_env = eq_home + '/scripts/eq36cfg' 					#	To set environemntal variables for eqpt/eq3/eq6 running
eq_db_p = eq_home + '/db' 								#	Location of data0/data1/dat1f for eqpt/eq3/eq6 running
eqs_exe_p = '/Users/tuckerely/EQ3_6v8.0a/bin' 			#	executables folder with runeqpt/runeq3/runeq6/xcon3/xcon6

# python Cockpit.v1.py test_project dslop15_v2.1_Ti.dat 500_0-300 j43_tde.3i 1_Gale_DB_seg_ME.csv

# python Cockpit.v1.py test_project dslop15_v2.2_Ti.dat 500_0-300 j43_tde.3i 1_Gale_DB_seg_ME_Ti_0.csv

# python Cockpit.v1.py test_project dslop15_v2.2_Ti.dat 500_0-300 peters_bsw_tde.3i 1_Gale_DB_seg_ME_Ti_0.csv

#################################################################
########################    Main    #############################

def main(argv):

	# python Cockpit.v1.py test_project dslop15_v2.1_Ti.dat 500_0-300 j43_tde.3i

	project_folder_name, dslop_file, con_file, initial_SW_file, db, old_data0_tde, old_data0_tdf = process_arguments(argv, pwd) 	#	pricess user arguments
	project_path = pwd + '/' + project_folder_name + '/' 						#	Full path to created project folder


	
	""" Call DBCreate """
	# call_DBCreate('tde', old_data0_tde, dslop_file, con_file, project_path) 		#	build data0.tde at conditions in con file
	# call_DBCreate('tdf', old_data0_tdf, dslop_file, con_file, project_path)		#	build data0.tdf at conditions in con file


	""" Call EQPT """
	# call_EQPT('tde', project_path)
	# call_EQPT('tdf', project_path)


	"""   Stage 1:	Call EQ3 on initial water file.   """
	### Add run title and text into header of file.   chang call rutine to the wait function
	# call_EQ3('tde', project_path, initial_SW_file)


	"""    Stage 2:	 Build rxn.6i files. run v7 conversion error fixes   """
	### Temporarily obsolete, untill i cen figure out what information is in the new 6i files. For now they are constructed as version7, and then converted, to let XCON6 do the work
	# six_i_file_names = build_6i_files('tde', project_path, project_folder_name, initial_SW_file, '6i_v7_template.6i', db) 								#	build version 7 copy of all 6i files from rock database
	
	### 6i spot fixes related to xcif6 errors and v7-v8 conversion faults
	# six_i_syntax_fix_1(project_path, six_i_file_names)			#	add '    qgexsh=        F' to 3p fluid reactant.
	# six_i_syntax_fix_2(project_path, six_i_file_names) 			# 	jtemp to 1 from 2
	# six_i_syntax_fix_3(project_path, six_i_file_names) 			#	replaces mistranslations in iopt1-10
	# six_i_syntax_fix_4(project_path, six_i_file_names) 			#	replaces mistranslations in iopt11-20
	# six_i_syntax_fix_5(project_path, six_i_file_names) 			#	remove Ti from rock 
	# six_i_syntax_fix_6(project_path, six_i_file_names) 			#	remove P from rock 
	# six_i_syntax_fix_7(project_path, six_i_file_names) 			#	remove Mn from rock 

	"""   Stage 3:	Run all 6i files   """
	# six_i_file_names, six_i_file_lists = read_inputs('6i', project_path + 'rxn_6i') 						# 	Redundnat addition, inse this stage wants to be run as standalone from those above it.
	# run_6i_rxn('tde', project_path, six_i_file_names)														#	run all 6i files created above
	# check_xi_completion(project_path)
	# generate_final_rock_stats(project_path, project_folder_name)


	"""   Stage 4: 	Process vent fluid and build mix files with separated SW   """
	# basis_ele, redox_active_ele, basis_composite_list = determine_redox_ele_and_sp('tdf', project_path)		# 	determin redox active elements, and a list of all aq species using that element.  For e-separated data0
	# build_vfl_3i_redox_sep_files(project_path, basis_ele, redox_active_ele, basis_composite_list)			#	build e-separated vfl.3i files from rxn.6o files
	# run_vfl_3i(project_path, 'tdf') 																		# 	Run all vfl.3i files
	# check_vfl(project_path)																				#	check for errors and xi completion


	"""   Stage 5:	Construct and run e-separated sw   """
	# sw_tdf_3i = build_sw_tdf(project_path, initial_SW_file,  basis_ele, redox_active_ele, basis_composite_list) # 	build e-separated SW.
	# sw_tdf_3i = 'peters_bsw_tde_sw_tdf.3i'
	# call_EQ3('tdf', project_path, sw_tdf_3i) 																	#	Call EQ3 on sw_tdf.3i


	"""   Stage 6:	Build and run vfl/sw mixes.   """
	# build_mix_files(project_path, sw_tdf_3i, 'mix_template.6i')
	# six_i_file_names, six_i_file_lists = read_inputs('6i', project_path + 'mix_6i') 	# 	Redundnat addition, inse this stage wants to be run as standalone from those above it.
	# run_6i_mix('tdf', project_path, six_i_file_names) 								#	Run all mix.6i files
	# check_mix_xi_completion(project_path)					 							#	check completion


	##############################  Stage  7  Forward  ##############################


	"""   Stage 7:  grab species info for metabolic reactions.   """
	t_slices = 54 			#	number of xi/t steps to ta ke from each 6o file. this includes step 1 (pure vfl, + however many steps counting from the back)
	six_o_file_names, six_o_file_lists = read_inputs('6o', project_path + 'mix_6o') 	# 	grab mix outputs. These need to have been processed by the above check_mix_xi_completion()

	""" Xi/T determination """
	### (1) determine Xi path steps to be grabbed and their associated T; Cycle through all files, varifying correct Xi Steps, removing all other steps (this should cut down on file size to search through at later steps)
	average_xi, ln_n_array, average_t, average_p = build_xi_t_path(project_path, six_o_file_names, t_slices)

	""" run supcrt  / grab K"""
	### (3) build universal con file
	# supcrt_con_file = build_con_file(project_path, project_folder_name, average_t, average_p)
	### (4) run all supcrt.rxn files
	# call_supcrtgrid(project_path, dslop_file, old_data0_tde, supcrt_con_file)
	
	### (5) grab all needed K data from supcrt.tab files
	rxn_file_names, rxn_list, rxn_k_value_list, rxn_ID_list = load_supcrt_rxn_and_k_data(project_path)
	e_list = grab_e_transfered(e_list_file)																#	grab the associated e- tranfered n meta_reactions. index is the same as rxn list

	""" build Q values """
	### (6) grab unique speceis from .tab files (built already)
	total_unique_species, unique_nonsolids, unique_solids = check_aq_species(rxn_list)

	### (7) grab species activities, and concentrations. Determine Q, and limiting
	species_activity_array, species_conc_array = grab_rxn_species_from_6o(project_path, six_o_file_names, unique_nonsolids, ln_n_array, average_t)
	Q_value_array = determine_Q_values(project_path, six_o_file_names, average_t, species_activity_array, rxn_list, rxn_ID_list, unique_nonsolids, unique_solids)
	solid_conc = solid_conc_determination(average_xi, unique_solids)										#	determine effetive solids concentrations for bottom seawater, by diluting best guess seawater values for those minerals.
	limiting_array = determine_limiting_species(project_path, six_o_file_names, average_t, rxn_ID_list, rxn_list, species_conc_array, unique_nonsolids, unique_solids, solid_conc)

	""" build Affinities and Vis"""
	### (8) Determine Affinity and E/ml
	rock_affinity_array, rock_energy_array = determine_affinity_and_Energy(project_path, six_o_file_names, rxn_k_value_list, average_t, e_list, rxn_ID_list, rxn_list, limiting_array, Q_value_array)
	array_to_csv(project_path, 'rock_energy_array', rock_energy_array, six_o_file_names, average_t, rxn_ID_list)  				#	(project_path, name, array, index_1, index_2, index_3)		
	array_to_csv(project_path, 'rock_affinity_array', rock_affinity_array, six_o_file_names, average_t, rxn_ID_list)			  	#	(project_path, name, array, index_1, index_2, index_3)		

	### (8.1) Convert the above per file arrays, to per reaction.
	rxn_affinity_array, rxn_energy_array = convert_rock_to_rxn_array(six_o_file_names, average_t, rxn_ID_list, rock_affinity_array, rock_energy_array)
	### (8.2) Tack on rock stats to the bottom of the rxn arrays
	combined_energy_array, combined_affinity_array, combined_index = load_rock_stats(project_path, project_folder_name, six_o_file_names, rxn_ID_list, average_t, rxn_affinity_array, rxn_energy_array)

	### (7) build output sufficient for vis.
	array_to_csv(project_path, 'rxn_energy_array', combined_energy_array, rxn_ID_list, combined_index, six_o_file_names)  			#	(project_path, name, array, index_1, index_2, index_3)		
	array_to_csv(project_path, 'rxn_affinity_array', combined_affinity_array, rxn_ID_list, combined_index, six_o_file_names)		#	(project_path, name, array, index_1, index_2, index_3)		


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
	if len(args)<5:												#	if fewer than 5 arguments are given, FIX THAT SHIT !
		help() 													#	FIX THAT SHIT ! !  Notify user of correct arguments


	### Argument 1: project folder in pwd
	if not os.path.isdir(pwd + '/' + args[0]): 					#	Check for existance of folder containing the inputs
		print 'Project directry not present in current folder.'
		sys.exit()
	project_folder_name = args[0] 								#	set project folder name


	### Argument 2: dslop file to be used
	if not os.path.isfile(slop_db_p  + args[1]):				#	check for dslop file in slop working directory
		print 'Associated dslop file not present in:'
		print '		' + slop_db_p 
		sys.exit() 	
	dslop_file = args[1] 											#	set dslop file for DBCreate


	### Argument 3: CON file to be used (setting T-P conditions)
	if not os.path.isfile(pwd + '/' + project_folder_name + '/' + args[2]):	#	check for CON file in project directory
		print 'Associated CON file not present in:'
		print '		' + project_folder_name
		sys.exit() 	
	con_file = args[2] 											#	set SW file


	### Argument 4: Initial SW file.3i
	if not os.path.isfile(pwd + '/' + project_folder_name + '/' + args[3]):	#	check for sw file named
		print 'Associated SW.3i file not present in:'
		print '		' + project_folder_name
		sys.exit() 	
	initial_SW_file = args[3] 									#	set SW file

	### Argument 5: rock database to be used
	if not os.path.isfile(pwd + '/' + args[4]):					#	check for database called
		print 'Rock database file not present:'
		print '		' + args[4]
		sys.exit() 	
	db = args[4] 												#	set db



	### Grab paths of most current data0.tde and data0.tdf in data0 working db
	data0s_in_db = [f for f in os.listdir(data0_storage_p) if os.path.isfile(os.path.join(data0_storage_p, f))] 	#	list of all data0's in db
	tde_list = [f for f in data0s_in_db if 'tde' in f]			#	all data0.tde.R?? in db
	tdf_list = [f for f in data0s_in_db if 'tdf' in f] 			#	all data0.tdf.R?? in db
	
	old_data0_tde = find_highest_data0(tde_list, data0_storage_p, 'tde') 		#	determine file path and name of most current data0 release
	old_data0_tdf = find_highest_data0(tdf_list, data0_storage_p, 'tdf') 		#	determine file path and name of most current data0 release


	return project_folder_name, dslop_file, con_file, initial_SW_file, db, old_data0_tde, old_data0_tdf	#	Grab all that shit !

def call_DBCreate(l, old_file, dslop_file, con_file, project_path):
	print 'calling DBCreate on ', l
	### l = tde or tdf, old file is path + file of old data0, project path is the cirrent project folder
	# old_file is the whole old path
	os.chdir(pwd) 								#	move to pwd
	### Clean up the home folder
	mk_check_del_file('logK.grid')	
	mk_check_del_file('new_data0')	
	mk_check_del_file('old_data0')	
	mk_check_del_file('species.txt')	
	mk_check_del_file('spxNotFound.txt')
	mk_check_del_file('supcrt.log')
	mk_check_del_file('crtDBfile.log')	
	mk_check_del_file('crtRXNfile.log')	
	mk_check_del_file(con_file)	
	mk_check_del_file(dslop_file)	

	### Move needed run sepcific files
	shutil.copy(old_file, pwd + '/old_data0')					#	copy data0 out of data0 workign_db, and rename generaically 'old'data0'		
	shutil.copy(slop_db_p + dslop_file, pwd)					#	copy dslop out of slop workign_db 
	shutil.copy(project_path + con_file, pwd)					#	copy con file out of project dir

	### DBCreate fails if done in csh from here. Oddly, it works in csh manually.
	# args = ['/bin/sh', './DBCreate', '-c', con_file, '-s', 'old_data0', '-o', 'new_data0', '-d', dslop_file, '-eq36'] 											#	Arguments to call
	# with open(os.devnull, 'w') as fp: 								#	use of devnull will allow me to supress written output, and hopefull speed things up.
	# 	Popen(args, stdout=fp).wait()								#	wait should give time for process to complete 
	
	args = ['/bin/sh', './DBCreate', '-c', con_file, '-s', 'old_data0', '-o', 'new_data0', '-d', dslop_file, '-eq36'] 		#	Arguements to be run inside csh. ./runeq3 data0_suffix .3i file
	runs_script_and_wait(args, 3.0, ['logK.grid', 'new_data0', 'old_data0', 'species.txt', 'spxNotFound.txt', 'supcrt.log', 'crtDBfile.log', 'crtRXNfile.log']) 			#	call args, and list files to wait for
	sleep(2.0)


	ck_for_empty_file('new_data0')									#	check to make sure that the output in fact worked.

	### Deal with Outputs. ALl DBCreate info needed is stored instide the prodject folder, inside its
	### its own suffix specific folder. IE, each project has its own local copies of its data0 version.
	### The origional is written on the first life of the file. This means that those copies in data0/working_db
	### are the origional masters which get altered.
	try:
		DB_path = project_path + '/local_' + l + '_data0'
		mk_check_del_directory(DB_path)								#	check for/build output directory to store aux DBCreate file, incase theya re needed for trouble shooting
		shutil.move('new_data0', DB_path + '/data0_local.' + l)			#	move output back. If this doesnt work, it implies that the output was not made, and thus that DBCreate failed.
		shutil.move('species.txt', DB_path)
		shutil.move('spxNotFound.txt', DB_path)
		### Clean up DBCreate files generated upon running.
		os.remove(dslop_file)
		os.remove('old_data0')
		os.remove(con_file)
		os.remove('crtDBfile.log')
		os.remove('crtRXNfile.log')
		os.remove('logk.grid')
		os.remove('supcrt.log')
	## If any of the outputs are not correct, notify of failure.
	except:
		print 'DBCreate failed to produce output on' + old_file
		sys.exit()

def call_EQPT(l, project_path):
	print 'calling eqpt on ', l
	DB_path = project_path + '/local_' + l + '_data0' 				#	path to speciifc project local data0 folder
	shutil.copy(DB_path + '/data0_local.' + l, eq_db_p +'/data0.' + l) 	#	Move local data0 to eq36 eqpt db folder and generically rename 'data0' for eqpt purposes


	###	run eqpt
	os.chdir(eq_home + '/db') 										#	eqpt must be called from the db folder
	args = ['/bin/csh', eqs_exe_p + '/runeqpt', l] 					#	Arguements to be run inside csh. ./runeqpt 'suffix'
	runs_script_and_wait(args, 2.0, ['data1','data1f']) 					#	call args, and list files to wait for

	try:															#	code still in /db
		### move / rename output
		shutil.move(eq_db_p +'/data1', eq_db_p +'/data1.' + l) 		#	rename data1 with suffix so eq3/6 can use it, stroe locallcaly in /db
		mk_check_del_file(DB_path + '/data1.' + l)					#	delet faile in copy destination if it already exists.
		shutil.copy(eq_db_p +'/data1.' + l, DB_path) 				#	copy renamed data1 back to project folder 
		mk_check_del_file(DB_path + '/data1f')						#	delet faile in copy destination if it already exists.
		shutil.move(eq_db_p +'/data1f', DB_path) 					#	move data1f to project folder.
		### Clean up eq_db_p
		os.remove(eq_db_p +'/slist')
		os.remove(eq_db_p +'/output')

	# ## If any of the outputs are not correct, notify of failure.
	except:
		print 'EQPT failed to produce output on data0.' + l
		sys.exit()

def call_EQ3(l, path, water):
	print 'calling EQ3 on ', water, ' using ', l
	os.chdir(path) 												#	bring to local folder 
	args = ['/bin/csh', eqs_exe_p + '/runeq3', l, water] 		#	Arguements to be run inside csh. ./runeq3 data0_suffix .3i file
	runs_script_and_wait(args, 0.5, ['output','pickup']) 			#	call args, and list files to wait for

	try:
		### rename output
		shutil.move('output', water[:-1] + 'o')					#	local calling of outut should work given code execution location
		shutil.move('pickup', water[:-1] + 'p')
	
	# ## If any of the outputs are not correct, notify of failure.
	except:
		print 'EQ3 failed to produce output on ' + initial_SW
		sys.exit()

def build_6i_files(l, project_path, project_folder_name, initial_SW_file, template, db):
	print 'building rxn.6i files'
	#####  Temporary work around  ######
	### run 1_6i_build.p on version 7 formate.
	### proess all one at a time:
	###		(1) convert to version 8
	###		(2) remove pickup from bottom
	###		(3) rename aq species in reactant block, as xcif6 assumes syntax here, and enters species names incompatable with my data0's
	### 	(4) re attach new version 8 SW pickup
	
	os.chdir(pwd) 															#	run python code from its home (pwd)
	six_i_folder = project_path + 'rxn_6i' 										#	storage location for all 6i files
	mk_check_del_directory(six_i_folder)								#	chk, del old if nessecary, then build 

	### Run 1_6i_build.py on version 7 template (with placeholder pickup file). This generates all 6i files from the rock db used
	args = [sys.executable, six_i_build_code, project_folder_name, initial_SW_file, template, db] 			#	Arguements to be run inside csh. python 1_6i_build.v2.py, sys.executable is used here, as POpen does not recognize 'python' being called.
	with open(os.devnull, 'w') as fp:  																		#	use of devnull will allow me to supress written output, and hopefull speed things up.
		Popen(args, stdout=fp).wait()
	os.chdir(six_i_folder) 													#	move into populated 6i folder
	if os.listdir(six_i_folder)=='':										#	check that 6i folder was populated
		print '6i storage folder empty'
		sys.exit()	

	### Rename all 6i files. This avoids complication slater when all files in the file are called, and items other than 6i (such as DS.stroe) are picked up	
	six_i_file_names, six_i_file_lists = read_inputs('6i', six_i_folder)

	### Convert all version 7 6i's just built into version 8 with xcif6
	for file in six_i_file_names: 		 									#	cycle through all version 7 6i files in 6i storage
		args = ['/bin/csh', eqs_exe_p +'/xcif6', '7.0', '8.0', 'W', file]	#	call xcif6 to convert to verion 8
		# with open(os.devnull, 'w') as fp:
		# 	Popen(args, stdout=fp).wait()
		runs_script_and_wait(args, 0.01, [file, 'newin', 'ixcon']) 				#	call args, and list files to wait for
		os.remove(file)														#	remove old version 7 file locally
		shutil.move('newin', file) 											#	rename the new version 8 output (newin) with the origional 6i name
	os.remove('ixcon')														#	remove the ixcon file generated and left whereever xcif is run


	### read sw.3p in memory for use in replacing standin pickup file (this is outside the loop as it only needs to happen once)
	sw_3p =  project_path + '/' + initial_SW_file[:-1] + 'p' 
	p = open(sw_3p, "r")
	p_lines = p.readlines()													#	read the sw.3p files lines into memory
	p.close()


	### find and remove standin pickup file from all 6i's
	os.chdir(six_i_folder)
	for file in six_i_file_names: 											#	open file one at a time
		f = open(file, "r")
		lines = f.readlines()												#	read the fiels lines into memory
		f.close()

		### find index of beginning of stand-in pickup file within 6i.
		x = len(lines) - 1													#	begin at end of file			
		while not re.findall('^\*\-{20}', lines[x]):						#	find last occurance of this string, x +1 will not be the index of which everything after is deleted.
			x -= 1

		### remove all stand-in pickup file lines
		f = open(file,'w+')
		f.writelines([item for item in lines[:(x+1)]])						#	copy all lines up to the old puckup file
		f.close 										
	

	### Rename aq species names autofilled by xcif6 with the correct data0 syntax.
	### in the space immediately below the title '^\* Reaction$' xicf6 adds aq blocks corresponding to each species inthe solid reactant.
	### It also tacks on H2O, O2(g), and H+. These can be ignored as my syntax is the same.
	### (1) Grab list of elemental comp components
	rnt_ele_list = []															#	list to store reactant oxide element names
	rnt_basis_list = []															#	list of associated basis species
	rnt_old_syntax = []															#	list of old syntax aq species to be replaced
	g = open(six_i_file_names[0],'r')											#	grab info from the first file in the list. They are all the saem becasue they were generatd formt eh same db
	g_lines = g.readlines()
	g.close()
	x = 0
	while not re.findall('^\* Elemental composition$', g_lines[x]):				#	search for beginning of oxide compoents list
		x += 1
	x += 1 																		#	move index to first element
	while not re.findall('   endit\.$', g_lines[x]):							#	untill you reach the end of the oxide descriptoin
		if not re.findall('   O', g_lines[x]):									#	exclude O
			a = str(g_lines[x]).split()[0]
			rnt_ele_list.append(a)
			x += 1
		else:
			x += 1 																#	to skip the '   O'
	### build list of old syntax aq species to be replaced
	x += 3 																		#	this skips to the reaction block with the old syntax to be replaced. It assumes only one reactant.
	while not re.findall('   H\+', g_lines[x]):									#	search untill the end of the regular aq components of the reaction block
			a = str(g_lines[x]).split()[0]
			rnt_old_syntax.append(a)
			x += 1
	### (2) Search for the syntax of assciated basis species in data0 being used.
	h = open(project_path + '/local_' + l + '_data0/data0_local.' + l) 			# search local projects data0.tde file 
	h_lines = h.readlines()
	h.close()
	x = 0
	while not re.findall('^basis species', h_lines[x]): 						# 	brings you to the top of basis species list
		x += 1
	for ele in rnt_ele_list:													#	cycle through all elements
		y = x 																	#	reset index to top of basis
		while not re.findall('^auxiliary basis species', h_lines[y]): 			#	stop looking if aux is reached
			if re.findall(' ' + ele, h_lines[y]):								#	find the line with the elemetn of interest
				while not re.findall('^\+\-{20}', h_lines[y]):					#	back up to being of species block
					y -= 1
				y += 1
				a = str(h_lines[y]).split()										
				rnt_basis_list.append(a[0])										#	grab species name
				break 															# 	go to next element
			else:
				y += 1
	line_len = 27 																#	total potential spces used in species names in pickup file
	old_syntax = []																#	build list for old syntax being replaced with correct spaces
	new_syntax = [] 															#	build list for new syntax to replace the old, with correct spaces
	for item in rnt_old_syntax:													#	build the list
		old_syntax.append('   ' + str(item) + ' '*(line_len - len(item)))
	for item in rnt_basis_list:													#	build the list
		new_syntax.append('   ' + str(item) + ' '*(line_len - len(item)))
	replacements = dict(zip(old_syntax, new_syntax)) 							#	weave together for input into file
	### (3) Cycle through all files replace them.
	z = 0
	while z < len(six_i_file_names):
		inn = open(six_i_file_names[z]).read() 										#	open a read copy
		out = open(six_i_file_names[z], 'w')											#	open a write copy
		for i in replacements.keys():											# 	cycle through replace emnts
			inn = inn.replace(i, replacements[i])
		out.write(inn)
		out.close
		z += 1
	

	### Seperated the following job out of the previous for file loop, as it was not behaving properly
	### Attach proper pickup file (sw.3p) to end of 6i files
	for file in six_i_file_names: 								#	open file one at a time
		with open(file, 'a') as build:
			for p_line in p_lines: 										#	write 3p file line by line onto the end of the template
				build.write(p_line)
		build.close()

	return six_i_file_names

def six_i_syntax_fix_1(project_path, six_i_file_names):
	###	This patch is currently in place untill i can dtermine a better fix (or the real root of the problem:)
	### my eq3 generated pickup files contains no line:
	###     qgexsh=        F
	###
	###  imediately following the line:
	###
	### * Ion exchangers

	six_i_folder = project_path + 'rxn_6i' 										#	storage location for all 6i files
	os.chdir(six_i_folder) 													#	run python code from older with 6i files

	old_syntax = '* Ion exchangers' 										#	insert folloowing   '* Ion exchanger'
	new_syntax = '* Ion exchangers\n    qgexsh=        F' 					#	with insert

								
	z = 0
	while z < len(six_i_file_names): 												#	cycle through all 6i files
		inn = open(six_i_file_names[z]).read() 									#	open a read copy
		out = open(six_i_file_names[z], 'w')										#	open a write copy
		inn = inn.replace(old_syntax, new_syntax)
		out.write(inn)
		out.close
		z += 1

def six_i_syntax_fix_2(project_path, six_i_file_names):
	###	This patch is currently in place untill i can dtermine a better fix (or the real root of the problem:)
	###     chang jtemp=   2 to 1

	six_i_folder = project_path + 'rxn_6i' 										#	storage location for all 6i files
	os.chdir(six_i_folder) 													#	run python code from older with 6i files

	old_syntax = '     jtemp=  2' 	 										#	insert folloowing   '* Ion exchanger'
	new_syntax = '     jtemp=  1' 											#	with insert
	
	z = 0
	while z < len(six_i_file_names): 										#	cycle through all 6i files
		inn = open(six_i_file_names[z]).read() 								#	open a read copy
		out = open(six_i_file_names[z], 'w')								#	open a write copy
		inn = inn.replace(old_syntax, new_syntax)
		out.write(inn)
		out.close
		z += 1

def six_i_syntax_fix_3(project_path, six_i_file_names):
	###	This patch is currently in place untill i can dtermine a better fix (or the real root of the problem:)
	###     chang jtemp=   2 to 1

	six_i_folder = project_path + 'rxn_6i' 										#	storage location for all 6i files
	os.chdir(six_i_folder) 													#	run python code from older with 6i files

	old_syntax = '  iopt1-10=     1    1    0    1    0    1    0    0    0    0' 	 										#	insert folloowing   '* Ion exchanger'
	new_syntax = '  iopt1-10=     1    0    0    1    1    0    0    0    0    0' 											#	with insert
	
	z = 0
	while z < len(six_i_file_names): 										#	cycle through all 6i files
		inn = open(six_i_file_names[z]).read() 								#	open a read copy
		out = open(six_i_file_names[z], 'w')								#	open a write copy
		inn = inn.replace(old_syntax, new_syntax)
		out.write(inn)
		out.close
		z += 1

def six_i_syntax_fix_4(project_path, six_i_file_names):
	###	This patch is currently in place untill i can dtermine a better fix (or the real root of the problem:)
	###     chang jtemp=   2 to 1

	six_i_folder = project_path + 'rxn_6i' 										#	storage location for all 6i files
	os.chdir(six_i_folder) 													#	run python code from older with 6i files

	old_syntax = ' iopt11-20=     0    0    0    0    0   -1    0    0    0    0' 	 										#	insert folloowing   '* Ion exchanger'
	new_syntax = ' iopt11-20=     0    0    0    0    0   -1    0   -1    0    0' 											#	with insert
	
	z = 0
	while z < len(six_i_file_names): 										#	cycle through all 6i files
		inn = open(six_i_file_names[z]).read() 								#	open a read copy
		out = open(six_i_file_names[z], 'w')								#	open a write copy
		inn = inn.replace(old_syntax, new_syntax)
		out.write(inn)
		out.close
		z += 1

def six_i_syntax_fix_5(project_path, six_i_file_names):
	###	This patch is currently in place untill i can dtermine a better fix. it removes Ti from reactant
	six_i_folder = project_path + 'rxn_6i' 										#	storage location for all 6i files
	os.chdir(six_i_folder) 														#	run python code from older with 6i files

	for file in six_i_file_names:												#	read in file
		f = open(file, "r")
		lines = f.readlines()
		f.close()

		f = open(file, "w")														#	write all parts of file accept the foolowing two lines
		for line in lines:
			if not re.findall('^   Ti', line):
				if not re.findall('^   Ti\+4', line):
					f.write(line)

def six_i_syntax_fix_6(project_path, six_i_file_names):
	###	This patch is currently in place untill i can dtermine a better fix. it removes Ti from reactant
	six_i_folder = project_path + 'rxn_6i' 										#	storage location for all 6i files
	os.chdir(six_i_folder) 														#	run python code from older with 6i files

	for file in six_i_file_names:												#	read in file
		f = open(file, "r")
		lines = f.readlines()
		f.close()

		f = open(file, "w")														#	write all parts of file accept the foolowing two lines
		for line in lines:
			if not re.findall('^   P', line):
				if not re.findall('^   HPO4-2', line):
					f.write(line)

def six_i_syntax_fix_7(project_path, six_i_file_names):
	###	This patch is currently in place untill i can dtermine a better fix. it removes Ti from reactant
	six_i_folder = project_path + 'rxn_6i' 										#	storage location for all 6i files
	os.chdir(six_i_folder) 														#	run python code from older with 6i files

	for file in six_i_file_names:												#	read in file
		f = open(file, "r")
		lines = f.readlines()
		f.close()

		f = open(file, "w")														#	write all parts of file accept the foolowing two lines
		for line in lines:
			if not re.findall('^   Mn', line):
				if not re.findall('^   Mn\+2', line):
					f.write(line)


def run_6i_rxn(l, project_path, six_i_file_names):
	### Old way of calling eq6. abandonded because it does not wait for output.
	# os.environ[eq_env]
	# args = ['/bin/csh', eqs_exe_p + '/runeq6', l, 'AARR1_j43.6i'] #	Arguements to be run inside csh. ./runeq3 data0_suffix .3i file
	# with open(os.devnull, 'w') as fp:  							#	use of devnull will allow me to supress written output, and hopefull speed things up.
	# 	p = Popen(args, stdout=fp)


	os.chdir(project_path) 												#	move to project folder
	mk_check_del_directory('rxn_6o')									#	build output directory
	mk_check_del_directory('rxn_6p')									#	build pickup directory

	six_i_folder = project_path + 'rxn_6i' 	
	os.chdir(six_i_folder)												#	move into 6i folder
	

	for file in six_i_file_names:										#	cycle through all 6i files
		print 'running ', file
		args = ['/bin/csh', eqs_exe_p + '/runeq6', l, file] 			#	Arguements to be run inside csh. ./runeq3 data0_suffix .3i file
		runs_script_and_wait(args, 4.0, ['output','pickup']) 				#	call args, and list files to wait for
		try:															#	try is entered here, becuase if EQ6 fails, it will produce a partial 6o, but no 6p
			shutil.move('output', project_path + 'rxn_6o/' + file[:-1] + 'o') 	#	rename and move output
			shutil.move('pickup', project_path + 'rxn_6p/' + file[:-1] + 'p')	#	rename and move pickup
		except:
			continue 													#	continue even if previous run failed
	os.remove('data1')

def check_xi_completion(project_path):
	print 'Checking Xi Completion'
	### Chec for xi completion on all files, and move if nessecary.
	os.chdir(project_path)
	mk_check_del_directory('rxn_xi_incomplete') 					# 	Xi incomplete storage directory

	six_o_folder = project_path + 'rxn_6o'
	six_p_folder = project_path + 'rxn_6p'							
	os.chdir(six_o_folder) 										#	move to output folder.
	six_o_file_names, six_o_file_lists = read_inputs('6o', six_o_folder)

	for file in six_o_file_names:
		six_o = open(file, "r") 								#	Read current file
		six_o_lines = six_o.readlines() 						#	stored 6p file lines
		six_o.close()
	
		x = 0
		while x < len(six_o_lines):			
			if re.findall('^     xistti=', six_o_lines[x]):	 	#	beginning of line with ximax contains xistti
				a = str(six_o_lines[x])
				b = a.split()
				xi_max = float(b[3]) 							#	the 4 item on the xistti / ximaxi line is the value for ximax
				break
			else:
				x += 1
		### then search from the bottom up unill 
		### check xi completion
		x = len(six_o_lines) -1
		while x > 0: 							#	count down from end of file
			if re.findall('^ {20}Xi=', six_o_lines[x]):
				a = str(six_o_lines[x])
				b = a.split()		
				xi_reached = float(b[1])	
				break 							#	this should leave x where i want it, on the line number of the last reaction step
			else:
				x -= 1		
		### Exit if xi_max is not reached
		if xi_reached != xi_max:
			shutil.move(file, project_path + 'rxn_xi_incomplete') 									#	move out 6o
			try:
				shutil.move(six_p_folder + '/' + file[:-1] + 'p', project_path + 'rxn_xi_incomplete') 	#	move 6p as well
			except:
				continue
def determine_redox_ele_and_sp(l, project_path):
	### Load species list from e_inhibited data0 file	
	
	os.chdir(project_path + 'local_' + l + '_data0') 			#	move to e-separated data0 filder
	tdf = open('data0_local.' + l, "r") 						#	open e-separated data0 for reading into lines
	tdf_lines = tdf.readlines()	 								#	READ THAT SHIT !
	tdf.close() 												#	All Done


	### Build list of redox active elements (redox_active_ele)
	redox_active_ele = []
	x = 0
	while not re.findall('^elements', tdf_lines[x]): 			#	search for top of elements section.
		x += 1
	x += 2 														#	found top of element, move to first one
	while not re.findall('^basis species', tdf_lines[x]):		#	stop if basis species are encountered
		a = str(tdf_lines[x])									
		b = a.split()
		if re.findall('[0-9]', b[0]): 							#	look for numbers in the element name
			redox_active_ele.append(b[0])
			x += 1
		else:
			x += 1
	### Should now be at the top of the Basis Species

	### build list of all basis species in data0, with associated element
	basis_ele = [] 														#	build_list for basis species in data0, [species, element component, sto]
	x += 1 																#	move to first basis (h2o) 
	while x < len(tdf_lines):
		if re.findall('\+-----', tdf_lines[x]): 						#	This should find the fist as well, as the x indexed begins on top of the +---- above h2o
			x += 1
			a = str(tdf_lines[x])
			b = a.split()
			spe = b[0] 													#	basis species name of this block
			### Exit upon hitting the auxiliary set
			if re.findall('^auxiliary', spe):
				x = len(tdf_lines) 										#	skip to exit condition
			### deal with H2O, H+, and O2(g)
			elif re.findall('^[hH]2[oO]$', spe) or re.findall('^[hH]\+$', spe) or re.findall('^[oO]2\([gG]\)$', spe):							#	if it is h2o
				ele = 'NONE' 											#	default to an element of none
				sto = 'NONE'
				basis_ele.append([spe, ele, sto])
				x += 1 													# 	back to while loop
			### All other basis species
			else:														#	seach for elements other than H or O in remainging basis species
				x += 1	
				while not re.findall('\+-----', tdf_lines[x]): 			#	stop if next basis is reached
					if re.findall('element\(s\)', tdf_lines[x]): 		#	data of interest is ononeline bloew this line
						x += 1
						a = str(tdf_lines[x])
						b = a.split()
						c = [item for item in b[1::2] if item != 'H' and item != 'h' and item != 'O' and item != 'o'] 	#	for element that is not O or H, should only be 1
						d = re.search(c[0], a)
						e = d.start()
						sto = a[(e - 7) : (e - 1)] 						#	grab element stochiometry (dont knwo if i will end oup needing this)
						basis_ele.append([spe, c[0], sto]) 				#	append this basis species info to list
						x += 1
					else:
						x += 1
		else:
			x += 1


	### Determine what aqueous species are built from each basis species in data0
	basis_composite_list = [] 											#	this is the list of lists (aqueous species linked a specific basis species) Same index as the basis_search list.
	for item in redox_active_ele: 										#	for each unique e-separated element in basis_search set by user in edit section above
		species_list = generate_redox_species_list(tdf_lines, item) 	#	find all non-solids which contain the basis element (item)
		basis_composite_list.append(species_list)   					#	Add item species list to 


	return basis_ele, redox_active_ele, basis_composite_list,
def generate_redox_species_list(tdf_lines, item):
	###	Search the lines of data0.tdf.Rxxx
	ss = ' ' + item 									#	search string
	species_list = [] 									#	item (ele) specific list of associated species in data0

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

def generate_final_rock_stats(project_path, project_folder_name):
	print 'Adding aq and run data to rock stats'
	os.chdir(pwd)

	### the default rock db is hardwired into Search_6o.py. So are the aq species being sought. They need to be changed there.
	args = [sys.executable, 'Search_6o.py', project_folder_name, 'rxn', 'final'] 						#	Arguements to be run inside csh. python 1_6i_build.v2.py, sys.executable is used here, as POpen does not recognize 'python' being called.
	with open(os.devnull, 'w') as fp:  																	#	use of devnull will allow me to supress written output, and hopefull speed things up.
		Popen(args, stdout=fp).wait()
		sleep(3.0)

def build_vfl_3i_redox_sep_files(project_path, basis_ele, redox_active_ele, basis_composite_list):
	fail_list = [] 																#	for failed files
	os.chdir(project_path) 														#	move to project folder
	mk_check_del_directory('vfl_3i')											#	build storage for vfl.3i
	six_o_folder = project_path + 'rxn_6o'							
	os.chdir(six_o_folder) 														#	move to pickup folder.
	six_o_file_names, file_list = read_inputs('6o', six_o_folder) 				#	make new list of inputs, as all 6i's on previously used last might not have run thorugh to completion


	### make basis composite list re compatible.
	basis_composite_list_re = copy.deepcopy(basis_composite_list) 	 			#	deepcopy is needed here to keep alterations downstream from affect both lists. As 						#	copy an re compatible version
	for item in basis_composite_list_re: 										#	all [sto, sp] pairs within a given basis species
		for species in item:
			if '+' in species[1]:
				species[1] = species[1].replace('+', '\+') 


	for file in six_o_file_names:
		six_o = open(file, "r") 												#	Read current file
		six_o_lines = six_o.readlines() 										#	stored 6p file lines
		six_o.close()

		###Grab needed info from 6o
		t, pH, fO2, six_o_basis, basis_composite_list_individual = grab_6o_piecies(six_o_lines, file, basis_composite_list, basis_composite_list_re)

		### Prepare the new basis species list for .3i file from combination of .6o and the basis_search items listed by user
		new_basis_value = prep_species_list(redox_active_ele, basis_ele, basis_composite_list_individual, six_o_basis)


		### To remove all 0.0's
		# new_basis_value_lim = [item for item in new_basis_value if item[1] != str(format_e(0.0,5))]
		### Instead of removing basis speceis that have valiues of 0.0, they need to be made extremely low.
		### Seeing as the comp. of the vfl (being in the pickup position), needs to be all-incompassing of the system composistion.
		for item in new_basis_value:
			if item[1] == str(format_e(0.0,5)):
				item[1] = str(format_e(1.0E-99,5))


		### Write to new .3i file
		write_3i(project_path, file, new_basis_value, t, fO2, pH)
		os.chdir(six_o_folder) 
def grab_6o_piecies(six_o_lines, file, basis_composite_list, basis_composite_list_re):
	
	### defaults to be replaced
	t = 'NaN'
	pH = 'NaN'
	fO2 = 'NaN'
	six_o_basis = []
	basis_composite_list_individual = copy.deepcopy(basis_composite_list) 	#	create 6o_specific basis composite list, since they will all oave different molalities add

	### search from the bottom up unill last xi is encountered
	x = len(six_o_lines) -1
	while x > 0: 							#	count down from end of file
		if re.findall('^ {20}Xi=', six_o_lines[x]):
			a = str(six_o_lines[x])
			b = a.split()		
			xi_reached = float(b[1])	
			break 							#	this should leave x where i want it, on the line number of the last reaction step
		else:
			x -= 1		

	### Grab needed components from last step (x should already be at the correct index)
	while not re.findall('^ {16}--- Distribution of Aqueous Solute Species ---', six_o_lines[x]): 				#	Now count back to the end of the file
		### Grab t
		if re.findall('^ Temperature=', six_o_lines[x]):
			a = str(six_o_lines[x])
			b = a.split()		
			t = "%.4f" % float(b[1]) 							#	needs to be formatted for 4 decimal places.
			x += 1
		### Grab pH
		elif re.findall('^ NBS pH scale', six_o_lines[x]):
			a = str(six_o_lines[x])
			b = a.split()		
			pH = float(b[3])
			x += 1
		### Grab log fO2
		elif re.findall('^ {16}Log oxygen fugacity=', six_o_lines[x]):
			a = str(six_o_lines[x])
			b = a.split()		
			fO2 = float(b[3])
			x += 1
		### Grab aq element totals
		elif re.findall('^ {11}--- Elemental Composition of the Aqueous Solution ---', six_o_lines[x]):	
			x += 4 													#	skip to first element
			while not re.findall('^\n', six_o_lines[x]): 			#	untill the end of the data block
				a = str(six_o_lines[x])
				b = a.split()
				six_o_basis.append([str(b[0]),str(b[2])]) 			#	[name, molality] 
				x += 1		
		else:
			x += 1
	

	### Grab aq species values
	while not re.findall(' {16}--- Distribution of Aqueous Solute Species ---', six_o_lines[x]):	#	find top of aqueous solutes
		x += 1
	### Should now be at the top of aq species.
	x += 5
	top = x 										#	mark top of aq species
	ba_el = 0										#	index within basis_composite_list
	### Find the end of the aq sections
	while not re.findall('      --- Major Species by Contribution to Aqueous Mass Balances ---', six_o_lines[x]):
		x += 1
	end_aq_line = x
	
	### search for aq species.
	x = top																		#	put back on top of the aq stack
	while ba_el < len(basis_composite_list): 									#	iterate through basis elements in basis_search / basis_composite_list
		aq = 0 																	# 	index within basis_composite_list[ba_el]
		while aq < len(basis_composite_list[ba_el]):							#	iterate trough all aq species of a given basis element in basis_search
			x = top
			while x < end_aq_line: 												#	stop if species is not found id aq block
				if re.findall('^ ' + str(basis_composite_list_re[ba_el][aq][1]) + ' ', six_o_lines[x]): 	#	the trailing space should keep partial hits from happening.
					a = str(six_o_lines[x])
					b = a.split()
					basis_composite_list_individual[ba_el][aq].append(b[1]) 	#	1 grabs molality
					x = end_aq_line 											#	if species found, send to end of section
				else:
					x += 1
			if len(basis_composite_list_individual[ba_el][aq]) != 3:			#	then no value was found, set to 0
				basis_composite_list_individual[ba_el][aq].append('0.0')
				aq += 1
			else:
				aq += 1
		ba_el += 1

	return t, pH, fO2, six_o_basis, basis_composite_list_individual
def prep_species_list(redox_active_ele, basis_ele, basis_composite_list_individual, six_o_basis):
	### build new_basis_value for new basis sp. in redox_active_ele from 6o files aq section
	new_basis_value = []											#	list of eventual basis_sepcies names and values: same index as redox_active_ele
	ba = 0 															#	index for redox_active_ele
	while ba < len(redox_active_ele):
		summ = 0.0 													#	float for summing together the basis composite species
		s = 0 														#	index for composite species within a given basis element
		while s < len(basis_composite_list_individual[ba]):
			try:													#	this try statement was added to deal with syntax errors in 6o files, in which the scientific E in a number is consumed if the exponent has 3 digits.
				summ = summ + basis_composite_list_individual[ba][s][0]*float(basis_composite_list_individual[ba][s][2])		#	sum togeher of the of the species which contain the redox_active_ele element, * the sto of the basis element in them
			except:
				pass
			s += 1			
	
		summ = format_e(summ, 5)		
		new_basis_value.append(summ) 								#	add tot he basis_elemet the total molol conc. of its summed composite species
		ba += 1
	

	### Conbine redox_active_ele elements with their derived values
	s = 0
	while s < len(redox_active_ele):
		new_basis_value[s] = [redox_active_ele[s], new_basis_value[s]]
		s += 1 


	###	find the ele underlying the redox active elements, as listed in the original 6o output
	original_ele_to_redox_active = []										#	list of all elements from origional data0 which were seperated into mutiple oxidation states
	for ele in redox_active_ele:											#	for all redox active elements in list
		### srip numbers out and everything poast the number
		a = re.search('[0-9]', ele)											#	find number in all redox active elements
		b = a.start() 														#	find where the numbers start
		ele = ele[:b]														#	grab everything from beginning to there				
		original_ele_to_redox_active.append(ele)
	original_ele_to_redox_active = list(set(original_ele_to_redox_active)) 	#	delete duplicates


	### build new elements list from redox_active_ele elemetns, and remaining new_basis_values form 6o (Excluding H and O, and those representative of formly redox inactive elements). 
	s = 0
	while s < len(six_o_basis):
		if six_o_basis[s][0] not in original_ele_to_redox_active and six_o_basis[s][0] != 'O' and six_o_basis[s][0] != 'H':  		#	start with basis_searhc elements, then add all item in 6o element, if they are not already present in the basis search
			new_basis_value.append(six_o_basis[s])
			s += 1
		else:
			s += 1


	### replace the elemtns of new_bais_values with the correct associated basis_species, for input in .3i file
	r = 0
	while r < len(new_basis_value):
		s = 0
		while s < len(basis_ele):
			if new_basis_value[r][0] == basis_ele[s][1]:
				new_basis_value[r][0] = basis_ele[s][0]
				break
			else:
				s += 1
		r += 1
	return new_basis_value
def write_3i(project_path, file, new_basis_value_lim, t, fO2, pH):
	os.chdir(project_path + 'vfl_3i/')  						#	move to vfl_3i folder
	template = 'vfl_template.3i'
	new_file = file[:-3] + '_vfl.3i' 							#	new file title
	shutil.copyfile(project_path + template, new_file)			# 	Copy template and put it vfl_3i folder							#	copy 3i template for editing
	title ='EQ3NR input file name= ' 							#	line 1 old header
	newtitle = title + file[:-3] 								#	line 1 new header
	fO2 = format_e(fO2, 5) 										#	convert to scientific notation
	t = format_e(float(t), 5) 									#	convert to scientific notation
	for line in fileinput.input(new_file, inplace=1):			#	lines to replace
		line = re.sub(title, newtitle, line.rstrip())
		line = re.sub('     tempc=', '     tempc=  ' + str(t), line.rstrip())   									#	print t:  13 spaces total to get the end of the 4 decimals places alotted to t
		line = re.sub('    fo2lgi=.{20}', '    fo2lgi= ' + str(fO2) + (19 - len(str(fO2)))*' ', line.rstrip()) 		# 	print fep: 15 spaces
		print(line)																									#	write into file
	fileinput.close()															


	### reopen and write from the end of the file
	with open(new_file, 'a') as build:
		line_1 = 'species= '
		line_2 = '   jflgi=  0    covali=  ' 								#	0 for total_molality

		for item in new_basis_value_lim:
			species_add = line_1 + item[0] + '\n' + line_2 + str(item[1]) + '\n'
			build.write(species_add) 
		pH_add = line_1 + 'H+' + '\n' + '   jflgi= 20    covali=  ' + str(pH) + '\n'		#	Add pH
		build.write(pH_add)
	
		end = ('endit.' + '\n' +
			'* Ion exchangers' + '\n' +
			'    qgexsh=        F' + '\n' +
			'       net=   0' + '\n' +
			'* Ion exchanger compositions' + '\n' +
			'      neti=   0' + '\n' +
			'* Solid solution compositions' + '\n' +
			'      nxti=   0' + '\n' +
			'* Alter/suppress options' + '\n' +
			'     nxmod=   0' + '\n' +
			'* Iopt, iopg, iopr, and iodb options' + '\n' +
			'*               1    2    3    4    5    6    7    8    9   10' + '\n' +
			'  iopt1-10=     0    0    0    0    0    0    0    0    0    0' + '\n' +
			' iopt11-20=     0    0    0    0    0    0    0    0    3    0' + '\n' +
			'  iopg1-10=     0    0    0    0    0    0    0    0    0    0' + '\n' +
			' iopg11-20=     0    0    0    0    0    0    0    0    0    0' + '\n' +
			'  iopr1-10=     0    0    0    0    0    0    0    0    0    0' + '\n' +
			' iopr11-20=     0    0    0    0    0    0    0    0    0    0' + '\n' +
			'  iodb1-10=     2    0    2    2    0    0    0    0    0    0' + '\n' +
			' iodb11-20=     0    0    0    0    0    0    0    0    0    0' + '\n' +
			'* Numerical parameters' + '\n' +
			'     tolbt=  0.00000E+00     toldl=  0.00000E+00' + '\n' +
			'    itermx=   0' + '\n' +
			'* Ordinary basis switches' + '\n' +
			'    nobswt=   0' + '\n' +
			'* Saturation flag tolerance' + '\n' +
			'    tolspf=  0.00000E+00' + '\n' +
			'* Aqueous phase scale factor' + '\n' +
			'    scamas=  1.00000E+00')
		build.write(end) 								# 	Add ending and close
	build.close()

def run_vfl_3i(project_path, l):
	os.chdir(project_path) 														#	move to project folder
	mk_check_del_directory('vfl_3p')											#	build storage for vfl.3p
	mk_check_del_directory('vfl_3o')											#	build storage for vfl.3o
	vfl_3i_folder = project_path + 'vfl_3i'
	os.chdir(vfl_3i_folder) 													#	move to project folder
	file_name, file_list = read_inputs('3i', vfl_3i_folder) 					#	make new list of inputs, as all 6i's on previously used last might not have run thorugh to completion

	for file in file_name:														#	cycle through all vfl.3i files
		args = ['/bin/csh', eqs_exe_p + '/runeq3', l, file] 					#	Arguements to be run inside csh. ./runeq3 data0_suffix .3i file
		runs_script_and_wait(args, 0.5, ['output','pickup']) 					#	call args, time to wait, and list files to wait for
		try:																	#	try here, becuase a failed run may not produce outputs or pickups
			shutil.move('output', project_path + 'vfl_3o/' + file[:-1] + 'o') 		#	rename and move output
			shutil.move('pickup', project_path + 'vfl_3p/' + file[:-1] + 'p')		#	rename and move pickup
		except:
			continue
	os.remove('data1')

def check_vfl(project_path):
	### get rid of all 6o files, that dont have an associated 6p file, as they did not succeed.
	### this way, when the next function calls the list of 6o files for 6i build, it will not include incomplete 6o's
	os.chdir(project_path)
	mk_check_del_directory('vfl_3i_incomplete') 								# 	fail storage directory


	vfl_3p = project_path + 'vfl_3p' 											#	vfl pickup directory to look in
	vfl_3o = project_path + 'vfl_3o' 											#	vfl output directory to look in
	file_name_3p, file_list = read_inputs('3p', vfl_3p)  						#	generate list of all outputs in that directory
	file_name_3o, file_list = read_inputs('3o', vfl_3o)  						#	generate list of all outputs in that directory

	stay = [file[:-1] + 'o' for file in file_name_3p]							#	convert 3p files to 3o name
	move = [file for file in file_name_3o if file not in stay] 					#	files to move to vfl_3i incomplete folder

	for file in move:
		shutil.move(vfl_3o + '/' + file, project_path + 'vfl_3i_incomplete')

def build_sw_tdf(project_path, initial_SW_file,  basis_ele, redox_active_ele, basis_composite_list):
	os.chdir(project_path) 										#	move to e-separated data0 filder
	sw_3o = initial_SW_file[:-1] + 'o'
	three_o = open(initial_SW_file[:-1] + 'o', "r") 					#	open e-separated data0 for reading into lines
	sw_3o_lines = three_o.readlines()	 						#	READ THAT SHIT !
	three_o.close() 											#	All Done


	### make basis composite list re compatible.
	basis_composite_list_re = copy.deepcopy(basis_composite_list) 	 			#	deepcopy is needed here to keep alterations downstream from affect both lists. As 						#	copy an re compatible version
	for item in basis_composite_list_re: 										#	all [sto, sp] pairs within a given basis species
		for species in item:
			if '+' in species[1]:
				species[1] = species[1].replace('+', '\+') 


	t, pH, fO2, three_o_basis, basis_composite_list_individual = grab_3o_piecies(sw_3o_lines, sw_3o, basis_composite_list, basis_composite_list_re)

	### Prepare the new basis species list for .3i file from combination of .6o and the basis_search items listed by user
	new_basis_value = prep_species_list(redox_active_ele, basis_ele, basis_composite_list_individual, three_o_basis)

	### Instead of removing basis speceis that have valiues of 0.0, they need to be made extremely low.
	### Seeing as the comp. of the vfl (being in the pickup position), needs to be all-incompassing of the system composistion.
	for item in new_basis_value:
		if item[1] == str(format_e(0.0,5)):
			item[1] = str(format_e(1.0E-99,5))

	### Write to new .3i file
	new_file = write_sw_tdf_3i(project_path, initial_SW_file, new_basis_value, t, fO2, pH)
	return new_file
def grab_3o_piecies(sw_3o_lines, file, basis_composite_list, basis_composite_list_re):
	
	### defaults to be replaced
	t = 'NaN'
	pH = 'NaN'
	fO2 = 'NaN'
	three_o_basis = []
	basis_composite_list_individual = copy.deepcopy(basis_composite_list) 	#	create 3o_specific basis composite list


	### Grab needed components from last step (x should already be at the correct index)
	x = 0
	while not re.findall('^ {16}--- Distribution of Aqueous Solute Species ---', sw_3o_lines[x]): 				#	Now count back to the end of the file
		### Grab t
		if re.findall('^ Temperature=', sw_3o_lines[x]):
			a = str(sw_3o_lines[x])
			b = a.split()		
			t = "%.4f" % float(b[1]) 							#	needs to be formatted for 4 decimal places.
			x += 1
		### Grab pH
		elif re.findall('^ NBS pH scale', sw_3o_lines[x]):
			a = str(sw_3o_lines[x])
			b = a.split()		
			pH = -float(b[3])
			x += 1
		### Grab log fO2
		elif re.findall('^ {16}Log oxygen fugacity=', sw_3o_lines[x]):
			a = str(sw_3o_lines[x])
			b = a.split()		
			fO2 = float(b[3])
			x += 1
		### Grab aq element totals
		elif re.findall('^ {11}--- Elemental Composition of the Aqueous Solution ---', sw_3o_lines[x]):	
			x += 4 													#	skip to first element
			while not re.findall('^\n', sw_3o_lines[x]): 			#	untill the end of the data block
				a = str(sw_3o_lines[x])
				b = a.split()
				three_o_basis.append([str(b[0]),str(b[4])]) 			#	[name, molality] 
				x += 1		
		else:
			x += 1
	

	### Grab aq species values
	while not re.findall(' {16}--- Distribution of Aqueous Solute Species ---', sw_3o_lines[x]):	#	find top of aqueous solutes
		x += 1
	### Should now be at the top of aq species.
	x += 5
	top = x 										#	mark top of aq species
	ba_el = 0										#	index within basis_composite_list
	### Find the end of the aq sections
	while not re.findall('      --- Major Species by Contribution to Aqueous Mass Balances ---', sw_3o_lines[x]):
		x += 1
	end_aq_line = x

	### search for aq species.
	x = top																		#	put back on top of the aq stack
	while ba_el < len(basis_composite_list): 									#	iterate through basis elements in basis_search / basis_composite_list
		aq = 0 																	# 	index within basis_composite_list[ba_el]
		while aq < len(basis_composite_list[ba_el]):							#	iterate trough all aq species of a given basis element in basis_search
			x = top
			while x < end_aq_line: 												#	stop if species is not found id aq block
				if re.findall('^ ' + str(basis_composite_list_re[ba_el][aq][1]) + ' ', sw_3o_lines[x]): 	#	the trailing space should keep partial hits from happening.
					a = str(sw_3o_lines[x])
					b = a.split()
					basis_composite_list_individual[ba_el][aq].append(b[1]) 	#	1 grabs molality
					x = end_aq_line 											#	if species found, send to end of section
				else:
					x += 1
			if len(basis_composite_list_individual[ba_el][aq]) != 3:			#	then no value was found, set to 0
				basis_composite_list_individual[ba_el][aq].append('0.0')
				aq += 1
			else:
				aq += 1
		ba_el += 1

	return t, pH, fO2, three_o_basis, basis_composite_list_individual
def write_sw_tdf_3i(project_path, file, new_basis_value_lim, t, fO2, pH):
	template = 'sw_tdf_template.3i'
	new_file = file[:-3] + '_sw_tdf.3i' 									#	new file title
	shutil.copyfile(project_path + template, new_file)					# 	Copy template and put it vfl_3i folder							#	copy 3i template for editing
	title ='EQ3NR input file name= ' 									#	line 1 old header
	newtitle = title + file[:-3] 								#	line 1 new header
	fO2 = format_e(fO2, 5)
	t = format_e(float(t), 5)
	for line in fileinput.input(new_file, inplace=1):			#	lines to replace
		line = re.sub(title, newtitle, line.rstrip())
		line = re.sub('     tempc=', '     tempc=  ' + str(t), line.rstrip())   								#	print t:  13 spaces total to get the end of the 4 decimals places alotted to t
		line = re.sub('    fo2lgi=.{20}', '    fo2lgi= ' + str(fO2) + (19 - len(str(fO2)))*' ', line.rstrip()) 		# 	print fep: 15 spaces
		print(line)															#	write into file
	fileinput.close()														#	


	### reopen and write from the end of the file
	with open(new_file, 'a') as build:
		line_1 = 'species= '
		line_2 = '   jflgi=  0    covali=  ' 								#	0 for total_molality

		for item in new_basis_value_lim:
			species_add = line_1 + item[0] + '\n' + line_2 + str(item[1]) + '\n'
			build.write(species_add) 

		pH_add = line_1 + 'H+' + '\n' + '   jflgi= 20    covali=  ' + str(-pH) + '\n'		#	Add pH
		build.write(pH_add)
	
		end = ('endit.' + '\n' +
			'* Ion exchangers' + '\n' +
			'    qgexsh=        F' + '\n' +
			'       net=   0' + '\n' +
			'* Ion exchanger compositions' + '\n' +
			'      neti=   0' + '\n' +
			'* Solid solution compositions' + '\n' +
			'      nxti=   0' + '\n' +
			'* Alter/suppress options' + '\n' +
			'     nxmod=   0' + '\n' +
			'* Iopt, iopg, iopr, and iodb options' + '\n' +
			'*               1    2    3    4    5    6    7    8    9   10' + '\n' +
			'  iopt1-10=     0    0    0    0    0    0    0    0    0    0' + '\n' +
			' iopt11-20=     0    0    0    0    0    0    0    0    3    0' + '\n' +
			'  iopg1-10=     0    0    0    0    0    0    0    0    0    0' + '\n' +
			' iopg11-20=     0    0    0    0    0    0    0    0    0    0' + '\n' +
			'  iopr1-10=     0    0    0    0    0    0    0    0    0    0' + '\n' +
			' iopr11-20=     0    0    0    0    0    0    0    0    0    0' + '\n' +
			'  iodb1-10=     2    0    2    2    0    0    0    0    0    0' + '\n' +
			' iodb11-20=     0    0    0    0    0    0    0    0    0    0' + '\n' +
			'* Numerical parameters' + '\n' +
			'     tolbt=  0.00000E+00     toldl=  0.00000E+00' + '\n' +
			'    itermx=   0' + '\n' +
			'* Ordinary basis switches' + '\n' +
			'    nobswt=   0' + '\n' +
			'* Saturation flag tolerance' + '\n' +
			'    tolspf=  0.00000E+00' + '\n' +
			'* Aqueous phase scale factor' + '\n' +
			'    scamas=  1.00000E+00')
		build.write(end) 								# 	Add ending and close



	build.close()
	return new_file

def build_mix_files(project_path, sw_tdf_3i, mix_template):
	print 'Constructing mix.6i files'
	os.chdir(project_path) 														#	move to project folder
	mk_check_del_directory('mix_6i')											#	build storage for mix.6i
	three_p_folder = project_path + 'vfl_3p'									#	3p folder holder all vfl.3p's
	file_name, file_list = read_inputs('3p', three_p_folder) 					#	make new list of inputs, as all 6i's on previously used last might not have run thorugh to completion

	### Load sw_tdf puckup
	pickup = sw_tdf_3i[:-1] + 'p'												#	pickup name
	sw_3p = open(pickup, "r") 													#	Read current file
	sw_3p_lines = sw_3p.readlines() 											#	store that shit !
	sw_3p.close()
	### seperate jsut the pickup file portion needed (this is onlt present if iopt19 was set to 1 in the sw.3i file)
	x = 0
	while not re.findall('^\*------------------', sw_3p_lines[x]): 				#	find the top of the  reactant pickup block
		x += 1
	x += 1
	start_sw = x 																	#	section of interest begins here
	while not re.findall('^\*------------------', sw_3p_lines[x]): 				#	find the top of the  reactant pickup block
		x += 1
	end_sw = x 																	#	section of interest ends here

	os.chdir(three_p_folder) 													#	move to pickup folder.
	for file in file_name: 														#	cycle through vfl.3p files
		vfl_3p = open(file, "r") 												#	Read current file
		vfl_3p_lines = vfl_3p.readlines() 										#	stored 6p file lines
		vfl_3p.close()
		x = len(vfl_3p_lines) - 1
		while not re.findall('^\*------------------', vfl_3p_lines[x]): 			#	need to sfind the first occurance from the bottom of the file
			x -= 1
		x += 1
		start_vfl = x  															#	this should be the whole end of the file.
		write_mix_6i(project_path, file, sw_3p_lines, vfl_3p_lines, start_sw, end_sw, start_vfl, mix_template)
		os.chdir(three_p_folder) 													#	move back to pickup folder (write function takes you out)
def write_mix_6i(project_path, file, sw_3p_lines, vfl_3p_lines, start_sw, end_sw, start_vfl, mix_template):
	os.chdir(project_path + 'mix_6i/')  						#	move to vfl_3i folder
	new_file = file[:-6] + 'mix.6i' 							#	new_name must contain ridge
	shutil.copyfile(project_path + mix_template, new_file)		# 	Copy template and put it vfl_3i
	title ='EQ3NR input file name= ' 							#	line 1 old header
	newtitle = title + file[:-3] 								#	line 1 new header
	for line in fileinput.input(new_file, inplace=1):			#	lines to replace
		line = re.sub(title, newtitle, line.rstrip())			#	replace title
		print(line)												#	write into file
	fileinput.close()															

	### reopen and write from the end of the file
	with open(new_file, 'a') as build:
		### Write in reactant block
		x = start_sw
		while x < end_sw:
			if not re.findall('^      morr=', sw_3p_lines[x]):
				build.write(sw_3p_lines[x])
				x += 1
			else:
				build.write('      morr=  1.00000E+10      modr=  0.00000E+00' + '\n')
				x += 1


		middle = ('*-----------------------------------------------------------------------------' + '\n' +
			'    xistti=  1.00000E-12    ximaxi=  5.00000E+01' + '\n' +
			'    tistti=  0.00000E+00    timmxi=  1.00000E+38' + '\n' +
			'    phmini= -1.00000E+38    phmaxi=  1.00000E+38' + '\n' +
			'    ehmini= -1.00000E+38    ehmaxi=  1.00000E+38' + '\n' +
			'    o2mini= -1.00000E+38    o2maxi=  1.00000E+38' + '\n' +
			'    awmini= -1.00000E+38    awmaxi=  1.00000E+38' + '\n' +
			'    kstpmx=          999' + '\n' +
			'    dlxprn=  5.00000E+00    dlxprl=  0.04000E+00' + '\n' +
			'    dltprn=  1.00000E+38    dltprl=  1.00000E+38' + '\n' +
			'    dlhprn=  1.00000E+38    dleprn=  1.00000E+38' + '\n' +
			'    dloprn=  1.00000E+38    dlaprn=  1.00000E+38' + '\n' +
			'    ksppmx=          999' + '\n' +
			'    dlxplo=  1.00000E+38    dlxpll=  1.00000E+38' + '\n' +
			'    dltplo=  1.00000E+38    dltpll=  1.00000E+38' + '\n' +
			'    dlhplo=  1.00000E+38    dleplo=  1.00000E+38' + '\n' +
			'    dloplo=  1.00000E+38    dlaplo=  1.00000E+38' + '\n' +
			'    ksplmx=        10000' + '\n' +
			'*               1    2    3    4    5    6    7    8    9   10' + '\n' +
			'  iopt1-10=     1    0    0    0    1    0    0    0    0    0' + '\n' +
			' iopt11-20=     0    0    0    0    0   -1    0   -1    0    1' + '\n' +
			'  iopr1-10=     0    0    0    0    0    0    0   -1    0    0' + '\n' +
			' iopr11-20=     0    0    0    0    0    0    0    0    0    0' + '\n' +
			'  iodb1-10=     0    0    0    0    0    0    0    0    0    0' + '\n' +
			' iodb11-20=     0    0    0    0    0    0    0    0    0    0' + '\n' +
		    '     nxopt=  1' + '\n' +
		    '    option= All    ' + '\n' +
		    '    nxopex=  0' + '\n' +
		    '      nffg=  0' + '\n' +
		    '    nordmx=   6' + '\n' +
		    '     tolbt=  0.00000E+00     toldl=  0.00000E+00' + '\n' +
		    '    itermx=   0' + '\n' +
		    '    tolxsf=  0.00000E+00' + '\n' +
		    '    tolsat=  0.00000E+00' + '\n' +
		    '    ntrymx=   0' + '\n' +
		    '    dlxmx0=  0.00000E+00' + '\n' +
		    '    dlxdmp=  0.00000E+00' + '\n' +
			'*-----------------------------------------------------------------------------' + '\n')
		build.write(middle) 
		### tack on vfl_pickup
		x = start_vfl
		while x < len(vfl_3p_lines):
			build.write(vfl_3p_lines[x])
			x += 1





	build.close()

def run_6i_mix(l, project_path, six_i_file_names):
	os.chdir(project_path) 												#	move to project folder
	mk_check_del_directory('mix_6o')									#	build output directory
	mk_check_del_directory('mix_6p')									#	build pickup directory

	six_i_folder = project_path + 'mix_6i' 	
	os.chdir(six_i_folder)												#	move into 6i folder
	

	for file in six_i_file_names:										#	cycle through all 6i files
		print 'Running ', file
		args = ['/bin/csh', eqs_exe_p + '/runeq6', l, file] 			#	Arguements to be run inside csh. ./runeq3 data0_suffix .3i file
		runs_script_and_wait(args, 3.0, ['output','pickup']) 				#	call args, and list files to wait for
		try:
			shutil.move('output', project_path + 'mix_6o/' + file[:-1] + 'o') 	#	rename and move output
			shutil.move('pickup', project_path + 'mix_6p/' + file[:-1] + 'p')	#	rename and move pickup
		except:
			continue
	os.remove('data1')

def check_mix_xi_completion(project_path):
	### Check for failed 6i runs first, then, of those remaining,check for xi completion.
	os.chdir(project_path)
	mk_check_del_directory('mix_xi_incomplete') 					# 	Xi incomplete storage directory
	six_o_folder = project_path + 'mix_6o'							#	mix.6o folder
	six_p_folder = project_path + 'mix_6p'							#	mix.6p folder
	
	### first remove all 6o files with no associaeted 6p.
	file_name_6p, file_list = read_inputs('6p', six_p_folder)  			#	generate list of all outputs in that directory
	file_name_6o, file_list = read_inputs('6o', six_o_folder)  			#	generate list of all outputs in that directory
	stay = [file[:-1] + 'o' for file in file_name_6p]				#	convert 3p files to 3o name
	move = [file for file in file_name_6o if file not in stay] 		#	files to move to vfl_3i incomplete folder
	for file in move:
		shutil.move(six_o_folder + '/' + file, project_path + 'mix_xi_incomplete')	#	move all 6o's to the incomplete folder, if they have no associated 6p
	
	### then, of the 6o files remaining, search for xi incomplete files
	six_o_file_names, file_list = read_inputs('6o', six_o_folder)  			#	new list of 6o's, now that the failures have been removed
	os.chdir(six_o_folder)
	for file in six_o_file_names:
		six_o = open(file, "r") 								#	Read current file
		six_o_lines = six_o.readlines() 						#	stored 6p file lines
		six_o.close()
		x = 0
		while x < len(six_o_lines):			
			if re.findall('^     xistti=', six_o_lines[x]):	 	#	beginning of line with ximax contains xistti
				a = str(six_o_lines[x])
				b = a.split()
				xi_max = float(b[3]) 							#	the 4 item on the xistti / ximaxi line is the value for ximax
				break
			else:
				x += 1
		### then search from the bottom up unill 
		### check xi completion
		x = len(six_o_lines) -1
		while x > 0: 							#	count down from end of file
			if re.findall('^ {20}Xi=', six_o_lines[x]):
				a = str(six_o_lines[x])
				b = a.split()		
				xi_reached = float(b[1])	
				break 							#	this should leave x where i want it, on the line number of the last reaction step
			else:
				x -= 1		
		### Exit if xi_max is not reached
		if xi_reached != xi_max:
			try:
				shutil.move(file, project_path + 'mix_xi_incomplete') 									#	move out 6o
				shutil.move(six_p_folder + '/' + file[:-1] + 'p', project_path + 'mix_xi_incomplete') 	#	move 6p as well
			except:
				continue


##############################  Stage  7  Forward  ##############################


def build_xi_t_path(project_path, six_o_file_names, t_slices):
	print 'Determining global Xi-T/P path'
	"""
	This function determines the needed step information in a list of 6o's

	demension of xi_array, t_array, ln_n_array, and p_array:
	cut down to the first 6o Xi step, and the last t-slice many
	
		xi_array[file][t-step]
				t- step  ---->
      6o     __________________________
     file   |
	  |	    |   xi values /
	  |	    |   t values / 
	  V	    |   p values / 
		    |   or ln numbers for the beginning
		    |   of xi step
		    |
            |

	"""
	os.chdir(project_path + 'mix_6o') 									#	move to the mix output folder
	f = 0
	while f < len(six_o_file_names):
		print six_o_file_names[f]
		### remove uneeded xi steps
		six_o = open(six_o_file_names[f], "r") 							#	Read current file
		six_o_lines = six_o.readlines() 								#	stored 6p file lines
		six_o.close()
		xi, ln_n, t, p = grab_xi_ln_t_p(six_o_lines)					#	grab all xi steps [xi, line number]
		if f == 0:														#	if this is the first file, use it for dimensioning
			xi_array = np.zeros((len(six_o_file_names), len(xi)))		#	array to add all info (xi, ln, or t) (rows, cols). set up on first file, as dimensions are not known before
			ln_n_array = np.zeros((len(six_o_file_names), len(ln_n)))
			t_array = np.zeros((len(six_o_file_names), len(t)))
			p_array = np.zeros((len(six_o_file_names), len(p)))
		xi_array[f] = xi 												#	add files to arrays		
		ln_n_array[f] = ln_n				
		t_array[f] = t
		p_array[f] = p
		f += 1
	grab_range = [0] + range(-t_slices,0)								#	grab first col of array, and last (t_slice) many
	xi_array = xi_array[:,grab_range]									#	reduce arrays to the above dimensions
	ln_n_array = ln_n_array[:,grab_range]								
	t_array = t_array[:,grab_range]
	p_array = p_array[:,grab_range]
	average_t = []														#	build the average p path for con file.
	average_p = []														#	build the average p path for con file.
	average_xi = []														#	build the average p path for con file.
	x = 0 																
	while x < (t_slices + 1):
		average_t.append(format_e(np.mean(t_array[:,x]), 4))			#	t_path for con file
		average_p.append(format_e(np.mean(p_array[:,x]), 4))			#	p_path for con file
		average_xi.append(format_e(np.mean(xi_array[:,x]), 4))			#	p_path for con file
		x += 1
	return average_xi, ln_n_array, average_t, average_p 
def grab_xi_ln_t_p(lines):
	xi = []
	ln_n = [] 					#	line numbers for the arrays
	t = []
	p = []
	x = 0
	while x < len(lines): 		# builds t and zi in master file,   per .6o
		if re.findall('^ {20}Xi=', lines[x]):
			a = str(lines[x])
			b = a.split()
			xi.append(float(b[1])) 					#	[xi_value, line number]
			ln_n.append(x) 							#	for later use in possibly cutting un-needed steps out of the 6o
			x += 1
		if re.findall('^ Temperature= [^\(]', lines[x]): 	#	[^\()] needed here to keep the t_equation near the top of 6o from being grabbed
			a = str(lines[x])
			b = a.split()
			t.append(float(b[1]))
			x += 1
		if re.findall('^ Pressure=', lines[x]):
			a = str(lines[x])
			b = a.split()
			p.append(float(b[1]))
			x += 1
		else:
			x += 1      
	return xi, ln_n, t, p;

def build_con_file(project_path, project_folder_name, average_t, average_p):
	print 'Building con file for supcrtgrid'
	os.chdir(project_path)
	sample_count = str(len(average_t))
	supcrt_con_file = project_folder_name + '.con'
	mk_check_del_file(project_path + supcrt_con_file) 				#	if con already exists, then delet
	with open(supcrt_con_file, 'w') as con:		
		header = (' Line 1 (free format): isat, iopt, iplot, univar, noninc' + '\n'
			' Lines i=2.. 6 (free format): oddv1(i), oddv2(i)' + '\n'
			'******************************************************************' + '\n'
			'******************************************************************' + '\n')
		con.write(header)
		sto_string = '   0   2   2   0   ' +  sample_count
		con.write(sto_string)
		con.write('\n')
		x = 0
		while x <len(average_t):
			string = '     ' + str(average_t[x]) + '    ' + str(average_p[x])
			con.write(string)
			con.write('\n')
			x += 1
	con.close()
	return supcrt_con_file

def call_supcrtgrid(project_path, dslop_file, old_file, supcrt_con_file):
	print 'Calling supcrtgrid'
	os.chdir(pwd) 													#	bring to local folder 
	meta_dir = pwd + '/5_Meta_rxn_files' 							#	directory housing the supcrt rxn files
	supcrt_out = project_path + '/supcrt_out' 						#	directory to house the supcrt rxn files
	mk_check_del_directory(supcrt_out) 								#	chk, then buidl supcrt_out directory
	rxn_file_names, rxn_file_lists = read_inputs('rxn', meta_dir) 	# 	grab mix outputs
	for file in rxn_file_names: 									#	copy rxn files to pwd
		mk_check_del_file(file)
		shutil.copy(meta_dir + '/' + file, pwd)

	shutil.copy(old_file, dbc_exe_p + '/old_data0')					#	copy data0 out of data0 workign_db, and rename generaically 'old'data0'		
	shutil.copy(project_path + supcrt_con_file, pwd)				#	copy supcrt_con into local directory		
	shutil.copy(slop_db_p + dslop_file, dbc_exe_p)					#	copy dslop out of slop workign_db 

	args = ['/bin/sh', './Call_supcrtgrid', '-c', supcrt_con_file, '-s', 'old_data0', '-o', 'new_data0', '-d', dslop_file, '-eq36'] 											#	Arguments to call		
	with open(os.devnull, 'w') as fp: 								#	use of devnull will allow me to supress written output, and hopefull speed things up.
		Popen(args, stdout=fp).wait()								#	wait should give time for process to complete 
		sleep(1.0) 													#	give supcrtgrid a little time to finish
	for file in rxn_file_names: 									#	copy rxn files to pwd
		try:
			shutil.move(file[:-3] + 'tab', project_path + 'supcrt_out')
			# ### rename output
			# for item in out_list:
			# 	shutil.move(item, supcrt_out)					#	local calling of outut should work given code execution location
			# ## If any of the outputs are not correct, notify of failure.
		except:
			print 'supcrtgrid failed. . . fucking typical.'
			sys.exit()

	### remove temporarily needed files
	os.remove('old_data0')
	os.remove(supcrt_con_file)
	os.remove(dslop_file)
	os.remove('supcrt.log')

def load_supcrt_rxn_and_k_data(project_path):
	print 'Loading metabolic reaction LogK values from supcrtgrid output'
	"""
	This function cycles through the supcrt outputs of all metabolic reactions:
	it builds:

		rxn_list['rxn number']: A series where:
							 		 index = values
								rxn_number = [[sto1, sp1], [sto2, sp2], ...]

					for all reaction components (sp = species), with correct signs on sto (stochiometry).


		rxn_k_value_list['520']: A series where:
									index = values
								rxn_number = [k-value(t1), k-value(t2), . . . ]

					therefore rxn_k_value_list['520'][4] = the log k value for reactio 520, 
					at the 4th t-step of those grabbed. (xi = 1, plus the last t-slice many)

		rxn_ID_list: is jsut a list containing all reaction ID numbers as strings

		rxn_file_names: is just a list of the supcrt output file names containing all fo the reaction data.

	"""
	supcrt_out = project_path + 'supcrt_out' 							#	directory housing the supcrt rxn files
	os.chdir(supcrt_out)												#	move into rxn file directory
	rxn_file_names, rxn_file_lists = read_inputs('tab', supcrt_out) 	# 	load all reaction files
	### defining the total rxn list out here should allow for all rxn_files to be joined together
	### however, they are now dependent on the order in which they are read in, for total reaction order, hmm./
	### This shouldnt matter, because I will use the same file order, and thus overall rxn order, to look for Q values
	rxn_ID_list = []
	rxn_list = pd.Series([]) 	
	rxn_k_value_list = pd.Series([])
	### loop through rxn files, grabbing redox reaction stoichiometry and species info
	w = 0											#	rxn file index
	while w < len(rxn_file_names):
		f = open(rxn_file_names[w], "r") 				# 	Read current file
		lines = f.readlines()
		f.close()
		x = 0
		while x <len(lines):												#	Find all reactions in order in current file
			if re.findall('^ REACTION STOICHIOMETRY:', lines[x]):			#	while inside a given reaction	
				rxn_ID = str(lines[x-3]).split()[0]
				rxn_ID_list.append(rxn_ID)
				temp_rxn_list = []											#	temporarily stores rxn specific stoichi and species
				x += 4
				while not re.findall('^\n', lines[x]):
					a = str(lines[x])
					b = a.split()
					temp_rxn_list.append([b[0],b[1]]) 						#	tuple [stoich, species]
					x += 1
				temp_rxn_list = np.array(temp_rxn_list)
				rxn_list = rxn_list.append(pd.Series([temp_rxn_list], index=[rxn_ID]))
				# rxn_list.append(temp_rxn_list) 								#	Add local reaction to global rxn list
				x += 1

			### Find the spample specific K-values following the reaction tallied above.
			### given that x is counting throught he lines of the filel, this will keep the index of
			### all reactions and the associated K values the same.
			if re.findall('^ STANDARD STATE PROPERTIES OF THE REACTION AT ELEVATED TEMPERATURES AND PRESSURES', lines[x]):
				temp_rxn_k_value_list = []
				x += 7
				while x < len(lines) and not re.findall('^\n', lines[x]): 		#	x < len()list is here to deal with operstepping the end of a rxn file
					a = str(lines[x])
					b = a.split()
					temp_rxn_k_value_list.append(b[3]) 							#	tuple [stoich, species]
					x += 2 														#	2 here, becuase the grid of K data skips lines
				# rxn_k_value_list.append(temp_rxn_k_value_list) 				#	Add local reaction to global rxn list
				rxn_k_value_list = rxn_k_value_list.append(pd.Series([temp_rxn_k_value_list], index=[rxn_ID]))
				### There is no x +=1 sitting here, becuase it casues an over-reach at the end of a given rxn file
			else:
				x += 1    		
		w += 1
	return rxn_file_names, rxn_list, rxn_k_value_list, rxn_ID_list

def grab_e_transfered(e_list_file):
	os.chdir(pwd)
	f = open(e_list_file, "r") 				# 	Read current file
	lines = f.readlines()
	f.close()
	e_list = []
	for line in lines:
		e_list.append(float(str(line).split()[0]))
	return e_list

def check_aq_species(rxn_list):
	"""
	This function grabs all unique species from the list of all metabolic reations being considered.

	"""
	total_species = [] 										#	total species in metabolic reactions
	for reaction in rxn_list:
		for species in reaction:
			total_species.append(species[1]) 
	total_unique_species = list(set(total_species))
	unique_nonsolids = ['H2O']	
	add_items = [item for item in total_unique_species if re.findall(',AQ',item) or re.findall('[-\+]', item)]
	add_items.sort() 										#	item order is H2O, then the rest sorted alphabetically
	unique_nonsolids = unique_nonsolids + add_items 		# 	AQ species List
	unique_solids = [item for item in total_unique_species if item not in unique_nonsolids]
	unique_solids.sort()
	return total_unique_species, unique_nonsolids, unique_solids 

def grab_rxn_species_from_6o(project_path, six_o_file_names, unique_nonsolids, ln_n_array, average_t):
	print "Grabbing all needed speceis from all 6o files"
	""" 
		This function populates the below array, with all aqueous species log activities needed to compute
		the metabolisms identified in meta_rxn_list

		six_o_species_array: a three demnsional array. First index cycles through 6o files.
		for every six_o_species_array[w] = file
		exists the sub 2D array:
			six_o_species_array[w][t-step][species] = log activity of that species

			species of interest  ---->
      t      __________________________
    steps   |
	  |	    |   array of log activities
	  |	    |   for each species at
	  V	    |   each time point identified
		    |   as the first 6o Xi step, and
		    |   the last t-slice many
		    |
            |

	"""
	os.chdir(project_path + 'mix_6o')														#	step into the mix.6o folder to cycle through outputs
	species_activity_array = np.zeros((len(six_o_file_names),len(average_t), len(unique_nonsolids))) 		#	array with all species values for all files at all grabbed t-steps. The order here corresponds to the indexing order.
	species_conc_array = np.zeros((len(six_o_file_names),len(average_t), len(unique_nonsolids) - 1)) 		#	array with all species values for all files at all grabbed t-steps. The n-1 is to shorten this list to account for the lack of H2O.


	### make unique non_solids list re compatible.
	unique_nonsolids_re = copy.deepcopy(unique_nonsolids) 	 			#	deepcopy is needed here to keep alterations downstream from affect both lists. As 						#	copy an re compatible version
	
	i = 0
	while i < len(unique_nonsolids_re): 								#	make re compatible.
		if '+' in unique_nonsolids_re[i]:
			unique_nonsolids_re[i] = unique_nonsolids_re[i].replace('+', '\+')
		i += 1	


	w = 0
	while w < len(six_o_file_names):					#	cycle through files w
		six_o = open(six_o_file_names[w], "r") 			#	Read current file
		lines = six_o.readlines() 						#	stored 6o file lines
		six_o.close()
		t = 0
		while t < len(average_t): 						#	cycle through t-steps t
			x = int(ln_n_array[w][t])					#	set line index to beginning of step t		
			while not re.findall(' {14}Log activity of water=', lines[x]):
				x += 1
			a = str(lines[x])
			b = a.split()
			species_activity_array[w][t][0] = float(b[4]) 	#	index 0 here is for H2O in unique_nonsolids
			x += 1		
			while not re.findall(' {16}--- Distribution of Aqueous Solute Species ---', lines[x]):					#	find top of aq distro list
				x += 1
			x += 4 
			start = x
			while not re.findall(' {6}--- Major Species by Contribution to Aqueous Mass Balances ---', lines[x]):	#	find bottom of aq distro list
				x += 1
			end = x
			n = 1 										# 	begin at 1, as H2O has already been filled.
			###	grab all unique_non_solids
			while n < len(unique_nonsolids):
				x = start 								#	put line index back at begging of aq distro.
				while x < end:
					if re.findall('^ ' + unique_nonsolids_re[n] + ' ', lines[x]):	#	the sapce after the name should insure no partial fits at teh beginning of the string
						a = str(lines[x])
						b = a.split()
						species_activity_array[w][t][n] = float(b[4])				#	b[4] is log act
						species_conc_array[w][t][n - 1] = float(b[1])				#	b[1] is conc. n - 1 here, to offset the fact that no water is in this list as is in position 0 in the activities list.
						break
					else:
						x += 1
				n += 1
			t += 1
		w += 1

	return species_activity_array, species_conc_array

def determine_Q_values(project_path, six_o_file_names, average_t, species_activity_array, rxn_list, rxn_ID_list, unique_nonsolids, unique_solids):
	print 'Calculating logQ'
	Q_value_array = np.empty((len(six_o_file_names),len(average_t), len(rxn_ID_list)), dtype=float) 		#	empty, to eventually house the limiting concentration floats.
	""" 
		This function populates the below array, with all rxn Q values needed to compute
		metabolic return energies.

		six_o_species_array: a three demnsional array. First index cycles through 6o files.
		for every six_o_species_array[w] = file
		exists the sub 2D array:
			six_o_species_array[w][t-step][rxn] = Q value for rxn

			rxn of interest  ---->
      t      __________________________
    steps   |
	  |	    |   array of limting species
	  |	    |   identities
	  V	    |   
		    |   
		    |   
		    |

	"""
	w = 0
	while w < len(six_o_file_names):						#	cycle through files
		print '. . . ', six_o_file_names[w]
		t = 0
		while t < len(average_t): 							#	cycle through t-steps					
			r = 0
			while r < len(rxn_ID_list):						#	cycle through reactions
				build = float(0.0) 							#	Q build: begin with 0, as these are log activities beign added or subtracted together, not multiplied, (as with raw activities)
				c = 0
				while c < len(rxn_list[rxn_ID_list[r]]):	#	cycle through rxn components
					sto = rxn_list[rxn_ID_list[r]][c][0]	#	stochiometry
					spe = rxn_list[rxn_ID_list[r]][c][1]	#	species
					if spe in unique_solids: 				#	skip solids, as activity is 1
						c += 1
					elif spe in unique_nonsolids:
						ind = unique_nonsolids.index(spe)
						log_act = species_activity_array[w][t][ind]
						build = build + float(sto)*float(log_act)				#	i think that the sign on the sto negates me needing to differentiate between products and reactants
						c += 1
					else:		
						print 'Species ' + rxn_list[rxn_ID_list[r]][c][1] + ' in rxn ' + rxn_ID_list[r]
						print 'Not found in species log activity array'
						sys.exit()
				Q_value_array[w][t][r] = build 				#	add build value to Q array
				r += 1
			t += 1
		w += 1
	return Q_value_array

def solid_conc_determination(average_xi, unique_solids):
	print "Determining the concentration of dilute solids in titrant"
	""" 
		This function produces a list of series, one series for each time-step.
		Each series index is the list of solids, the value returned is their concentration at the step.

		Therefore solid_conc[t-step]['sp_name'] = conc.
	"""
	### This works, its a list of series,  there
	starting_conc = pd.Series([float(0.0), float(1*10**(-10)), float(0.0), float(1*10**(-10)), float(1*10**(-10)), float(1*10**(-10)), float(1*10**(-10)), float(1*10**(-10)), float(0.0), float(0.0)],index=unique_solids) 			#	
	solid_conc = []											#	list of solid conc series, one for each t-step
	for xi in average_xi:
		step = float(xi)*starting_conc / (float(xi) + 1.0)	#	( (mols/kg)*sw_mass = mols)  / (total mass of fluid = (sw mass + 1 kg vent fluid)
		solid_conc.append(step)
	###replace solid_conc[0] with 0.0's, as xi does not actually =0, but 10E-12. However, there is no sw.

	solid_conc[0] = pd.Series([float(0.0)]*len(unique_solids),index=unique_solids) 	
	return solid_conc

def determine_limiting_species(project_path, six_o_file_names, average_t, rxn_ID_list, rxn_list, species_conc_array, unique_nonsolids, unique_solids, solid_conc):
	print 'Determining limiting species and their concentrations'
	limiting_array = np.empty((len(six_o_file_names),len(average_t), len(rxn_ID_list)), dtype=float) 		#	empty, to eventually house the limiting concentration floats.
	"""
		This function determines, for every time step of every file, for every reaction,
		what the limiting reactant would be, for purposes of later coverting the 
		affinities to Cal/ml of fluid.

			limiting_array: a three demnsional array. First index cycles through 6o files.
		for every six_o_species_array[w] = file
		exists the sub 2D array:
			limiting_array[w][t-step][rxn] = limiting species identity for that rxn

			rxn of interest  ---->
      t      __________________________
    steps   |
	  |	    |   array of limting species
	  |	    |   identities
	  V	    |   
		    |   
		    |   
		    |

	"""
	w = 0
	while w < len(six_o_file_names):						#	cycle through files
		print '. . . ', six_o_file_names[w]
		t = 0
		while t < len(average_t): 							#	cycle through t-steps					
			r = 0
			while r < len(rxn_ID_list):						#	cycle through reactions
				# print '\n', 'reaction = ', rxn_ID_list[r]
				temp_values = []
				c = 0 										#	reaction components
				while re.findall('^-', rxn_list[rxn_ID_list[r]][c][0]):		#	consider only reactants, not products. They contain a negative out front.
					if rxn_list[rxn_ID_list[r]][c][1] != 'H2O':				#	skip water
						if rxn_list[rxn_ID_list[r]][c][1] in unique_solids:									#	if its a solid, look to the solid_conc series
							com = solid_conc[t][rxn_list[rxn_ID_list[r]][c][1]] / -float(rxn_list[rxn_ID_list[r]][c][0]) 	# conc. value, divided by the negative of the  sto. negative here to counter the negative already present in the sto.
							temp_values.append(com)											#	add this value to temp list for 'lowest' testing
							# print rxn_list[rxn_ID_list[r]][c][1], solid_conc[t][rxn_list[rxn_ID_list[r]][c][1]], com
						else: 																#	this should catch all non-colids
							ind = unique_nonsolids.index(rxn_list[rxn_ID_list[r]][c][1])-1 					#	index of the species in unique_nonsolids, which is the same as its index in the species conc. array. -1 here shfts index as water is not in conc_array, as it is in act_array.
							com = species_conc_array[w][t][ind] / -float(rxn_list[rxn_ID_list[r]][c][0])	# conc. value, divided by the negative of the  sto. negative here to counter the negative already present in the sto.
							temp_values.append(com)			#	add this value to temp list for 'lowest' testing
							# print rxn_list[rxn_ID_list[r]][c][1], species_conc_array[w][t][ind], com
						c += 1
					else:									#	skip water
						c += 1
				temp_values.sort() 							#	order the concn. to get the lowest value
				limiting_array[w][t][r] = temp_values[0]	#	if the lowest numeric is in position 0, this should add its value (but not identity) to the array.
				r += 1
			t += 1
		w += 1
	return limiting_array

def determine_affinity_and_Energy(project_path, six_o_file_names, rxn_k_value_list, average_t, e_list, rxn_ID_list, rxn_list, limiting_array, Q_value_array):
	print 'Calculating affinity, and energy (cal/ml)'
	affinity_array = np.empty((len(six_o_file_names),len(average_t), len(rxn_ID_list)), dtype=float) 		#	empty, to eventually house the limiting concentration floats.
	energy_array = np.empty((len(six_o_file_names),len(average_t), len(rxn_ID_list)), dtype=float) 		#	empty, to eventually house the limiting concentration floats.	
	""" 
		This function populates the below arrays, with all affinities, or energies (in cal/ml)
		for all metabolic rxn's

		six_o_species_array: a three demnsional array. First index cycles through 6o files.
		for every six_o_species_array[w] = file
		exists the sub 2D array:
			six_o_species_array[w][t-step][rxn] = affinity of reaction (or energy in cal/ml)

			rxn of interest  ---->
      t      __________________________
    steps   |
	  |	    |   array of affinities/energys
	  |	    |   
	  V	    |   
		    |   
		    |   
		    |

		affinity = R * T * ln(K/Q)    /  e_transfered 							#	cal / mol e transfered
				 = R(cal) * T(K) * 2.303*log(K/Q)   /  e_transfered
				 = R(cal) * T(K) * 2.303*( log(K) - log(Q) )  /  e_transfered



	"""
	w = 0
	while w < len(six_o_file_names):						#	cycle through files
		print '. . . ', six_o_file_names[w]
		t = 0
		while t < len(average_t): 							#	cycle through t-steps					
			r = 0
			while r < len(rxn_ID_list):								#	cycle through reactions
				
				# print 'reaction  = ', rxn_ID_list[r]
				# print rxn_list[rxn_ID_list[r]]
				lim_conc = limiting_array[w][t][r] 					#	limiting concn. molal. for energy determination
				logQ = float(Q_value_array[w][t][r]) 				#	reaction Q
				logK = float(rxn_k_value_list[rxn_ID_list[r]][t]) 	#	reaction log K, 
				t_k = float(average_t[t]) + kelvin					#	temperature in kelvin at this time step
				e_trans = e_list[r] 								#	
				# print 'lim conc. = ', lim_conc 
				# print 'logQ = ', logQ
				# print 'logK = ', logK
				# print 'T = ', t_k
				affinity = ( R*t_k*2.303*( logK - logQ ) )			#	affinity, in cal / mol     (mol of what?)			
				# print 'affinity = ', affinity
				affinity_e = affinity / float(e_trans)				#	affinity, in cal / mol e transfered				
				# print 'affinity per e= ', affinity_e
				affinity_array[w][t][r] = affinity_e 				#	add to affinity array
				energy_array[w][t][r] = lim_conc*affinity 			#	add to energy array. I thnk i can use affinity here, as e-tranfered would cancle out from before.
				# print 'energy = ', lim_conc*affinity 				#	divsion by 1000 here takes units from cal/kg, to ~ cal/ml (with the orugh estimate of 1kg ~ 1L)
				# print '\n'
				r += 1
			t += 1
		w += 1
	###	convert to log values
	blank = np.zeros((len(six_o_file_names),len(average_t), len(rxn_ID_list)), dtype=float)		#	zeros to use as stand in for -afinities (which I cannot convert to log values)
	affinity_array = np.where(affinity_array>0, np.log10(affinity_array), blank)				#	if array value > 0, use its log10 value, otherwise, use the value from the blank array.
	energy_array = np.where(energy_array>0, np.log10(energy_array), blank)
	return affinity_array, energy_array

def convert_rock_to_rxn_array(six_o_file_names, average_t, rxn_ID_list, rock_affinity_array, rock_energy_array):
	print 'Convert rock based rxn outputs to Rxn based outputs'
	""" 
		This function converts the existing file(rock) based arrays, to arrays based on rxns:

		IE, the first index now cycles through Rxn's, not files
		six_o_species_array[w] = rxn
		exists the sub 2D array:
			six_o_species_array[rxn][t-step][w] = affinity of reaction (or energy in cal/ml)

			host_rock of interest  ---->
      t      __________________________
    steps   |
	  |	    |   array of affinities/energys
	  |	    |   
	  V	    |   
		    |   
		    |   
		    |
                     -  Added below  -
     rock   
     stats  |
	  |	    |   array of host rock oxide components
	  |	    |   
	  V	    |   
		    |   
		    |   
	

	"""
	rxn_affinity_array = np.empty((len(rxn_ID_list),len(average_t), len(six_o_file_names)), dtype=float) 		#	empty, to eventually house the energy and affinity outputs.
	rxn_energy_array = np.empty((len(rxn_ID_list),len(average_t), len(six_o_file_names)), dtype=float) 			#	empty, to eventually house the energy and affinity outputs.	
	###  all col[0] from each array, will be stacked together, to give a new array
	r = 0
	while r < len(rxn_ID_list):
		w = 0
		while w < len(six_o_file_names):
			t = 0
			while t <len(average_t):
				rxn_affinity_array[r][t][w] = rock_affinity_array[w][t][r]
				rxn_energy_array[r][t][w] = rock_energy_array[w][t][r]
				t += 1
			w += 1
		r += 1
	return rxn_affinity_array, rxn_energy_array

def load_rock_stats(project_path, project_folder_name, six_o_file_names, rxn_ID_list, average_t, rxn_affinity_array, rxn_energy_array):
	print 'Postfix all rock stats into rxn based arrays'
	### Postfix Segment info onto all rxn_arrays.
	rock_stats = project_path + '/' + project_folder_name + '_rock_stats.csv' 			#	build on path
	os.chdir(pwd)
	rxn_array_row_names = []

	print rock_stats

	with open(rock_stats, 'rU') as csvfile:
		csvfile.seek(0)
		reader = csv.reader(csvfile)
		### Convert six_o_names to jsut the perceeding segment name.
		seg_names = []
		for f in six_o_file_names:
			seg_names.append(f[:(f.find('_'))])				#	Grab the segment name which precedes the  _  in the six_o_name


		for row in reader:
			if row[0] == '': 					#	grab what will become row-titles
				a = row[1::]
				x = 0
				while x < len(a):
					rxn_array_row_names.append(a[x])
					x += 1

	seg_info = [] 										#	file indexed list of lists(all data associated with a given file to be added to the arrays)
	w = 0
	while w < len(seg_names):							#	find data rows associated with all file names
		with open(rock_stats, 'rU') as csvfile:
			csvfile.seek(0)
			reader = csv.reader(csvfile)
			for row in reader:
				name = seg_names[w]
				if row[0] == name:
					seg_info.append(row[1::])
		w += 1

	### build array to be added wholesale onto each rxn_array
	rock_stat_array = np.empty((len(rxn_array_row_names), len(seg_names)), dtype=float)
	row = 0
	while row < len(rxn_array_row_names):
		seg = 0
		while seg < len(seg_names):
			rock_stat_array[row][seg] = seg_info[seg][row]
			seg += 1
		row += 1

	### new final build_array with correct dimensions of the final product
	combined_energy_array = np.zeros((len(rxn_ID_list), (len(average_t) + len(rxn_array_row_names)), len(six_o_file_names)), dtype=float) 			#	empty array with the dimensions needed for adding the rock stats
	combined_affinity_array = np.zeros((len(rxn_ID_list), (len(average_t) + len(rxn_array_row_names)), len(six_o_file_names)), dtype=float) 		#	empty array with the dimensions needed for adding the rock stats
	combined_index = average_t + rxn_array_row_names		#	new index, to use as row names in the final dataframe

	### populate the new combined array
	w = 0
	while w < len(rxn_ID_list):															#	For each sub array  (t, host_rock), add the rock stat array built above
		combined_energy_array[w] = np.vstack((rxn_energy_array[w], rock_stat_array))	#	stack together the origiona data array, n top of the new rock stat array
		combined_affinity_array[w] = np.vstack((rxn_energy_array[w], rock_stat_array))
		w += 1

	return combined_energy_array, combined_affinity_array, combined_index

def array_to_csv(project_path, name, array, index_1, index_2, index_3):
	print 'Writing to csv'
	### index 1 = csv file name (a separate csv is created for each index step here)
	### index 2 = t stpes
	### index 3 = reaction number
	""" 
		This function populates produces csv files for each index_1 step in a 3D array.

		Designed here to output Affinity grids and energy grads, one for each 6o file assessed.

		Inputs:
			array: 3 D array, such as energy_array and affinity_array

			array[index_1][index_2][index_3]

			index_1: are the ordered names of otehrmost index. file names
			index_2: temperature steps
			index_3: reaction numbers

			Therefore within each file 'array[w]'' exists the array:


			rxn of interest  ---->
      t      __________________________
    steps   |
	  |	    |   array of affinities/energys
	  |	    |   
	  V	    |   
		    |   
		    |   
		    |
	"""
	os.chdir(project_path)
	out_folder = project_path + '/' + name 							#	directory to house energy.csv used for vis
	mk_check_del_directory(out_folder) 
	os.chdir(out_folder) 											#	Move to output folder
	w = 0
	while w < len(index_1):
		out_name = index_1[w] + '_' + name + '.csv'
		df = pd.DataFrame(array[w], index=index_2, columns=index_3) # INDEX = row titles
		df.to_csv(out_name)
		w += 1




#################################################################
#####################    Helper Functions     ###################
def help():
	### Displayed to user of any of the arguments are enetered wrong
	print '5 arguements needed:' 
	print '		(1) Project folder name. Must be in local directory'
	print '		(2) dslop file name. Must be in ~1_SLOP/working_db'
	print '		(3) CON file name. Must be in project forlder identified in argument 1'
	print '		(4) SeaWater.3i file to be used for initial water rock reactions'
	print '		(5) Rock database to use for 6i generation'	

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

def find_highest_data0(file_list, db, l):
	###	called with the list of all data0's of a type (tde or tdf), database path (db), and the sublist to be searched (tde or tdf currently)
	temp = [] 												#	list of ints dervied from file suffixes
	for f in file_list:
		temp.append(f[(f.find('R')+1):])					#	Grab all file name characters after the 'R'
	temp = [int(item) for item in temp] 					#	convert to ints from the initial strings
	path = db + 'data0.' + l + '.R' + str(max(temp)) 		#	find the largest and rebuild path
	return path

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

def ck_for_empty_file(f):
	if 	os.stat(f).st_size == 0:
		print 'file: ' + str(f) + ' is empty.'
		sys.exit()

def runs_script_and_wait(args, wait_time, wait_file):	#	arguments sent to pOpen, file the appearence of whic is being waited for to keep going.
	### Puts the code to sleep while waiting for an output from a Popen call.
	### The wait file must be in the local directory set before this function was called.
	# with open(os.devnull, 'w') as fp:  							#	use of devnull will allow me to supress written output, and hopefull speed things up.
		# p = Popen(args, stdout=fp).wait()
	# call(args)
	with open(os.devnull, 'w') as fp: 								#	use of devnull will allow me to supress written output, and hopefull speed things up.
		Popen(args, stdout=fp).wait()								#	wait should give time for process to complete 

	sleep(wait_time)

	# for file in wait_file:
	# 	while not os.path.isfile(file):
			# sleep(0.1)

	return

def format_e(n, p):
	####     n = number, p = precisions
    return "%0.*E"%(p,n)

    # return b.split('E')[0].rstrip('0').rstrip('.') + 'E' + b.split('E')[1]



#################################################################
#####################     Call Main    ##########################


if __name__=='__main__':
    main(sys.argv[1:]) 						#	send through all arguemens except for the program name itself, sitting in position 0



#################################################################
##################   Troubleshooting Notes    ###################


