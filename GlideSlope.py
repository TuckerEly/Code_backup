########################################################################################################################################
######    GLideSlope.py:  To dynamically call EQ3/6 for largescale projects.  (here, for DCO/C-DEBI gale database work)    #############
#########################################      Tucker Ely, 10 Jan 2017   ###############################################################
########################################################################################################################################
###  	
### 	v1.   Just deal with calling and renaming outputs. Future version will be able to call other .py
###			files to deal with file building.
###
###
###				(0) Source correct environemnt 							(no dealt with yet, this must stll be done manually)
### 			(1) find all 3i/6i files of interest  					(done)
###				(2) Run them through EQ3/6 one at a time 				(done)
###				(3) Check for errors returned from csh 					# 	for the time being error checking here is ignored
###				(4) check for errors returned from EQ3/6 				# 	for the time being error checking is not handeled, other than to print.
###				(5) Rename outputs, continue. 							(done)
###
###


import os, re, sys, shutil, fileinput, getopt, argparse
import numpy as np
import pandas as pd
from itertools import *
from collections import *
from subprocess import * 
from time import gmtime, strftime


"""	EDIT BELOW	"""

EQ36_parent = "/Users/tuckerely/EQ3_6v8.0a/"

"""	EDIT ABOVE	"""




#################################################################
########################    Main    #############################
def main(argv):
	pwd = os.getcwd()													#	what is the present working directory											
	environemnt = EQ36_parent + "scripts/eq36cfg" 						#	environemental file to be sourced
	exe_directory = EQ36_parent + "bin/" 								#	location of executable scripts

	data0, file_type, input_folder, out = process_arguments(argv, pwd) 	#	pricess user arguments

	# my_env = os.environ.copy()[environemnt]

	file_name, file_list = read_inputs(file_type, pwd, input_folder) 	#	grab all 3i/6i files of interest
	fail_list = [] 														#	list of failed input files

	### Set correct executable
	if file_type == '.6i':
		exe = 'runeq6'
	if file_type == '.3i':
		exe = 'runeq3'


	"""Run all 3i/6i files"""
	w = 0
	while w < len(file_name): 											#	cycle through input files
		# print file_name[w]
		p_out = file_type[:-1] + 'p' 									#	pickup suffix
		o_out = file_type[:-1] + 'o' 									#	output suffix
		args = ['/bin/csh', exe_directory + exe, data0, input_folder + '/' + file_name[w]] 	#	Arguements to be run inside csh.
		with open(os.devnull, 'w') as fp: 								#	use of devnull will allow me to supress written output, and hopefull speed things up.
			Popen(args, stdout=fp).wait()								
		try:																				#	try to move file (if output is present)
			shutil.move('output', out + '/' + file_name[w][:-3] + o_out)		#	move and rename output for given 3i/6i run, to the output folder out				#	rename generaic 'output' to input.3/6o
			shutil.move('pickup', out + '/' + file_name[w][:-3] + p_out)		#	move and rename output for given 3i/6i run, to the output folder out				#	rename generaic 'output' to input.3/6o
		except:																				# 	if output is not present, add its parent name to the failed file list
			fail_list.append(file_name[w])
			print file_name[w], 'failed'
			w += 1 
			continue 
		w += 1 																				#	go to next file



	""" Deal with failed 3i / 6i files"""
	if len(fail_list) > 0:
		print 'The following inputs failed  ', fail_list
		fail_dir = out + '/fail' 					#	fail directory name within out folder
		os.makedirs(fail_dir)						#	make failed directory
		for item in fail_list: 						#	proform file move on all failed files, if they exist, hence the pass options
			try: 									#	try to move output. this cannot be assumed, but tried, becuase not all EQ36 failures produce and o-file.
				shutil.move(out + '/' + item[:-3] + o_out, fail_dir + '/' + item[:-3] + o_out)
			except:
				pass
			try: 									#	try to move output. this cannot be assumed, but tried, becuase not all failures produce and o-file.
				shutil.move(out + '/' + item[:-3] + p_out, fail_dir + '/' + item[:-3] + p_out)
			except:
				pass
	else:
		print 'All inputs seccessful'




#################################################################
#####################    Functions     ##########################


def read_inputs(file_type, pwd, input_folder):
	### this function can find all .6o files in all downstream folders, hence the additio of file_list
	file_name = [] 						#	file names
	file_list = []						# 	file names with paths
	for root, dirs, files in os.walk(pwd + '/' + input_folder):# this works. 
		for file in files:
			if file.endswith(file_type):
				file_name.append(file)
				file_list.append(os.path.join(root, file)) 					
	
	return file_name, file_list

def process_arguments(argv, pwd):
	### Grab Arguments
	try: 														#	test for the correct arguemnts in command line
		opts, args = getopt.getopt(argv,"h",["help"]) 								
	except getopt.GetoptError: 									#	if an unrecognized flag is passed, an exception is raised. right now only -h is recognized
		print 'Unrecognized argement.'
		sys.exit(2) 											#	2 here is the exit code and means (command line syntax errors)


	### Check for exceptions and errors in arguments.
	for opt, arg in opts: 										#	If user called help, then give it to them
		if opt in ("-h", "--help"):
			help()
			sys.exit()
	if len(args)<3:												#	if fewer than 3 arguments are give, FIX THAT SHIT !
		help() 													#	FIX THAT SHIT ! !  Notify user of correct arguments
	if not os.path.isfile(EQ36_parent + 'db/data1.' + args[0]):	#	check for data1 file assiciated with 3 letter suffix called
		print 'Associated data1 file not present in:'
		print '		' + EQ36_parent + 'db/'
		sys.exit() 												
	if args[1] != '3i' and args[1] != '6i': 					#	Check for 3i or 6i in input file type argument
		print 'Incorrect input file suffix:'					#	nope !
		print '		Must be 3i or 6i.'
		sys.exit()
	if not os.path.isdir(pwd + '/' + args[2]): 					#	Check for existance of folder containing the inputs
		print 'Input file directry not present in current folder.'
		sys.exit()


	### Set Arguements
	data0 = args[0] 											#	data0 suffix: data1 must exist in /db
	file_type = '.' + args[1]									#	3i or 6i file?
	input_folder = args[2]	 									#	folder within local: where are all of the input files stroed
	run_time = strftime("%Y-%m-%d_%H-%M-%S", gmtime()) 			#	runtime, for unique folder nameing
	out = input_folder + '_out_' + run_time 					#	new output folder name
	if not os.path.exists(out):									#	check if desired output folder exists
	    os.makedirs(out)										#	build 3/6op output folder 

	return data0, file_type, input_folder, out

def help():
	print '3 arguements needed:' 
	print '		(1) Data0 3 letter suffix such as ymp,'
	print '		(2) Two letter file type (3i or 6i)'
	print '		(3) Local folder name containing the input files,'


#################################################################
#####################     Call Main    ##########################


if __name__=='__main__':
    main(sys.argv[1:]) 						#	send through all arguemens except for the program name itself, sitting in position 0



#################################################################
##################   Troubleshooting Notes    ###################
###
###
###
### 		EQ36 trys to look in folders that dont exist for error messages:
###				Note: The following floating-point exceptions are signalling: IEEE_UNDERFLOW_FLAG
###				/usr/tmp/errchk27195: No such file or directory.
###				/usr/tmp/errchk27195: No such file or directory.
###				/usr/tmp/errchk27195: No such file or directory.
###				/usr/tmp/ifrun27195: No such file or directory.
###
###				These destinations only show up in the runeq3-6 scripts. Not in any of the .f or .h files
### 				So, they are either assumed by the scripts to be set up by the build process, or they
### 				are natural errorchecking components of csh, the shell t be used. Cant find anything
###					online however. Online documentation notes that /usr should not be written to, and that
###					usr/tmp is depreciated and no longer used.
###
###				Checked first use of EQ3 3o output with those provided by wolery, and they were identicle, except
###				for the 0.000's int he Temperature coefficients, which showed up as -0.000 in mine, and without the -
###				in his. Hmm. I should see if this effects the 6i files with T ranges, as it had no effect on the 3o 
###				file chemistries
###
###				Checked first use of EQ6 with this code to do all 6i test files provided. Oddly, the outputs are identicle
###				(i looked at rwitr.3o) in all aspects except the 'charge descrepancy' and the 'Relative charge descrepancy'
###				strangly though, no other charge issues are present, its as if the calculation of these two factors, from an 
###				identicle set of charges, is wrong. I should look at the fortran file that is responsible for this.
###
### 		Must find out relationship between parent and chid process, insofaras python interacts with csh:
### 		Can view the interaction has having happend in a single terminal session?
### 		Or does parent daughter only amply to single commands, (kind of like piping) unless i am also misunderstanding that
###



#################################################################
#########################   Sources    ##########################
### Dynamic handeling of csh calling:   
###			https://docs.python.org/2.6/library/subprocess.html
###			https://pymotw.com/2/subprocess/
###
### Error handeling:  
###			http://stackoverflow.com/questions/12335356/how-to-handle-error-exception-in-shell-script
###			Popen(['/bin/csh','-c','source', environemnt])




