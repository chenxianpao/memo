#!/usr/bin/env python
# encoding: utf-8
'''
This program used for test the unix2dos function

Author:QCQ
'''

import sys
import os
from protein_information import *
from unix2dos import *

file_directory = '/home/qcq/python/'
if os.path.exists(file_directory):
	for index_T in range(1 , len(KT_list_input) + 1):
	    file_directory_T = file_directory + str(index_T) + '/'
	    if os.path.exists(file_directory_T):
		for index_FA in FA_list:
		    file_directory_FA = file_directory_T + str(index_FA) + '/'
		    if os.path.exists(file_directory_FA):
			for index_MC in MC_list:
			    file_directory_MC = file_directory_FA + str(index_MC) + '/'
			    if os.path.exists(file_directory_MC):
			        for protein in want_to_run:
			            file_directory_protein = file_directory_MC + protein + '/'
			            if os.path.exists(file_directory_protein):
			                file_deal_with = []
			                file_directory_data = file_directory_protein + 'data.txt'
			                file_directory_information = file_directory_protein + 'information.txt'
			                file_directory_test = file_directory_protein + 'test.txt'
			                file_directory_temperature = file_directory_protein + 'temperature.txt'
			                file_deal_with.append(file_directory_data)
			                file_deal_with.append(file_directory_information)
			                file_deal_with.append(file_directory_test)
			                file_deal_with.append(file_directory_temperature)
			                for index_file in file_deal_with:
					    if os.path.exists(index_file):
			                        unix2dos(index_file)  #unix to dos
						#dos2unix(index_file)
