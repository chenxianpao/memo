#!/usr/bin/python  
# -*- coding: utf-8 -*- 
'''
Created on 2014-4-1

@author: qcq
'''
from protein_information import *
#import anneal
import time
import anneal_modify as anneal
import os

iteration_times = 2500
population_size = 200
file_path_directory_window = ''
file_path_directory_ubuntu = '/home/qcq/space_anneal_usr/'

time_start_of_the_code = time.clock()
for pose_clean_file, pose_name, pose_sequence, pose_fragment, number in zip(protein_clean_file_name, protein_name, sequence_name\
                                                                    , fragment_name, range(1, len(want_to_run) + 1)):

    if pose_name in want_to_run1:
        current_path = file_path_directory_ubuntu +  pose_name + '/'
        if not os.path.exists(current_path):
            os.mkdir(current_path)
        
        anneal.ConformationalAnnealing(population_num = population_size, file_directory = current_path, figure_num = number, iteration = iteration_times, \
                                    protein_clean = pose_clean_file, protein_name = pose_name,\
                                    protein_sequence = pose_sequence, protein_fragment = pose_fragment,\
									FA_num = 1, MC_num = len(pose_sequence) - 2,
                                    KT = 4.35)
        time_end_of_the_code = time.clock()
        hour = (time_end_of_the_code -time_start_of_the_code) / 3600
        minute = (hour - int(hour)) * 60
        second = (minute - int(minute)) * 60
        #print "This program used %d hours %d minutes"%(int(hour), minute)
        #used for time count for each protein's run time
        with open(file_path_directory_ubuntu + "time.txt","a") as f:
            f.write("%s used %d hours %d minutes %d seconds."%(pose_name, int(hour), minute , second))
            f.write("finished at %s\n"%(time.strftime('%Y-%m-%d-%H:%M.%S',time.localtime(time.time()))))
