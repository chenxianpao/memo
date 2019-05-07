#!/usr/bin/env python
# -*- coding: utf-8 -*-
from rosetta import *
init()
import sys
import os
reload(sys)
from USR import *
from perturbation import *
from pylab import *
import random
#from cleanfile import *
import time

pdb_num = 200                          #初始种群个数
mutatorFactor = 0.5                       #变异系数  (useless right now)
crossFactor = 0.2							#交叉系数 (useless right now)
times = 1       #program run times
generations = 1000  #the number of the population's generations

#The parameters used in this program which can be changed depends\
#-on different protein.
protein_pdb_file_path = "1GYZ.clean.pdb"
protein_sequence = "WIARINAAVRAYGLNYSTFINGLKKAGIELDRKILADMAVRDPQAFEQVVNKVKEALQVQ"
protein_fragment_file_path = "aat000_03_05.200_v1_3"

com_pose = pose_from_pdb(protein_pdb_file_path)         #standard pose object
sequence_length = com_pose.total_residue()       #get the protein's sequence length

#The parameters used in perturbation module
intermediate_pose = Pose()  #intermediary pose object in case change the original object
accepted_pose = Pose()      #accepted pose object in perturbation module
scorefxn=create_score_function("score3")      #score3设置
reject_times = sequence_length - 2      #custom defined parameters

fragset = ConstantLengthFragSet(3)  #fragment's length
fragset.read_fragment_file(protein_fragment_file_path)  #fragments' file

usr_child = UltrafastShapeRecognition(com_pose)
usr_com = UltrafastShapeRecognition(com_pose)
usr_vector = usr_com.all_in_one()     #get standard pose objeect's 12-dimension-descriptor vector

movemap_fragment = MoveMap()  #Used for fragment
movemap_fragment.set_bb(True)

movemap_initial = MoveMap()  #Used for initial pose object
movemap_initial.set_bb(False)

#KT = 0.0
short_inserts = 1
#ncycles = 200
mover_3mer = ClassicFragmentMover(fragset, movemap_fragment)
mover_3mer_initial = ClassicFragmentMover(fragset, movemap_initial)

Parent = [0 for I in range(0, pdb_num)]
Child = Pose()

for i in range(0, pdb_num):
	Parent[i] = pose_from_sequence(protein_sequence, "fa_standard")

#centroid<->fullatom
to_centroid = SwitchResidueTypeSetMover("centroid")   #centroid->fullatom convert
to_fullatom = SwitchResidueTypeSetMover("fa_standard")   #fullatom->centroid

time_all_start =time.clock()  #统计程序运行所有时间
current_usr_value = [0 for i in range(0, pdb_num)]  #store usr value
current_carmsd_value = [0 for i in range(0, pdb_num)]  #store CA_rmsd of pose object
current_energy_score = [0 for i in range(0, pdb_num)]  #store energy of pose object

#fullatom->backbone(centroid)
for i in range(0, pdb_num):
	to_centroid.apply(Parent[i])

#below code which commented useless temporarily
"""	
insert_short_frag = RepeatMover(mover_3mer,short_inserts)   #apply short_inserts times 3-mer fragment assmebly
folding_mover = SequenceMover()  #apply a sequence perturbation to protin conformation
folding_mover.add_mover(insert_short_frag)
mc=MonteCarlo(Child[1],scorefxn,KT)
trial = TrialMover(folding_mover,mc)
folding = RepeatMover(trial,ncycles)
"""

try:
	child_data = open("child_data.txt", "w")
	populations_data = open("populations_data.txt", "w")
	for times in range(0, 1):
                
		time1_start = time.clock()  #统计程序第一次运行200代所用的时间

		#initialization procedure
		for i in range(0, pdb_num):
			for j in range(1, sequence_length - 1):
                                movemap_initial.set_bb_true_range(j, j + 2)
                                mover_3mer_initial.apply(Parent[i])

                        current_carmsd_value[i] = CA_rmsd(com_pose, Parent[i])  #in case not to calculation CA_rmsd repeatedly
                        current_energy_score[i] = scorefxn(Parent[i])  #in case not to calculation energy repeatedly

                #program procedure
		for p in range(0, generations):
			print "%s %d times' %d-generation calculation's %d generation "%(protein_pdb_file_path.split(".")[0], \
                                                                                         times + 1,generations, p + 1),
			time1_in_start = time.clock()  #统计程序第一次运行1代所用的时间

			#DE procedure
			for i in range(0, pdb_num):
                                #cross procedure
				while True:
					a = random.randint(0, pdb_num-1)	
					b = random.randint(0, pdb_num-1)	
					#c=random.randint(0,pdb_num-1)
					if(a != b):  #and(a!=c)and(b!=c)and(a!=i)and(b!=i)and(c!=i):
						break

				#choose the range which will be crossed
                                while True:
					rand1 = random.randint(1, sequence_length)	
					rand2 = random.randint(1, sequence_length)	
					#c=random.randint(0,pdb_num-1)
					if(rand1 != rand2):  #and(a!=c)and(b!=c)and(a!=i)and(b!=i)and(c!=i):
						break
				range_min = min(rand1, rand2)
				range_max = max(rand1, rand2)
				Child.assign(Parent[a])
				for range_index in range(range_min, range_max + 1):
                                        Child.set_phi(range_index, Parent[b].phi(range_index))
                                        Child.set_psi(range_index, Parent[b].psi(range_index))
                                        Child.set_omega(range_index, Parent[b].omega(range_index))
                        
                                #mutation procedure
                                mover_3mer.apply(Child)

                                #FA+MC 
                                fragment_localMC_simple(Child, intermediate_pose, accepted_pose, scorefxn, reject_times, mover_3mer)

				#USR push
                                usr_child.set_protein(Child)
                                usr_child_vector = usr_child.all_in_one()
                                for index_parent in range(0,pdb_num):
                                        usr_child.set_protein(Parent[index_parent])
                                        current_usr_value[index_parent] = usr(usr_child_vector, usr_child.all_in_one())
				max_index = current_usr_value.index(max(current_usr_value))
				energy_score_child = scorefxn(Child)
				CA_rmsd_score_child = CA_rmsd(com_pose, Child)
				#only the child replace father then change energy and carmsd
				if current_energy_score[max_index] > energy_score_child:  
                                        Parent[max_index].assign(Child)
                                        current_carmsd_value[max_index] = CA_rmsd_score_child
                                        current_energy_score[max_index] = energy_score_child
                                        

                                #output all Child data
                                
                                child_data.write(str(CA_rmsd_score_child))
                                child_data.write("\t")
                                child_data.write(str(energy_score_child))
                                child_data.write("\n")
                                if(CA_rmsd_score_child < 6):
                                        child_output_name="CA_rmsd" + str(CA_rmsd_score_child) + \
                                                           "+" + str(energy_score_child) + "child_output.pdb"
                                        Child.dump_pdb(child_output_name)
                                
                                        
		time1_in_end = time.clock()
		print	"Time used:%f"%(time1_in_end - time1_in_start)
			
                #output populations data
		populations_data.write("The %d-generation"%(p + 1))
		populations_data.write("\n")
                for index_output in range(0,pdb_num):
                        populations_data.write(str(current_carmsd_value[index_output]))              
                        populations_data.write("\t")
                        populations_data.write(str(current_energy_score[index_output]))
                        populations_data.write("\n")
                        if(current_carmsd_value[index_output] < 6):
                                populations_output_name="CA_rmsd" + str(current_carmsd_value[index_output])\
                                                        + "+" + str(current_energy_score[index_output]) + "father.pdb"
                                Parent[i].dump_pdb(populations_output_name)
                                                                
		time1_end = time.clock()
		print "The %d times' %s all %d-generation used %f time to compleste"%(times, \
                                                                                       protein_pdb_file_path.split(".")[0], \
                                                                                       generations, time1_end - time1_start)
		'''
                figure(num=figure_num)
		scatter(Xaxes,Yaxes, 50, color ='green')
		for index_x,index_y in zip(Xaxes,Yaxes):
                        if index_x > 0.5 and index_y < 6:
                                scatter(index_x, index_y, 50, color = 'red')
		for index_x,index_y,index in zip(Xaxes,Yaxes,range(0,len(Xaxes))):
                        
                        annotate(str(index+1),xy=(index_x, index_y), xycoords='data')
		xlim(min(Xaxes)*0.9, max(Xaxes)*1.1)
		ylim(min(Yaxes)*0.9, max(Yaxes)*1.1)
		
		figname = "figure"+str(figure_num)
		figure_num = figure_num + 1
		savefig(figname)
		'''
	time_all_end = time.clock()
	print "This program used total time is:%f"%(time_all_end - time_all_start)
	#figure(figsize = (8,6),dpi=80)
finally:
	child_data.close()
	populations_data.close()
        	
			
			
			
		
	
	



	
		
