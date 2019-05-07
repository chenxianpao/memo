#!/usr/bin/env python
# -*- coding: utf-8 -*-
from rosetta import *
init()
import sys
import os
import shutil
reload(sys)
from USR import *
from perturbation import *
from pylab import *
import random
#from cleanfile import *
import time
from read_fragment import *

pdb_num = 200                          #初始种群个数
mutatorFactor = 0.5                       #变异系数  (useless right now)
crossFactor = 0.2							#交叉系数 (useless right now)
times = 1       #program run times
generations = 30000  #the number of the population's generations

#The parameters used in this program which can be changed depends\
#-on different protein.
protein_1ENH_clean='protein/1ENH/1ENH.clean1.pdb'
sequence_1ENH = 'RPRTAFSSEQLARLKREFNENRYLTERRRQQLSSELGLNEAQIKIWFQNKRAKI'
fragment_file_1ENH='protein/1ENH/aat000_03_05.200_v1_3'

protein_1GB1_clean='protein/1GB1/1GB1.clean1.pdb'
sequence_1GB1 = 'MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE'
fragment_file_1GB1='protein/1GB1/aat000_03_05.200_v1_3'

protein_1GYZ_clean='protein/1GYZ/1GYZ.clean1.pdb'
sequence_1GYZ = 'WIARINAAVRAYGLNYSTFINGLKKAGIELDRKILADMAVRDPQAFEQVVNKVKEALQVQ'
fragment_file_1GYZ = 'protein/1GYZ/aat000_03_05.200_v1_3'

protein_1I6C_clean='protein/1I6C/1I6C.clean1.pdb'
sequence_1I6C = 'KLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGNSSSG'
fragment_file_1I6C = 'protein/1I6C/aat000_03_05.200_v1_3'

protein_1VII_clean='protein/1VII/1VII.clean1.pdb'
sequence_1VII = 'MLSDEDFKAVFGMTRSAFANLPLWKQQNLKKEKGLF'
fragment_file_1VII='protein/1VII/aat000_03_05.200_v1_3'
#too short no fragment library
protein_1L2Y_output='protein/1L2Y/1L2Y.clean.pdb'  
protein_1L2Y_clean='1L2Y/1L2Y.clean1.pdb'
protien_1L2Y_input = '1L2Y/1L2Y.pdb'
sequence_1L2Y = 'NLYIQWLKDGGPSSGRPPPS'
fragment_file_1L2Y ='1L2Y/aat000_03_05.200_v1_3'

protein_4ICB_clean='protein/4ICB/4ICB.clean1.pdb'
sequence_4ICB = 'MKSPEELKGIFEKYAAKEGDPNQLSKEELKLLLQTEFPSLLKGPSTLDELFEELDKNGDGEVSFEEFQVLVKKISQ'
fragment_file_4ICB='protein/4ICB/aat000_03_05.200_v1_3'
#cxp's protein start
protein_1BBO_clean='protein/1BBO/1BBO.clean1.pdb'
sequence_1BBO = 'KYICEECGIRAKKPSMLKKHIRTHTDVRPYHCTYCNFSFKTKGNLTKHMKSKAHSKK'
fragment_file_1BBO='protein/1BBO/aat000_03_05.200_v1_3'

protein_1MN3_clean='protein/1MN3/1MN3.clean1.pdb'
sequence_1MN3 = 'SSLIKKIEENERKDTLNTLQNMFPDMDPSLIEDVCIAKKSRIEPCVDALLSLSE'
fragment_file_1MN3='protein/1MN3/aat000_03_05.200_v1_3'

protein_2IMU_clean='protein/2IMU/2IMU.clean1.pdb'
sequence_2IMU = 'FGFKDIIRAIRRIAVPVVSTLFPPAAPLAHAIGEGVDYLLGDEAQA'
fragment_file_2IMU='protein/2IMU/aat000_03_05.200_v1_3'

protein_2JUJ_clean='protein/2JUJ/2JUJ.clean1.pdb'
sequence_2JUJ = 'ATASPQLSSEIENLMSQGYSYQDIQKALVIAQNNIEMAKNILREFVSISSPAHVAT'
fragment_file_2JUJ='protein/2JUJ/aat000_03_05.200_v1_3'
#end cxp

protein_pdb_file_path = protein_1ENH_clean
protein_sequence = sequence_1ENH
protein_fragment_file_path = fragment_file_1ENH

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
current_usr_vector = []  #store usr vector
current_carmsd_value = [0 for i in range(0, pdb_num)]  #store CA_rmsd of pose object
current_energy_score = [0 for i in range(0, pdb_num)]  #store energy of pose object

#fullatom->backbone(centroid)
for i in range(0, pdb_num):
	to_centroid.apply(Parent[i])

#create new folder to hold protein's information
#-if exist then delete it
file_directory_information = protein_pdb_file_path.split(".")[0].split('/')[-1] + '/'
if os.path.exists(file_directory_information):
        shutil.rmtree(file_directory_information)
os.mkdir(file_directory_information)
#draw the figure
figure(num = 1)

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
	child_data = open(file_directory_information + "child_data.txt", "w")
	populations_data = open(file_directory_information + "populations_data.txt", "w")
	f_test = open(file_directory_information + "test.txt", "w")
	for times in range(0, 1):
                
		time1_start = time.clock()  #统计程序第一次运行200代所用的时间

		#initialization procedure
		for i in range(0, pdb_num):

                        read_fragment(Parent[i], protein_fragment_file_path)
                        current_carmsd_value[i] = CA_rmsd(com_pose, Parent[i])  #in case not to calculation CA_rmsd repeatedly
                        current_energy_score[i] = scorefxn(Parent[i])  #in case not to calculation energy repeatedly
                        #the initial populations' point drawing in green
                        scatter(current_carmsd_value[i],current_energy_score[i],1, color ='green')  
                        #usr_child.set_protein(Parent[i])
                        #current_usr_vector.append(usr_child.all_in_one())  #in case not to calculation usr vector repeatedly
                        f_test.write(str(current_carmsd_value[i]))
                        f_test.write('\t')
                        f_test.write(str(current_energy_score[i]))
                        f_test.write('\n')
                #use the below two value to draw figure        
                max_energy = int(max(current_energy_score)) + 1
                max_rmsd = int(max(current_carmsd_value)) + 1
                
                #program procedure
		for p in range(0, generations):
			print "%s %d times' %d-generation calculation's %d generation "%(protein_pdb_file_path.split(".")[0].split('/')[-1], \
                                                                                         times + 1,generations, p + 1),
			time1_in_start = time.clock()  #统计程序第一次运行1代所用的时间

			#DE procedure
			
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
                        #fragment_localMC_simple(Child, intermediate_pose, accepted_pose, scorefxn, reject_times, mover_3mer)

			#USR push
                        '''
                        usr_child.set_protein(Child)
                        usr_child_vector = usr_child.all_in_one()
                        for index_parent in range(0,pdb_num):
                                #usr_child.set_protein(Parent[index_parent])
                                current_usr_value[index_parent] = usr(usr_child_vector, current_usr_vector[index_parent])
			max_index = current_usr_value.index(max(current_usr_value))
			'''
                        max_index = current_energy_score.index(max(current_energy_score))
			energy_score_child = scorefxn(Child)
			CA_rmsd_score_child = CA_rmsd(com_pose, Child)
			#only the child replace father then change energy and carmsd
			if current_energy_score[max_index] > energy_score_child:  
                                Parent[max_index].assign(Child)
                                current_carmsd_value[max_index] = CA_rmsd_score_child
                                current_energy_score[max_index] = energy_score_child
                                #current_usr_vector[max_index] = usr_child_vector

                        #draw the child point in red                
                        scatter(CA_rmsd_score_child,energy_score_child,1, color ='red')
                        #output all Child data
                                
                        child_data.write(str(CA_rmsd_score_child))
                        child_data.write("\t")
                        child_data.write(str(energy_score_child))
                        child_data.write("\n")
                        if(CA_rmsd_score_child < 6):
                                child_output_name = file_directory_information + "CA_rmsd" + str(CA_rmsd_score_child) + \
                                                        "+" + str(energy_score_child) + "child_output.pdb"
                                Child.dump_pdb(child_output_name)
                                
                                        
                        time1_in_end = time.clock()
                        print	"Time used:%f"%(time1_in_end - time1_in_start)
			
                #output populations data
		#populations_data.write("The %d-generation"%(p + 1))
		#populations_data.write("\n")
                for index_output in range(0, pdb_num):
                        populations_data.write(str(current_carmsd_value[index_output]))              
                        populations_data.write("\t")
                        populations_data.write(str(current_energy_score[index_output]))
                        populations_data.write("\n")

                        #the last populations' point drawing in blue
                        scatter(current_carmsd_value[index_output],current_energy_score[index_output],1, color ='blue')
                        if(current_carmsd_value[index_output] < 6):
                                populations_output_name = file_directory_information + "CA_rmsd" + str(current_carmsd_value[index_output])\
                                                        + "+" + str(current_energy_score[index_output]) + "father.pdb"
                                Parent[index_output].dump_pdb(populations_output_name)
                                                                
		time1_end = time.clock()
		print "The %d times' %s all %d-generation used %f time to compleste"%(times, \
                                                                                       protein_pdb_file_path.split(".")[0].split('/')[-1], \
                                                                                       generations, time1_end - time1_start)

                #xticks([dexi*1 for dexi in range(0,25)])
                ax = gca()
                ax.spines['right'].set_color('none')
                ax.spines['top'].set_color('none')
                ax.xaxis.set_ticks_position('bottom')
                ax.spines['bottom'].set_position(('data',0))
                ax.yaxis.set_ticks_position('left')
                ax.spines['left'].set_position(('data',0))
                xticks([dexi*1 for dexi in range(0,max_rmsd)])
                yticks([dexi*1 for dexi in range(0,max_energy,50)])
                xlabel('CA_rmsd')
                #yab=str(energy_num)
                ylabel('energy')
                title('The comparison between %s ensemble and standard protein'%(protein_pdb_file_path.split(".")[0].split('/')[-1]),fontsize=15,color='red')
                grid(True)
		savefig(file_directory_information + 'hello1.png',dpi=90)
		
	time_all_end = time.clock()
	print "This program used total time is:%f"%(time_all_end - time_all_start)
	#figure(figsize = (8,6),dpi=80)
finally:
	child_data.close()
	populations_data.close()
	f_test.close()
        	
			
			
			
		
	
	



	
		
