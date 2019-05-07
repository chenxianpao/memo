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
from Draw import *
from pylab import *
from ConditionsJudge import *
import random
#from cleanfile import *
import time
from read_fragment import *

pdb_num = 10                          #��ʼ��Ⱥ����
mutatorFactor = 0.5                       #����ϵ��  (useless right now)
crossFactor = 0.2                            #����ϵ�� (useless right now)
times = 1       #program run times
generations = 10  #the number of the population's generations
usr_bound = [(0,100),(0,100),(0,100)] # usr_bound: this is a list, which items are tuple, in which contains every dimension
                    #'s lower and upper bound
#usr_bound = ([0,0,0,0],[100,100,100,100])
M=30000
ge=-1
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

protein_1FD4_clean='protein/1FD4/1FD4.clean1.pdb'
sequence_1FD4 = 'GIGDPVTCLKSGAICHPVFCPRRYKQIGTCGLPGTKCCKKP'
fragment_file_1FD4='protein/1FD4/aat000_03_05.200_v1_3'

protein_3GWL_clean='protein/3GWL/3GWL.clean1.pdb'
sequence_3GWL = 'GSHMLHWGPKYWRSLHLYAIFFSDAPSWKEKYEAIQWILNFIESLPCTRCQHHAFSYLTKNPLTLNNSEDFQYWTFAFHNNVNNRLNKKIISWSEYKNIYEQSILK'
fragment_file_3GWL='protein/3GWL/aat000_03_05.200_v1_3'
#end cxp

protein_pdb_file_path = protein_1ENH_clean
protein_sequence = sequence_1ENH
protein_fragment_file_path = fragment_file_1ENH

com_pose = pose_from_pdb(protein_pdb_file_path)         #standard pose object
sequence_length = com_pose.total_residue()       #get the protein's sequence length

#The parameters used in perturbation module
intermediate_pose = Pose()  #intermediary pose object in case change the original object
accepted_pose = Pose()      #accepted pose object in perturbation module
scorefxn=create_score_function("score3")      #score3����
reject_times = sequence_length - 2      #custom defined parameters

fragset = ConstantLengthFragSet(3)  #fragment's length
fragset.read_fragment_file(protein_fragment_file_path)  #fragments' file

usr_child = UltrafastShapeRecognition(com_pose)
usr_com = UltrafastShapeRecognition(com_pose)
usr_vector = usr_com.low_dimensional_geometric_projection()     #get standard pose objeect's 12-dimension-descriptor vector

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

time_all_start =time.clock()  #ͳ�Ƴ�����������ʱ��
current_usr_child_value = [0 for i in range(0, pdb_num)]
current_usr_value = [0 for i in range(0, pdb_num)]  #store usr value
current_usr_vector = []  #store usr vector
current_carmsd_value = [0 for i in range(0, pdb_num)]  #store CA_rmsd of pose object
current_energy_score = [0 for i in range(0, pdb_num)]  #store energy of pose object
convertedusr_vector_value = [np.array([[0,0,0,0]]) for i in range(0, pdb_num)]
SupportVector_value = [np.array([[0,0,0,0]]) for i in range(0, pdb_num)]
unit_matrix = AcquireRootMatrix(1000,usr_bound,M)
matrix_tree = []
leaf_num = 0 #�ϴ��������Ҷ�ӽڵ����
leaf_index=[]
temp_index=[]
dict_list=[]
divergency_flag=0

    
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
    child_rmsd_data = open(file_directory_information + "child_rmsd_data.txt", "w")
    populations_rmsd_data = open(file_directory_information + "populations_rmsd_data.txt", "w")
    child_usr_data = open(file_directory_information + "child_usr_data.txt", "w")
    populations_usr_data = open(file_directory_information + "populations_usr_data.txt", "w")
    f_test = open(file_directory_information + "test.txt", "w")
    f_test_usr = open(file_directory_information + "test_usr.txt", "w")
    for times in range(0, 1):
                
        time1_start = time.clock()  #ͳ�Ƴ����һ������200�����õ�ʱ��

        #initialization procedure
        for i in range(0, pdb_num):

                        read_fragment(Parent[i], protein_fragment_file_path)
                        current_carmsd_value[i] = CA_rmsd(com_pose, Parent[i])  #in case not to calculation CA_rmsd repeatedly
                        current_energy_score[i] = scorefxn(Parent[i])  #in case not to calculation energy repeatedly
                        #the initial populations' point drawing in green
                        #scatter(current_carmsd_value[i],current_energy_score[i],1, color ='green')  
                        usr_child.set_protein(Parent[i])
                        current_usr_vector.append(usr_child.low_dimensional_geometric_projection())  #in case not to calculation usr vector repeatedly
                        #in case not to calculation usr value repeatedly
                        #current_usr_value[i] = usr(usr_vector,current_usr_vector[i])
                        #output ca_rmsd and energy
                        f_test.write(str(current_carmsd_value[i]))
                        f_test.write('\t')
                        f_test.write(str(current_energy_score[i]))
                        f_test.write('\n')
                        #output usr and energy
                        f_test_usr.write(str(current_usr_value[i]))
                        f_test_usr.write('\t')
                        f_test_usr.write(str(current_energy_score[i]))
                        f_test_usr.write('\n')
                #use the below two value to draw figure        
               # max_energy = int(max(current_energy_score)) + 1
             #   max_rmsd = int(max(current_carmsd_value)) + 1
        for i in range(0, pdb_num):
            convertedusr_vector_value[i] = ConvertUsr(current_usr_vector[i], usr_bound)        
            SupportVector_value[i] = SupportVector(current_energy_score[i],convertedusr_vector_value[i],M)
        
        for i in range(0, pdb_num):
            print i
            if (i == 0):
                temp_matrix=MatrixReplace(unit_matrix, SupportVector_value[i])
                for j in range(0, 4):
                    flag_condition_one=ConditionOne(temp_matrix[j])
                    print temp_matrix[j]   
                    print flag_condition_one
                    if flag_condition_one == True:
                        matrix_tree.append(temp_matrix[j])
                        print temp_matrix[j]    
                        #leaf_index.append(matrix_tree.index(temp_matrix[j]))
                        
            else:
                  #for j in range(0,len(matrix_tree) ):
                  totallen=len(matrix_tree)
                  ge=0
                  while True:
                      if totallen>0:
                           #ge=ge+1
                           divergency_flag=0
                           flag_condition_two=ConditionTwo(matrix_tree[ge],SupportVector_value[i])
                           if flag_condition_two==False: 
                                  temp_matrix=MatrixReplace(matrix_tree[ge], SupportVector_value[i])
                                  #matrix_tree.remove(matrix_tree[ge])
                                  print "temp_matrix",temp_matrix
                                  for k in range(0,4):
                                      flag_condition_one=ConditionOne(temp_matrix[k])
                                      if flag_condition_one == True:
                                          #insert_index=leaf_index[j]+1
                                          divergency_flag+=1
                                          if(divergency_flag==1):
                                              matrix_tree.remove(matrix_tree[matrix_tree.index(matrix_tree[ge])])
                                            #  print "fist",matrix_tree[ge]
                                              #temptemp=matrix_tree[ge]
                                            #  matrix_tree.pop(ge)
                                              #print matrix_tree[ge]
                                          matrix_tree.insert(ge, temp_matrix[k])
                                          ge+=1
                           totallen-=1               
                      else:
                          break                    
                                          
                                     # temp_m={}
                                     # temp_m['index']=j+1
                                    #  temp_m['value']=temp_matrix[k]
                                    #  dict_list.append(temp_m)
                             
                                      #matrix_tree.insert(insert_index, temp_matrix[k])
                                      #temp_index.append(matrix_tree.index(temp_matrix[k]))
                              #del leaf_index[:]
                             # for k in range(0,len(temp_index)-1):
                          
                                 # leaf_index.append(temp_index[k])
                              #del temp_index[:]
                            #  break
finally:
    child_rmsd_data.close()
    populations_rmsd_data.close()
    child_usr_data.close()
    populations_usr_data.close()
    f_test.close()
    f_test_usr.close()                        
