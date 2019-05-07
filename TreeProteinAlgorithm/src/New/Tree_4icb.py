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

pdb_num = 150                          #��ʼ��Ⱥ����
mutatorFactor = 0.5                       #����ϵ��  (useless right now)
crossFactor = 0.2                            #����ϵ�� (useless right now)
times = 1       #program run times
generations = 10000  #the number of the population's generations
usr_bound = [(0,100),(0,100),(0,100)] # usr_bound: this is a list, which items are tuple, in which contains every dimension
                    #'s lower and upper bound
#usr_bound = ([0,0,0,0],[100,100,100,100])
M=30000

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

protein_pdb_file_path = protein_4ICB_clean
protein_sequence = sequence_4ICB
protein_fragment_file_path = fragment_file_4ICB

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
#unit_matrix = AcquireRootMatrix(1000,usr_bound,M)
#unit_energy1=1000+M
#unit_energy2=1000+M
#unit_matrix = np.array([[1.0 / unit_energy1, 0, unit_energy1],[0, 1.0 / unit_energy2, unit_energy2]])
matrix_tree = []
divergency_flag=0

    
#fullatom->backbone(centroid)
for i in range(0, pdb_num):
    to_centroid.apply(Parent[i])

#create new folder to hold protein's information
#-if exist then delete it
file_directory_information = protein_pdb_file_path.split(".")[0].split('/')[-1] + '/'
#if os.path.exists(file_directory_information):
   # shutil.rmtree(file_directory_information)
#os.mkdir(file_directory_information)
#draw the figure
figure(num = 1)

#below code which commented useless temporarily

try:
        
        child_rmsd_data = open(file_directory_information + "child_rmsd_data.txt", "w")
        populations_rmsd_data = open(file_directory_information + "populations_rmsd_data.txt", "w")
        #child_usr_data = open(file_directory_information + "child_usr_data.txt", "w")
        #populations_usr_data = open(file_directory_information + "populations_usr_data.txt", "w")
        f_test = open(file_directory_information + "test.txt", "w")
        #f_test_usr = open(file_directory_information + "test_usr.txt", "w")
        
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
                           # print usr_child.low_dimensional_geometric_projection()
                            current_usr_vector.append(usr_child.low_dimensional_geometric_projection())  #in case not to calculation usr vector repeatedly
                            #in case not to calculation usr value repeatedly
                            #current_usr_value[i] = usr(usr_vector,current_usr_vector[i])
                            #output ca_rmsd and energy
                            f_test.write(str(current_carmsd_value[i]))
                            f_test.write('\t')
                            f_test.write(str(current_energy_score[i]))
                            f_test.write('\n')
                            #output usr and energy
                            #f_test_usr.write(str(current_usr_value[i]))
                            #f_test_usr.write('\t')
                           # f_test_usr.write(str(current_energy_score[i]))
                           # f_test_usr.write('\n')
                    #use the below two value to draw figure        
            #max_energy = int(max(current_energy_score)) + 1
                 #   max_rmsd = int(max(current_carmsd_value)) + 1
            print max(current_energy_score)
           # raw_input()
            unit_energy1=max(current_energy_score)+M+500
            unit_energy2=max(current_energy_score)+M+500
            unit_energy3=max(current_energy_score)+M+500
            unit_energy4=max(current_energy_score)+M+500
            unit_matrix = np.array([[1.0 / unit_energy1, 0, 0, 0, unit_energy1],[0, 1.0 / unit_energy2, 0, 0, unit_energy2],[0,0,1.0/unit_energy3,0,unit_energy3],[0,0,0,1.0/unit_energy4,unit_energy4]])
            for i in range(0, pdb_num):
                convertedusr_vector_value[i] = ConvertUsr(current_usr_vector[i], usr_bound)        
                SupportVector_value[i] = SupportVector(current_energy_score[i],convertedusr_vector_value[i],M)
            #raw_input()
            for i in range(0, pdb_num):
                print "pdbnum",i
                if (i == 0):
                    temp_matrix=MatrixReplace(unit_matrix, SupportVector_value[i])
                    for j in range(0, 4):
                        flag_condition_one=ConditionOne(temp_matrix[j])
                        #print temp_matrix[j]   
                        print flag_condition_one
                        if flag_condition_one == True:
                            matrix_tree.append(temp_matrix[j])
                            #print temp_matrix[j]    
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
                                      for k in range(0,4):
                                          flag_condition_one=ConditionOne(temp_matrix[k])
                                          if flag_condition_one == True:
                                              #insert_index=leaf_index[j]+1
                                              divergency_flag+=1
                                              if(divergency_flag==1):
                                                  
                                                  matrix_tree.pop(ge)
                                                 # print matrix_tree[ge].all()
                                              matrix_tree.insert(ge, temp_matrix[k])
                                              ge+=1
                               else:
                                   ge+=1
                               totallen-=1
                               
                                                  
                          else:
                              print 'break',len(matrix_tree)
                              break         
            #raw_input()           
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
                Parent_index=a
                for range_index in range(range_min, range_max + 1):
                                    Child.set_phi(range_index, Parent[b].phi(range_index))
                                    Child.set_psi(range_index, Parent[b].psi(range_index))
                                    Child.set_omega(range_index, Parent[b].omega(range_index))
                            
                            #mutation procedure
                mover_3mer.apply(Child)
                usr_child.set_protein(Child)
                usr_child_vector = usr_child.low_dimensional_geometric_projection()
                convert_child_usr=ConvertUsr(usr_child_vector, usr_bound)
                
                for all_tree in range(0,len(matrix_tree)):
                    IsRegion_flag=IsRegion(matrix_tree[all_tree], convert_child_usr)
                    #print IsRegion_flag
                    #print matrix_tree[all_tree]
                    if(IsRegion_flag==True):
                        print "in"
                        Estimated_Value=EstimatedValue(matrix_tree[all_tree],convert_child_usr,M)
                        if current_energy_score[Parent_index]> Estimated_Value:
                              child_energy=scorefxn(Child)
                              support_vector_child=SupportVector(child_energy,convert_child_usr,M)
                              if current_energy_score[Parent_index]>child_energy:
                                  Parent[Parent_index].assign(Child)
                                  CA_rmsd_score_child = CA_rmsd(com_pose, Child)
                                  current_carmsd_value[Parent_index] = CA_rmsd_score_child
                                  current_energy_score[Parent_index] = child_energy
                                  
                                  if(CA_rmsd_score_child < 6):
                                      child_output_name = file_directory_information + "CA_rmsd" + str(CA_rmsd_score_child) + \
                                                        "+" + str(child_energy) + "child_output.pdb"
                                      Child.dump_pdb(child_output_name)
                                  #location_matrix=matrix_tree[alll]
                                  #temp_matrix=temp_matrix=MatrixReplace(matrix_tree[alll], convert_child_usr)
                                  totallen=len(matrix_tree)
                                  ge=0
                                  while True:
                                      if totallen>0:
                                           #ge=ge+1
                                           divergency_flag=0
                                           flag_condition_two=ConditionTwo(matrix_tree[ge],support_vector_child)
                                           if flag_condition_two==False: 
                                                  temp_matrix=MatrixReplace(matrix_tree[ge], support_vector_child)
                                                  #matrix_tree.remove(matrix_tree[ge])
                                                  for k in range(0,4):
                                                      flag_condition_one=ConditionOne(temp_matrix[k])
                                                      if flag_condition_one == True:
                                                          #insert_index=leaf_index[j]+1
                                                          divergency_flag+=1
                                                          if(divergency_flag==1):
                                                              
                                                              matrix_tree.pop(ge)
                                                            #  print matrix_tree[ge].all()
                                                          matrix_tree.insert(ge, temp_matrix[k])
                                                          ge+=1

                                           else:

                                               ge+=1
                                               
                                           totallen-=1               
                                      else:
                                          break 
                        break
                
                    #else:
                       # break                
                time1_in_end = time.clock()
                print    "Time used:%f"%(time1_in_end - time1_in_start)
                CA_rmsd_score_child = CA_rmsd(com_pose, Child)
                child_rmsd_data.write(str(CA_rmsd_score_child))
                child_rmsd_data.write("\t")
                child_energy=scorefxn(Child)
                child_rmsd_data.write(str(child_energy))
                child_rmsd_data.write("\n")
                
                        #if(current_carmsd_value[index_output] < 6):
                         #       populations_output_name = file_directory_information + "CA_rmsd" + str(current_carmsd_value[index_output])\
                          #                              + "+" + str(current_energy_score[index_output]) + "father.pdb"
                           #     Parent[index_output].dump_pdb(populations_output_name)
                                                                
        for index_output in range(0, pdb_num):
                        populations_rmsd_data.write(str(current_carmsd_value[index_output]))              
                        populations_rmsd_data.write("\t")
                        populations_rmsd_data.write(str(current_energy_score[index_output]))
                        populations_rmsd_data.write("\n")   
        time_all_end = time.clock()
        print "This program used total time is:%f"%(time_all_end - time_all_start)  
        f_test.write("This program used total time is:%f"%(time_all_end - time_all_start))
        f_test.write('\n')                                 
finally:
    

        child_rmsd_data.close()
        populations_rmsd_data.close()
     #   child_usr_data.close()
       # populations_usr_data.close()
        f_test.close()
         
file_in=[]
file_in.append(file_directory_information + "test.txt")
file_in.append(file_directory_information + "child_rmsd_data.txt")
file_in.append(file_directory_information + "populations_rmsd_data.txt")
drawing_picture(file_in,file_directory_information + "rmsd", 1, 12, 500, protein_pdb_file_path)

