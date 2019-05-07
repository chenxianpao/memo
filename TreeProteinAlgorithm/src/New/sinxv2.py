#! /usr/bin/env python
# -*- coding: utf-8 -*- 
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
import numpy as np
x_array=[]
y_array=[]
num=50
M=100
generations=100
#why the range so bigger?
#usr_bound = [(0,50)]
usr_bound = [(0,13)]
matrix_tree=[]
#The next list do not need initialization
#convertedusr_vector_value=[0 for i in range(0, num)] 
convertedusr_vector_value = []
#maybe sth wrong here, who told you the support vector is 4 dimensions!?
#SupportVector_value = [np.array([[0,0,0,0]]) for i in range(0, num)]
SupportVector_value = []
unit_matrix = AcquireRootMatrix(2,usr_bound,M)
child_data = open("child_rmsd_data.txt", "w")
for i in range(0,num):
    x=random.uniform(0,13)
    #the fromat here is not match the requirement of the function in ConditionsJudge.py
    #,in which need a tuple not a float value.
    #x_array.append(x)
    x_array.append((x,))
    y=math.sin(x)
    y_array.append(y)
#maybe this line no need!
#for i in range(0,num):
    convertedusr_vector_value.append(ConvertUsr(x_array[i], usr_bound))
    SupportVector_value.append(SupportVector(y_array[i],convertedusr_vector_value[i],M))
for i in range(0, num):
    print i
    if (i == 0):
        temp_matrix=MatrixReplace(unit_matrix, SupportVector_value[i])
        for j in range(0, 2):
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
                          for k in range(0,2):
                              flag_condition_one=ConditionOne(temp_matrix[k])
                              if flag_condition_one == True:
                                  #insert_index=leaf_index[j]+1
                                  divergency_flag+=1
                                  if(divergency_flag==1):
                                      
                                      matrix_tree.pop(ge)
                                     # print matrix_tree[ge].all()
                                  matrix_tree.insert(ge, temp_matrix[k])
                                  #print ge,totallen
                                  ge+=1
                   totallen-=1               
              else:
                   print 'break'
                   break 
#need to re-write the code
figure(num = 1)
x_sin = np.linspace(0, 4 * np.pi, 200)
y_sin = np.sin(x)
new_x_array = []
new_y_array = []
for index in np.linspace(0, 13, 100):
    new_x_array.append((index,))
    new_vector = ConvertUsr(new_x_array[index], usr_bound)
    for all_tree in range(len(matrix_tree)):
        if IsRegion(matrix[all_tree], new_vector):
            print 'in'
            Estimated_Value = EstimatedValue(matrix_tree[all_tree], new_vector, M)
            new_y_array.append(Estimated_Value)
plot(x_sin, y_sin, color = 'green')
plot(new_x_array, new_y_array, color = 'red')
show()
            
    

'''     
for p in range(0,generations):
    new_x=random.uniform(0,13)
    #the same problem
    #new_vector=ConvertUsr(new_x, usr_bound)
    new_vector=ConvertUsr((new_x,), usr_bound)
    #here is wrong!
    #y=math.sin(x)
    y = math.sin(new_x)
    for all_tree in range(0,len(matrix_tree)):
        IsRegion_flag=IsRegion(matrix_tree[all_tree], new_vector)
        print IsRegion_flag
        if(IsRegion_flag==True):
            print "in"
            Estimated_Value=EstimatedValue(matrix_tree[all_tree],new_vector,M)
            child_data.write(str(Estimated_Value))
            child_data.write("\t")
            child_data.write(str(new_x))
            child_data.write("\t")
            child_random=random.randint(0,13)
            if y_array[child_random]> Estimated_Value:
                  child_energy=math.sin(new_x)
                  if  y_array[child_random]>child_energy:
                      x_array[child_random]=new_x
                      #CA_rmsd_score_child = CA_rmsd(com_pose, Child)
                     #current_carmsd_value[Parent_index] = CA_rmsd_score_child
                      #current_energy_score[Parent_index] = child_energy
                      #if(CA_rmsd_score_child < 6):
                       #   child_output_name = file_directory_information + "CA_rmsd" + str(CA_rmsd_score_child) + \
                          #                  "+" + str(child_energy) + "child_output.pdb"
                       #   Child.dump_pdb(child_output_name)
                      #location_matrix=matrix_tree[alll]
                      #temp_matrix=temp_matrix=MatrixReplace(matrix_tree[alll], convert_child_usr)
                      totallen=len(matrix_tree)
                      ge=0
                      while True:
                          if totallen>0:
                               #ge=ge+1
                               divergency_flag=0
                               flag_condition_two=ConditionTwo(matrix_tree[ge],new_vector)
                               if flag_condition_two==False: 
                                      temp_matrix=MatrixReplace(matrix_tree[ge], new_vector)
                                      #matrix_tree.remove(matrix_tree[ge])
                                      for k in range(0,2):
                                          flag_condition_one=ConditionOne(temp_matrix[k])
                                          if flag_condition_one == True:
                                              #insert_index=leaf_index[j]+1
                                              divergency_flag+=1
                                              if(divergency_flag==1):
                                                  
                                                  matrix_tree.pop(ge)
                                                #  print matrix_tree[ge].all()
                                              matrix_tree.insert(ge, temp_matrix[k])
                                              ge+=1
                               totallen-=1               
                          else:
                              break 
            break
'''
child_data.close()   
