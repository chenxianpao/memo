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
#import math
x_array=[]
y_array=[]
num=1000
M=100
generations=100
#why the range so bigger?
#usr_bound = [(0,50)]
usr_bound = [(0,4*np.pi)]
matrix_tree=[]
#The next list do not need initialization
#convertedusr_vector_value=[0 for i in range(0, num)] 
convertedusr_vector_value = []
#maybe sth wrong here, who told you the support vector is 4 dimensions!?
#SupportVector_value = [np.array([[0,0,0,0]]) for i in range(0, num)]
SupportVector_value = []
#unit_matrix = AcquireRootMatrix(2,usr_bound,M)
unit_energy1=math.sin(1 * (usr_bound[0][1] - usr_bound[0][0]))+M
unit_energy2=math.sin(0)+M
unit_matrix=np.array([[1.0 / unit_energy1, 0, unit_energy1],[0, 1.0 / unit_energy2, unit_energy2]])

print unit_matrix
child_data = open("child_rmsd_data.txt", "w")
#for i in range(0,num):
for i,x in enumerate(np.linspace(0, 13, num, endpoint = True)):
    #x=random.uniform(0,13)
    #x = random.random() * 13
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
#print x_array
#for i in range(0, num):
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
                   else:
                       ge+=1
                   totallen-=1               
              else:
                   print 'break',len(matrix_tree)
                   break 
#need to re-write the code
figure(num = 1)
x_sin = np.linspace(usr_bound[0][0], usr_bound[0][1], 200, endpoint = True)
y_sin = np.sin(x_sin)
new_x_array = []
new_y_array = []

for index,content in enumerate(np.linspace(usr_bound[0][0], usr_bound[0][1], 1000)):
    new_x_array.append(content)
    new_vector = ConvertUsr((content,), usr_bound)
    for all_tree in range(len(matrix_tree)):
        #print all_tree
        if IsRegion(matrix_tree[all_tree], new_vector):
            print 'in'

            Estimated_Value = EstimatedValue(matrix_tree[all_tree], new_vector, M)
            
            new_y_array.append(Estimated_Value)
            break
#添加哨兵
x_array.append((100000,))
new_x_array.append(100000)
print len(new_x_array), len(new_y_array)
plot(x_sin, y_sin, color = 'green')
temp_x = []
temp_y = []
index1 = 0
index2 = 0
for i in range(len(x_array) + len(new_x_array) - 2):
    if new_x_array[index1] <= x_array[index2][0]:
        temp_x.append(new_x_array[index1])
        temp_y.append(new_y_array[index1])
        index1+=1
    else:
        temp_x.append(x_array[index2][0])
        temp_y.append(y_array[index2])
        index2+=1

print len(temp_x), len(temp_y)            
            
plot(temp_x, temp_y, color = 'red')
show()
            
    

child_data.close()   
