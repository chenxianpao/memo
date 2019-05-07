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
num=100
M=0
generations=100
#why the range so bigger?
#usr_bound = [(0,50)]
usr_bound = [(0,10)]
matrix_tree=[]
#The next list do not need initialization
#convertedusr_vector_value=[0 for i in range(0, num)] 
convertedusr_vector_value = []
#maybe sth wrong here, who told you the support vector is 4 dimensions!?
#SupportVector_value = [np.array([[0,0,0,0]]) for i in range(0, num)]
SupportVector_value = []
#unit_matrix = AcquireRootMatrix(2,usr_bound,M)
unit_energy1=1 + M
unit_energy2=1 + M
unit_matrix=np.array([[1.0/unit_energy1,0,unit_energy1],[0,1.0/unit_energy2,unit_energy2]])

print unit_matrix
child_data = open("child_rmsd_data.txt", "w")
#for i in range(0,100):
#for i,x in  enumerate(np.linspace(1, 9, 9,endpoint=True)):
   # x=random.uniform(0,10)
    #x=[1,2,3,4,5,6,7,8,9]
for x in range(1,10):
    print x
    #the fromat here is not match the requirement of the function in ConditionsJudge.py
    #,in which need a tuple not a float value.
    #x_array.append(x)
    x_array.append((x,))
    y=1
    y_array.append(y)
#maybe this line no need!
#for i in range(0,num):
    convertedusr_vector_value.append(ConvertUsr(x_array[x-1], usr_bound))
    SupportVector_value.append(SupportVector(y_array[x-1],convertedusr_vector_value[x-1],M))
for i in range(0, 9):
    print i
    if (i == 0):
        temp_matrix=MatrixReplace(unit_matrix, SupportVector_value[i])
        for j in range(0, 2):
            flag_condition_one=ConditionOne(temp_matrix[j])
            #print temp_matrix[j]   
            #print flag_condition_one
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
                   if i==1:
                       print flag_condition_two
                       print totallen
                       print ge
                       
                       print SupportVector_value[i]
                       print matrix_tree[0]
                       print matrix_tree[1]
                      # print matrix_tree[2]
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
x_sin = np.linspace(0, 4 * np.pi, 200)
y_sin = [1 for i in range(0,200)]
new_x_array = []
new_y_array = []

for index,content in enumerate(np.linspace(1, 9, 100)):
    new_x_array.append(content)
    new_vector = ConvertUsr((content,), usr_bound)
    for all_tree in range(len(matrix_tree)):
        #print all_tree
        if IsRegion(matrix_tree[all_tree], new_vector):
            print 'in'
       
            Estimated_Value = EstimatedValue(matrix_tree[all_tree], new_vector, M)
            
            new_y_array.append(Estimated_Value)
            break
        
            #print 'notin'
print len(new_x_array), len(new_y_array)
#print len(x_sin),len(y_sin)
plot(x_sin, y_sin, color = 'green')
plot(new_x_array, new_y_array[0:len(new_x_array)], color = 'red')
show()
            
  
child_data.close()   
