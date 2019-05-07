
import numpy as np
'''
root_matrix = []
root_matrix.append(np.array([[1,1,1,1],[1,1,1,1],[1,1,1,1],[1,1,1,1]]))
root_matrix.append(np.array([[2,2,2,2],[2,2,2,2],[2,2,2,2],[2,2,2,2]]))
root_matrix.append(np.array([[3,3,3,3],[3,3,3,3],[3,3,3,3],[3,3,3,3]]))
print root_matrix
root_matrix.remove(np.array([[2,2,2,2],[2,2,2,2],[2,2,2,2],[2,2,2,2]]))
print root_matrix
root_matrix.append(np.array([[1,1,1,1],[1,1,1,1],[1,1,1,1],[1,1,1,1]]))
root_matrix.append(np.array([[1,1,1,1],[1,1,1,1],[1,1,1,1],[1,1,1,1]]))
'''
'''
usr_matrix=np.array([[1,2],[3,4]])
lk_vector=[]
for index in range(usr_matrix.shape[0]):
    lk_vector.append(usr_matrix[index,:])
    print lk_vector[index][1] 
    
    temp_x.append(lk_vector[index][0]*13)
                    print temp_x(xx)
'''

from ConditionsJudgev10 import *
y=[0,0]

matrix=np.array([[1,1],[2,2]])
usr=np.array([[1,2]])
M=100
#print EstimatedValue(y, matrix,usr,M)



'''
a = np.array([[9.83754679e-03,1.37325169e-04,1.00251913e+02],[0.00000000e+00,1.00000000e+00,1.00000000e+02]])
print a
print a[0,0:-1] 
a[1,1] = 10

print a
print ConditionOne(a)
'''
'''
matrix=np.array([[1,2,3],[4,5,6]])
temporary_matrix=matrix
diagonal_vector = matrix.diagonal()
temporary_matrix.sort(axis = 0)
    #get the second small items in the temporary matrix to form a vector
    #second_small_vector = temporary_matrix[1:2,]
   # third_small_vector = temporary_matrix[2:3,]
third_small_vector = temporary_matrix[0:1, 0:-1]
print diagonal_vector
print third_small_vector
compare_vector = (diagonal_vector >=2 third_small_vector)           
print compare_vector.all()
'''
'''
support = np.mat([0.2/ 201, 0.8 / 201, 0.0])
matrix1 = np.mat([[4.97512438e-04,4.47761194e-03,2.01000000e+02],[0.00000000e+00,4.97512438e-03,2.01000000e+02]])
matrix2 = np.mat([[4.97512438e-03,   0.00000000e+00 ,  2.01000000e+02],[  4.97512438e-04,   4.47761194e-03 ,  2.01000000e+02]])
print ConditionTwo(matrix1,support)
print ConditionTwo(matrix2,support)
'''
convertedusr_vector=np.mat([[1,2]])
#convertedusr_vector=[1,2]
diagonal_vector=[]
diagonal_vector.append(1)
diagonal_vector.append(2)
print convertedusr_vector
print diagonal_vector
print (convertedusr_vector[0, 0] / diagonal_vector[1])
'''
a=[]

a.append(1)
a.append(2)
a.pop(1)
a.insert(5,1)
print a
'''
a=np.mat([[0.2,0.8,1],[0.1,0.9,1]])
vecotor=np.array([[ 0.10808081 ,0.89191919]])
print IsRegion(a,vecotor)