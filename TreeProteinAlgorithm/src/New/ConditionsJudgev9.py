#!/usr/bin/python  
# -*- coding: utf-8 -*- 
'''
Created on 2014-5-5

@author: qcq
'''
import rosetta
import numpy as np
import math
import random
from USR import UltrafastShapeRecognition

def AcquireRootMatrix(energy_user, usr_bound, M):
    '''
    This function used for computing the Root matrix of the whole tree
    
    parameters:
        energy_user: the energy given by user as big as possible
        usr_bound: this is a list, which items are tuple, in which contains every dimension
                    's lower and upper bound
        M: The setted M value, which setted by user 
        
    return:
        root_matrix: the whole tree Root Matrix
    '''
    #===========================================================================
    # The thinking of this function---each dimension of the usr has its low and
    # high limits, with these information can construct n(n == usr's dimension)
    # vector, in which only one element is the high limits, others' elements are
    # low limits. like (hign, low, low, ...),(low, high, low ...)...
    #===========================================================================
    #get the dimension of the usr
    dimension = len(usr_bound)
    #set usr_matrix root_matrix to empty matrix
    usr_matrix = np.mat(np.empty((dimension + 1, dimension)))
    root_matrix = np.mat(np.empty((dimension + 1, dimension + 1)))
    #index stand for the row of the matrix, j stand for the column of the matrix
    for index, bound in enumerate(usr_bound):
        for j in range(usr_matrix.shape[1]): 
            if index == j:
                usr_matrix[index, j] =  bound[1]
            else:
                usr_matrix[index, j] = bound[0]
                
    for index in range(usr_matrix.shape[1]):
        usr_matrix[usr_matrix.shape[0] - 1, index] = usr_bound[index][1]
    for index in range(usr_matrix.shape[0]):
        root_matrix[index,:] = SupportVector(energy_user, ConvertUsr(usr_matrix.tolist()[index], usr_bound), M)
        
    return root_matrix    

def ConditionOne(matrix):
    '''
    This function used for computing whether matrix comply to the Condition One
    
    parameters:
        matrix: the input matrix, which should be a squre matrix N * N
        
    return:
        boolean: if matrix comply to Condition One, then return True, otherwise False 
    '''
    #copy the matrix to temporary matrix in order to not change the original matrix
    temporary_matrix = matrix.copy()
    #sort the temporary matrix by column, 
    temporary_matrix.sort(axis = 0)
    #get the diagonal items to form a vector
    diagonal_vector = matrix.diagonal()
    #get the second small items in the temporary matrix to form a vector
    #second_small_vector = temporary_matrix[1:2,]
    third_small_vector = temporary_matrix[2:3,]
    #judge if the diagonal vector's items less than the second small vector
    #===========================================================================
    #The reason to deal vector like this is because, In every column of matrix 
    #may be have two or more items are smallest items. So, if compare the smallest items formed
    #vector with diagonal items formed vector can occur the problem in which the diagonal
    #items formed vector is not the only smallest items in every column
    #用上述的第二小元素组成的向量参与比较是为了防止在矩阵的每一列之中有两个相同的极小值，从而达不到对角线元素组成的向量之中的元素是
    #矩阵每一列之中唯一的最小值
    #===========================================================================
     
    compare_vector = (diagonal_vector >= third_small_vector) 
               
    return compare_vector.all()

def ConditionTwo(matrix, support_vector):
    '''
    This function used for computing whether the matrix comply to the Condition Two
    
    parameters:
        matrix: the input matrix, which should be a squre matrix N * N
        support_vector: the computed support vector 1 * N
        
    return:
        boolean: if matrix comply to Condition Two, then return True, otherwise False 
    '''
    #get the diagonal items to form a vector
    diagonal_vector = matrix.diagonal()
    compare_vector = (diagonal_vector >= support_vector)
    
    return not compare_vector.all()

def IsRegion(matrix, convertedusr_vector):
    '''
    This function used for judging whether the convertedusr_vector within region of the
    matrix
    
    parameters:
        energy: This is list contains value coressponding row of the column
        matrix: the input matrix, which should be a squre matrix N * N
        convertedusr_vector: The input usr-convert vector 1 * N
        
    return:
        boolean: if element contains in matrix, then return True, otherwise False 
    '''
    flag = True
    #i is the column of the matrix
    #j is the row of the matrix
    last_column = max(matrix.shape)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[0]):
            if i != j:
                #将支撑向量转换到转换后的解空间
                xjj = matrix[j, j] * matrix[j, last_column - 1]
                xi = convertedusr_vector[0, i]
                xj = convertedusr_vector[0, j]
                xji = matrix[j, i] * matrix[j, last_column - 1]
                #change!
                if xjj * xi < xji * xj:
                    flag = False
                    return flag
                    
    return flag

def EstimatedValue(matrix, convertedusr_vector, M):
    '''
    This function used for getting the estimated value of the protein
    
    parameters:
        energy: this is a list, which items are f(x) value
        matrix: the input matrix, which should be a squre matrix N * N
        convertedusr_vector: The input usr-convert vector 1 * N
        M: The setted M value, which setted by user
        
    return:
        float: return the estimated value
    '''
    #i is the row of the matrix
    #j is the column of the matrix
    #===========================================================================
    # The code changed
    #===========================================================================
    '''
    max_temp_vector = []
    for i in range(matrix.shape[0]):
        min_temp_vector = []
        for j in range(matrix.shape[1]):
            min_temp_vector.append(convertedusr_vector[0, i] / matrix[i, j])
        #get the miniature value of the list
        max_temp_vector.append(min(min_temp_vector))
    #get the maximum value of the list
    estimate_value = max(max_temp_vector)
    
    return estimate_value - M
    '''
    #===========================================================================
    # new code ,added after
    #===========================================================================
    #get the diagonal items to form a vector
    diagonal_vector = matrix.diagonal()
    value_list = []
    for item in range(max(diagonal_vector.shape)):
        value_list.append(convertedusr_vector[0, item] / diagonal_vector[0, item])
        
    return max(value_list) - M
    
    

def MatrixReplace(matrix, support_vector):
    '''
    This function used for replace one column of the matrix with the support_vector
    
    parameters:
        matrix: the input matrix, which should be a squre matrix N * N
        support_vector: the computed support vector 1 * N
        
    retrun:
        list: return one list contains matrices which each row replacement by support_vector
                once
    '''
    
    result_matrix = []
    for i in range(matrix.shape[0]):
        #used the copied one, in order to keep the original matrix
        temp_matrix = matrix.copy()
        #let the matrix's i row equals support_vector
        temp_matrix[i,:] = support_vector
        result_matrix.append(temp_matrix)
        
    return result_matrix

def ConvertUsr(usr, usr_bound):
    '''
    This function convert usr list to convertedusr_vector
    
    parameters:
        usr: the input usr, which is a list contains muti-dimensional data of one protein
        usr_bound: this is a list, which items are tuple, in which contains every dimension
                    's lower and upper bound
                    
    return:
        matrix: return a matrix(vector) which is converted usr vector
    '''
    bound_sum = 0.0
    #initial the convertedusr_vector to 1 * (len(usr) + 1)
    convertedusr_vector = np.zeros((1, len(usr) + 1))
    for i in usr_bound:
        bound_sum = bound_sum + (i[1] - i[0])    
    for i, bound in zip(range(len(usr)), usr_bound):
        convertedusr_vector[0, i] = (usr[i] - bound[0]) * 1.0 / bound_sum
    #the last element of convertedusr_vector is the sum of the elements before this element.
    convertedusr_vector[0, len(usr)] = 1 - convertedusr_vector[:,:-1].sum()
    
    return convertedusr_vector
 
def SupportVector(energy, convertedusr_vector, M):
    '''
    This function used for computing the support vector
    
    parameters:
        energy: the protein energy calculated by the score function
        convertedusr_vector: the usr converted vector 1 * N
        M: The setted M value, which setted by user 
        
    return:
        matrix: return support matrix(vector)
    '''
    #initial the support vector
    support_vector = np.zeros((1, max(convertedusr_vector.shape)))
    temp_energy = energy + M
    for i in range(max(convertedusr_vector.shape)):
        #The added code
        support_vector[0, i] = convertedusr_vector[0, i] * 1.0 / temp_energy 
        #=======================================================================
        # The below code is changed to new one
        #=======================================================================
        '''
        if convertedusr_vector[0, i] == 0:
            #if the divisor is zero, then set the corresponding value to M * 100
            support_vector[0, i] = M * 100
        else:
            support_vector[0, i] = temp_energy * 1.0 / convertedusr_vector[0, i]
        '''
    
    return support_vector
        
    
    