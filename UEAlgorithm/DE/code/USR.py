#!/usr/bin/env python
# -*- coding: utf-8 -*-
import math
from rosetta import *
def usr(protein1,protein2):
	"""This funtion is used for USR(Ultrafast shape recognition) score between two given proteins
   
   求解两个蛋白质之间的相识度得分
   parameter:
   protein1:the usr encoded vector
   protein2: the USR encoded vector

   #Author: QCQ
        """
	temp = [0 for i in range(0,12)]
	sum_temp = 0.0
	for j in range(0,12):
		temp[j] = protein1[j] - protein2[j]
		sum_temp += temp[j]*temp[j]
	score = 1.0/(1 + 1.0/12 *(math.sqrt(sum_temp)))
	return score
	
class UltrafastShapeRecognition:
   """This class used for USR implement
   synx:
   protein:the given protein
   """
   def __init__(self,protein):
           self.protein=Pose()
   	   self.protein.assign(protein)
   
   def set_protein(self,protein):
           #self.protein = Pose()
   	   self.protein.assign(protein)
   	   
   def centroid(self):
	"""this function is used for calculating the centroid of the given self.protein
	
	求解蛋白质分之的质心
	parameter:
	self.protein: the given self.protein which used to calcuate the cenroid
	"""
	centroid = [0.0,0.0,0.0]
	for i in range(1,self.protein.total_residue() + 1):
		N_xyz = self.protein.residue(i).xyz("N")
		#print N_xyz
		CA_xyz = self.protein.residue(i).xyz("CA")
		C_xyz = self.protein.residue(i).xyz("C")
		O_xyz = self.protein.residue(i).xyz("O")
		centroid[0] += N_xyz.x + CA_xyz.x + C_xyz.x + O_xyz.x
		centroid[1] += N_xyz.y + CA_xyz.y + C_xyz.x + O_xyz.y
		centroid[2] += N_xyz.z + CA_xyz.z + C_xyz.z + O_xyz.z
	centroid[0] = centroid[0] / (self.protein.total_residue() * 4)
        centroid[1] = centroid[1] / (self.protein.total_residue() * 4)
        centroid[2] = centroid[2] / (self.protein.total_residue() * 4)
        centroid_vector = numeric.xyzVector_double(centroid[0], centroid[1], centroid[2])#将列表转换成rosetta里面的Vector类型
        return centroid_vector
        
   def nearest_farest_atom_from_centroid(self,centroid):
	"""this function is used for calculating the atom's cartesian coordinate(笛卡尔坐标) which nearest the centroid 
	
	求解离质心最近和最远原子的坐标,返回有两个元素的链表，第一个是最近，第二个是最远
	parameter:
	pose: the given portein
	centroid: the given self.protein's centroid
	"""
	nearest = self.protein.residue(1).xyz("N")
	farest = self.protein.residue(1).xyz("N")
	
	try:
	    #f = open("1.txt","w")  #这个文件的内容容易调试时查看得到的是否是最远和最近距离
	    for i in range(1, self.protein.total_residue()+1):
		N_xyz = self.protein.residue(i).xyz("N")
		CA_xyz = self.protein.residue(i).xyz("CA")
		C_xyz = self.protein.residue(i).xyz("C")
		O_xyz = self.protein.residue(i).xyz("O")
		N_distance = (N_xyz - centroid).norm
		CA_distance = (CA_xyz - centroid).norm
		C_distance = (C_xyz - centroid).norm
		O_distance = (O_xyz - centroid).norm
		current_near_distance = (nearest - centroid).norm  #得到当前的最近距离
	        current_far_distance = (farest - centroid).norm #得到当前的最远距离
		
		if(N_distance == min([N_distance, CA_distance, C_distance, O_distance])):
			if(N_distance < current_near_distance):
				 nearest = N_xyz   #用于记录当前离质心最近原子的坐标
				 flag_near = i      #调试使用，用于记录当前离质心最近原子所在的残基号码，下同
		elif(CA_distance == min([N_distance, CA_distance, C_distance, O_distance])):
			if(CA_distance < current_near_distance):
				nearest = CA_xyz
				flag_near = i
		elif(C_distance == min([N_distance, CA_distance, C_distance, O_distance])):
			if(C_distance < current_near_distance):
				nearest = C_xyz
				flag_near = i
		elif(O_distance == min([N_distance, CA_distance, C_distance, O_distance])):
			if(O_distance < current_near_distance):
				nearest = O_xyz    # calculate the nearest atom from centroid
				flag_near = i
		else:
			continue
		
		if(N_distance == max([N_distance, CA_distance, C_distance, O_distance])):
			if(N_distance > current_far_distance):
			        farest = N_xyz   #用于记录当前离质心最远原子的坐标
			        flag_far = i   #调试使用，用于记录当前离质心最远原子所在的残基号码，下同
		elif(CA_distance == max([N_distance, CA_distance, C_distance, O_distance])):
			if(CA_distance > current_far_distance):
			        farest = CA_xyz
			        flag_far = i
		elif(C_distance == max([N_distance, CA_distance, C_distance, O_distance])):
			if(C_distance > current_far_distance):
			        farest = C_xyz
			        flag_far = i
		elif(O_distance == max([N_distance, CA_distance, C_distance, O_distance])):
			if(O_distance > current_far_distance):
			        farest = O_xyz      # calculate the farest atom from centroid
			        flag_far = i
		else:
			continue
		#f.write("residue %d N(%f,%f,%f)CA(%f,%f,%f)C(%f,%f,%f)O(%f,%f,%f) N_distance: %f CA_distance: %f C_distance: %f O_distance: %f"%(i,N_xyz.x, N_xyz.y,N_xyz.z,CA_xyz.x, CA_xyz.y,CA_xyz.z,C_xyz.x,C_xyz.y,C_xyz.z,O_xyz.x, O_xyz.y, O_xyz.z,N_distance,CA_distance,C_distance,O_distance))
		#f.write("\n")
	finally:
	     pass
	    #print flag_near,(nearest - centroid).norm,flag_far,(farest - centroid).norm 	
            #f.close()	
	return [nearest,farest]

   def farest_from_farest_atom_from_centroid(self,farest_atom_from_centroid):
	"""this function is used for calculating the the atom's cartesian coordinate(笛卡尔坐标) that farest from the atom which farest from the centroid
	求距离质心最远原子最远原子的笛卡尔坐标
	
	parameter:
	self.protein: the given self.protein
	farest_atom_from_centroid:距离质心最远的原子的坐标
	"""
	farest = self.protein.residue(1).xyz("N")
	try:
	        #f = open("2.txt","w")
	        for i in range(1, self.protein.total_residue()+1):
	        	N_xyz = self.protein.residue(i).xyz("N")
		        CA_xyz = self.protein.residue(i).xyz("CA")
		        C_xyz = self.protein.residue(i).xyz("C")
		        O_xyz = self.protein.residue(i).xyz("O")
		        N_distance = (N_xyz - farest_atom_from_centroid).norm
		        CA_distance = (CA_xyz - farest_atom_from_centroid).norm
		        C_distance = (C_xyz - farest_atom_from_centroid).norm
		        O_distance = (O_xyz - farest_atom_from_centroid).norm
	                current_far_distance = (farest - farest_atom_from_centroid).norm #得到当前的最远距离
		
		        if(N_distance == max([N_distance, CA_distance, C_distance, O_distance])):
		        	if(N_distance > current_far_distance):
		        		farest = N_xyz   #用于记录当前离质心最远原子最远原子的坐标
				        flag_far = i      #调试使用，用于记录当前离质心最远原子最远原子所在的残基号码，下同
		        elif(CA_distance == max([N_distance, CA_distance, C_distance, O_distance])):
			        if(CA_distance > current_far_distance):
				        farest = CA_xyz
				        flag_far = i
		        elif(C_distance == max([N_distance, CA_distance, C_distance, O_distance])):
			        if(C_distance > current_far_distance):
				        farest = C_xyz
				        flag_far = i
		        elif(O_distance == max([N_distance, CA_distance, C_distance, O_distance])):
			        if(O_distance > current_far_distance):
				        farest = O_xyz    # calculate the nearest atom from centroid
				        flag_far = i
		        else:
			        continue
	                #f.write("residue %d N(%f,%f,%f)CA(%f,%f,%f)C(%f,%f,%f)O(%f,%f,%f) N_distance: %f CA_distance: %f C_distance: %f O_distance: %f"%(i,N_xyz.x, N_xyz.y,N_xyz.z,CA_xyz.x, CA_xyz.y,CA_xyz.z,C_xyz.x,C_xyz.y,C_xyz.z,O_xyz.x, O_xyz.y, O_xyz.z,N_distance,CA_distance,C_distance,O_distance))		
		        #f.write("\n")
            
	finally:
		pass
		#print flag_far,current_far_distance
		#f.close()
	return farest

   def all_in_one_atom(self):
	"""this function combine the above three methods centroid,nearest__farest_atom_from_centroid,farest_from_farest_atom_from_centroid
	
	这个函数结合了以上三个函数用于一步求解出质心，离质心最近原子，离质心最远原子，离质心最远原子最远原子的坐标返回一个链表
	第一个元素是质心，第二个是离质心最近原子，第三个离质心最远，第四个离质心最远原子最远原子
	parameter:
	self.protein: The given self.protein
	"""
	centroid_atom = self.centroid()
	temp_atom = self.nearest_farest_atom_from_centroid(centroid_atom)
	nearest_atom = temp_atom[0]
	farest_atom = temp_atom[1]
	farest_farest_atom = self.farest_from_farest_atom_from_centroid(farest_atom)
	return [centroid_atom, nearest_atom, farest_atom, farest_farest_atom]

   def twelve_moment(self,four_reference_atom_coordinate):
	"""this function is uesd to calculate the 12 moment vector
	
	这个函数求解USR所需要的12维的向量
	parameter:
	self.protein:the given self.protein
	four_reference_atom_coordinate:四个要求的原子坐标（列表），
	第一个元素是质心，第二个是离质心最近原子，第三个离质心最远，第四个离质心最远原子最远原子
	"""
       	try:
		#f = open("3.txt","w")
		atom_xyz_list = []
		distance_centroid_list = []
		distance_nearest_list = []
		distance_farest_list = []
		distance_farest_farest_list = []
		
	
                average_distance_centroid = 0.0 #求各个原子与质心的距离之和的平均值
		average_distance_nearest = 0.0
		average_distance_farest = 0.0
		average_distance_farest_farest = 0.0
			
		variance_cenroid = 0.0 #求解各原子相较于质心的方差
		variance_nearest = 0.0
		variance_farest = 0.0
		variance_farest_farest = 0.0
		
		skewness_cenroid = 0.0 #求解各原子相较于质心的偏度
		skewness_nearest = 0.0
		skewness_farest = 0.0
		skewness_farest_farest = 0.0
		
	        for i in range(1, self.protein.total_residue()+1):
	        	N_xyz = self.protein.residue(i).xyz("N")
		        CA_xyz = self.protein.residue(i).xyz("CA")
		        C_xyz = self.protein.residue(i).xyz("C")
		        O_xyz = self.protein.residue(i).xyz("O")
		        atom_xyz_list.append(N_xyz)    #将值进行储存，以利于下边求取方法不用再次做相同的运算
		        atom_xyz_list.append(CA_xyz)
		        atom_xyz_list.append(C_xyz)
		        atom_xyz_list.append(O_xyz)
		        
		        N_distance_centroid = (N_xyz - four_reference_atom_coordinate[0]).norm  #求解各原子到质心的距离
		        CA_distance_centroid = (CA_xyz - four_reference_atom_coordinate[0]).norm
		        C_distance_centroid = (C_xyz - four_reference_atom_coordinate[0]).norm
		        O_distance_centroid = (O_xyz - four_reference_atom_coordinate[0]).norm
		        distance_centroid_list.append(N_distance_centroid)
		        distance_centroid_list.append(CA_distance_centroid)
		        distance_centroid_list.append(C_distance_centroid)
		        distance_centroid_list.append(O_distance_centroid)
		        
		        N_distance_nearest = (N_xyz - four_reference_atom_coordinate[1]).norm  #求解各原子到最近原子的距离
		        CA_distance_nearest = (CA_xyz - four_reference_atom_coordinate[1]).norm
		        C_distance_nearest = (C_xyz - four_reference_atom_coordinate[1]).norm
		        O_distance_nearest = (O_xyz - four_reference_atom_coordinate[1]).norm
		        distance_nearest_list.append(N_distance_nearest)
		        distance_nearest_list.append(CA_distance_nearest)
		        distance_nearest_list.append(C_distance_nearest)
		        distance_nearest_list.append(O_distance_nearest)
		        
		        N_distance_farest = (N_xyz - four_reference_atom_coordinate[2]).norm  #求解各原子到最远原子的距离
		        CA_distance_farest = (CA_xyz - four_reference_atom_coordinate[2]).norm
		        C_distance_farest = (C_xyz - four_reference_atom_coordinate[2]).norm
		        O_distance_farest = (O_xyz - four_reference_atom_coordinate[2]).norm
		        distance_farest_list.append(N_distance_farest)
		        distance_farest_list.append(CA_distance_farest)
		        distance_farest_list.append(C_distance_farest)
		        distance_farest_list.append(O_distance_farest)
		        
		        N_distance_farest_farest = (N_xyz - four_reference_atom_coordinate[3]).norm  #求解各原子到最近原子的距离
		        CA_distance_farest_farest = (CA_xyz - four_reference_atom_coordinate[3]).norm
		        C_distance_farest_farest = (C_xyz - four_reference_atom_coordinate[3]).norm
		        O_distance_farest_farest = (O_xyz - four_reference_atom_coordinate[3]).norm
		        distance_farest_farest_list.append(N_distance_farest_farest)
		        distance_farest_farest_list.append(CA_distance_farest_farest)
		        distance_farest_farest_list.append(C_distance_farest_farest)
		        distance_farest_farest_list.append(O_distance_farest_farest)
		        

		        
		average_distance_centroid = sum(distance_centroid_list) / len(distance_centroid_list) #求各个原子与质心的距离之和的平均值
		average_distance_nearest = sum(distance_nearest_list) / len(distance_nearest_list)
		average_distance_farest = sum(distance_farest_list) / len(distance_farest_list)
		average_distance_farest_farest = sum(distance_farest_farest_list) / len(distance_farest_farest_list)
		
		sum_variance_centroid = 0.0   #获取求解方差过程中的和值
		sum_variance_nearest = 0.0
		sum_variance_farest =0.0
		sum_variance_farest_farest = 0.0
		
		sum_skewness_centroid = 0.0   #获取求解偏度过程中的和值
		sum_skewness_nearest = 0.0
		sum_skewness_farest =0.0
		sum_skewness_farest_farest = 0.0
	        for i,j,x,y in zip(distance_centroid_list,distance_nearest_list,distance_farest_list,distance_farest_farest_list):
	        	sum_variance_centroid += math.pow(i - average_distance_centroid,2)
	        	sum_variance_nearest += math.pow(j - average_distance_nearest,2)
	        	sum_variance_farest += math.pow(x - average_distance_farest,2)
	        	sum_variance_farest_farest += math.pow(y - average_distance_farest_farest,2)
	        	
	        	sum_skewness_centroid += math.pow(i - average_distance_centroid,3)
	        	sum_skewness_nearest += math.pow(j - average_distance_nearest,3)
	        	sum_skewness_farest += math.pow(x - average_distance_farest,3)
	        	sum_skewness_farest_farest += math.pow(y - average_distance_farest_farest,3)
	        	
	        	
	        variance_cenroid = sum_variance_centroid / len (distance_centroid_list)
	        variance_nearest = sum_variance_nearest / len (distance_nearest_list)
	        variance_farest = sum_variance_farest / len (distance_farest_list)
	        variance_farest_farest = sum_variance_farest_farest / len (distance_farest_farest_list)
	        
	        skewness_cenroid = sum_skewness_centroid / (len (distance_centroid_list) * math.pow(variance_cenroid,1.5))
	        skewness_nearest = sum_skewness_nearest / (len (distance_nearest_list) * math.pow(variance_nearest,1.5))
	        skewness_farest = sum_skewness_farest / (len (distance_farest_list) * math.pow(variance_farest,1.5))
	        skewness_farest_farest = sum_skewness_farest_farest / (len (distance_farest_farest_list) * math.pow(variance_farest_farest,1.5))
	        
	        result = []   #将计算得到的向量十二个维度写入一个列表
	        result.append(average_distance_centroid)
	        result.append(variance_cenroid)
	        result.append(skewness_cenroid)
	        result.append(average_distance_nearest)
	        result.append(variance_nearest)
	        result.append(skewness_nearest)
	        result.append(average_distance_farest)
	        result.append(variance_farest)
	        result.append(skewness_farest)
	        result.append(average_distance_farest_farest)
	        result.append(variance_farest_farest)
	        result.append(skewness_farest_farest)
	        return result
	        
	        
	        	
	        	
	        	
	        	
        finally:
        	pass
        	#f.close()
	        #print average_distance_centroid,average_distance_nearest,average_distance_farest,average_distance_farest_farest
	        #print len(atom_xyz_list),len(distance_centroid_list),len(distance_nearest_list),len(distance_nearest_list),len(distance_farest_farest_list)
	        #print self.protein.total_residue()
	        
   def all_in_one(self):
	"""直接由给定的蛋白质给出其12维的向量
	
	parameter:
	self.protein:the goven self.protein
	"""
	
	return self.twelve_moment(self.all_in_one_atom())

   def low_dimensional_geometric_projection(self):
   	"""this function used for calculation the simple vector comes from USR
   	
   	这个函数用于求解简化版的USR，只是取USR的四个坐标分量前三个坐标的第一个坐标分量
   	three projection coordinates are defined
        on each computed conformation: the mean atomic distance
        from the centroid (ctd), the mean atomic distance
        from the atom farthest from the centroid (fct), and the mean
        atomic distance from the atom farthest from fct (ftf)
        """
        vector_12 = self.twelve_moment(self.all_in_one_atom())
        vector_3 =[vector_12[0], vector_12[6], vector_12[9]]
        return vector_3
   
   def __str__(self):
        """return返回一个简化版USR的三维数字转变成字符串
        """

        None
           
       
        	
	
	
	        

		
