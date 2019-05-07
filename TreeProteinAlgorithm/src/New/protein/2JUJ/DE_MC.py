#!/usr/bin/env python
# -*- coding: utf-8 -*-
from rosetta import *
init()
import sys
import os
reload(sys)
from USR import *
from pylab import *
import random
#from cleanfile import *
import time
pdb_num=150                          #初始种群个数
mutatorFactor=0.5                       #变异系数
crossFactor=0.2							#交叉系数
times=1000
#代数
temp_pose=Pose()
temp_score=0.0
reject_times=0
com_pose=pose_from_pdb("2JUJ.clean.pdb")         #比对对象
scorefxn=create_score_function("score3")      #score3设置
fragset=ConstantLengthFragSet(3)
fragset.read_fragment_file("aat000_03_05.200_v1_3")
usr1 = UltrafastShapeRecognition(com_pose)
usr_com = UltrafastShapeRecognition(com_pose)
usr_vector = usr_com.all_in_one()
movemap=MoveMap()
movemap.set_bb(True)
KT=0.0
scorefxn_hv=ScoreFunction()
scorefxn_hv.set_weight(hbond_lr_bb,1.0)
scorefxn_hv.set_weight(vdw,1.0)
ncycles=50
mover_3mer=ClassicFragmentMover(fragset,movemap)
Parent=[0 for I in range(0,pdb_num)]
Child=[0 for J in range(0,pdb_num)]
for i in range(0,pdb_num):
	Parent[i]=pose_from_sequence("ATASPQLSSEIENLMSQGYSYQDIQKALVIAQNNIEMAKNILREFVSISSPAHVAT","fa_standard")  
for i in range(0,pdb_num):
	Child[i]=pose_from_sequence("ATASPQLSSEIENLMSQGYSYQDIQKALVIAQNNIEMAKNILREFVSISSPAHVAT","fa_standard")
switch=SwitchResidueTypeSetMover("centroid")   #质心switch
recoverer = ReturnSidechainMover(Parent[0])     #从质心recover成原来的状态
time_all_start =time.clock()  #统计程序运行所有时间
current_score = [0 for i in range(0,pdb_num)]  #用于存储当前的种群里面的能量值，用于求解最大能量值
sequence_size=Parent[1].total_residue()
small_mover=SmallMover(movemap,KT,5)
small_mover.angle_max(5)
for i in range(0,pdb_num):
	switch.apply(Parent[i])
	switch.apply(Child[i])	

try:
	if(os.path.isfile("F:/360Downloads/python/1GYZ/DE_buildup200.txt")):
		os.remove("F:/360Downloads/python/1GYZ/DE_buildup200.txt")
	fsock = open("DE_USR.txt","a")
	for generation in range(0,1):
              
		time1_start = time.clock()  #统计程序第一次运行200代所用的时间	
		for i in range(0,pdb_num):
			for j in range(1,Parent[1].total_residue()+1):
				a=random.uniform(-180,180)
				b=random.uniform(-180,180)
				#c=random.uniform(-180,180)
				Parent[i].set_psi(j,a)
				Parent[i].set_phi(j,b)
				Parent[i].set_omega(j,180)#产生一个随机种群
		for p in range(0,times):
			print "The %d times' 300-generation calculation's %d generation "%(generation, p),
			time1_in_start = time.clock()  #统计程序第一次运行1代所用的时间
			for i in range(0,pdb_num):
				while True:
					a=random.randint(0,pdb_num-1)	
					b=random.randint(0,pdb_num-1)	
					c=random.randint(0,pdb_num-1)
					if(a!=b)and(a!=c)and(b!=c)and(a!=i)and(b!=i)and(c!=i):
						break
			        rand = random.randint(1,Parent[1].total_residue())
				for j in range(1,Parent[1].total_residue()+1):
                                        if random.random() < crossFactor or j == Parent[1].total_residue():
                                                
                                                angle_phi=Parent[a].phi(rand)+mutatorFactor*(Parent[b].phi(rand)-Parent[c].phi(rand))
                                                angle_psi=Parent[a].psi(rand)+mutatorFactor*(Parent[b].psi(rand)-Parent[c].psi(rand))
                                                #angle_omega=Parent[a].omega(j)+mutatorFactor*(Parent[b].omega(j)-Parent[c].omega(j))		
                                
                                                Child[i].set_phi(rand,angle_phi)
                                                Child[i].set_psi(rand,angle_psi)
                                                #Child[i].set_omega(j,angle_omega)            #变异
                                        else :
                                                Child[i].set_phi(rand,Parent[i].phi(rand))
                                                Child[i].set_psi(rand,Parent[i].psi(rand))
                                                
                                        rand = (rand + 1)%Parent[1].total_residue()
                                        if(rand == 0):
                                                rand = Parent[1].total_residue()
                               
			for o in range(0,pdb_num):
                                #mc=MonteCarlo(Child[o],scorefxn_hv,KT)
				for f_t in range(0,10):
                                        mover_3mer.apply(Child[o])    
					

					#ran_num=random.random()
					#if(ran_num>crossFactor):
					
					#mc.boltzmann(Child[o])
                                mc=MonteCarlo(Child[o],scorefxn,KT)
                                for q in range(0,ncycles):
                                        mover_3mer.apply(Child[o])    
                                        mc.boltzmann(Child[o])
					
				
					
			for l in range(0,pdb_num):
                                #usr2 = UltrafastShapeRecognition(Parent[l])
                                usr_child = UltrafastShapeRecognition(Child[l])
                                usr_child_vector = usr_child.all_in_one()
                                for b in range(0,pdb_num):
                                        usr_child.set_protein(Parent[b])
                                        current_score[b] = usr(usr_child_vector, usr_child.all_in_one())
				#switch.apply(Parent[l])
				#switch.apply(Child[l])
                                #current_score[l] = CA_rmsd(Parent[l],Child[i])
				#current_score[l]=usr_child
				max_index=current_score.index(max(current_score))
				#Parent[current_score.index(max(current_score))].assign(Child[i])
				if(scorefxn(Parent[max_index])>scorefxn(Child[l])):
                                        Parent[max_index].assign(Child[l])
				#if(((scorefxn(Parent[l]))>scorefxn(Child[l]))):
					#Parent[l].assign(Child[l])                  #将对象child赋值给parent
				#recoverer.apply(Parent[l])						#选择						
				#recoverer.apply(Child[l])
			time1_in_end = time.clock()
			print	"Time used:%f"%(time1_in_end - time1_in_start)
               
                        for i in range(0,pdb_num):
                                                                        
                                #switch.apply(Parent[i])
                                        saveout=sys.stdout
                                        sys.stdout=fsock
                                        
                                        CA_rmsd_score = CA_rmsd(com_pose,Child[i])
                                        sys.stdout.write(str(CA_rmsd_score))
                              
                                        sys.stdout.write("\t")
                                        energy_score = scorefxn(Child[i])
                                        if(CA_rmsd_score<6):
                                                output_name="CA_rmsd"+str(CA_rmsd_score)+"+"+str(energy_score)+".pdb"
                                                Parent[i].dump_pdb(output_name)
                                        	
                                        sys.stdout.write(str(scorefxn(Parent[i])))
                                        sys.stdout.write("\n")
                                        sys.stdout=saveout
                                #fsock.close()			
                                #recoverer.apply(Parent[i])             #输出最终种群
		time1_end = time.clock()
		print "The %d times'2JUJ all 300-generation used %f time to compleste"%(generation,time1_end - time1_start)
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
	fsock.close()
        	
			
			
			
		
	
	



	
		
