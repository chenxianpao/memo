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
pdb_num=150                            #初始种群个数
mutatorFactor=0.5                       #变异系数
crossFactor=0.2							#交叉系数
times=1000
#代数
#clean_file("1GYZ.pdb","1GYZ.clean.pdb")
com_pose=pose_from_pdb("1GYZ.clean.pdb")         #比对对象
scorefxn=create_score_function("score3")      #score3设置
fragset=ConstantLengthFragSet(3)
fragset.read_fragment_file("aat000_03_05.200_v1_3")
usr1 = UltrafastShapeRecognition(com_pose)
usr_com = UltrafastShapeRecognition(com_pose)
usr_vector = usr_com.all_in_one()
movemap=MoveMap()
movemap.set_bb(True)
KT=1.0
scorefxn_hv=ScoreFunction()
scorefxn_hv.set_weight(hbond_lr_bb,1.0)
scorefxn_hv.set_weight(vdw,1.0)
figure_num=1;
mover_3mer=ClassicFragmentMover(fragset,movemap)
Parent=[0 for I in range(0,pdb_num)]
Child=[0 for J in range(0,pdb_num)]
for i in range(0,pdb_num):
	Parent[i]=pose_from_sequence("WIARINAAVRAYGLNYSTFINGLKKAGIELDRKILADMAVRDPQAFEQVVNKVKEALQVQ","fa_standard")  
for i in range(0,pdb_num):
	Child[i]=pose_from_sequence("WIARINAAVRAYGLNYSTFINGLKKAGIELDRKILADMAVRDPQAFEQVVNKVKEALQVQ","fa_standard")
switch=SwitchResidueTypeSetMover("centroid")   #质心switch
recoverer = ReturnSidechainMover(Parent[0])     #从质心recover成原来的状态
time_all_start =time.clock()  #统计程序运行所有时间
current_score = [0 for i in range(0,pdb_num)]  #用于存储当前的种群里面的能量值，用于求解最大能量值
for i in range(0,pdb_num):
	switch.apply(Parent[i])
	switch.apply(Child[i])	

try:
	if(os.path.isfile("F:/360Downloads/python/1GYZ/DE_buildup200.txt")):
		os.remove("F:/360Downloads/python/1GYZ/DE_buildup200.txt")
	fsock = open("DE_USR.txt","a")
	for generation in range(0,1):
                Xaxes=[]
                Yaxes=[]
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
                                
                                                Child[i].set_psi(rand,angle_phi)
                                                Child[i].set_phi(rand,angle_psi)
                                                #Child[i].set_omega(j,angle_omega)            #变异
                                        else :
                                                Child[i].set_psi(rand,Parent[i].psi(rand))
                                                Child[i].set_psi(rand,Parent[i].psi(rand))
                                                
                                        rand = (rand + 1)%Parent[1].total_residue()
                                        if(rand == 0):
                                                rand = Parent[1].total_residue()
                               
			for o in range(0,pdb_num):
                                #mc=MonteCarlo(Child[o],scorefxn_hv,KT)
				for m in range(0,40):
					#ran_num=random.random()
					#if(ran_num>crossFactor):
					mover_3mer.apply(Child[o])
					#mc.boltzmann(Child[o])
					'''
					Child[o].set_psi(m,Parent[o].psi(m))         
					Child[o].set_phi(m,Parent[o].phi(m))
					Child[o].set_omega(m,Parent[o].omega(m))
					'''    #交叉
				
					
			for l in range(0,pdb_num):
                                usr2 = UltrafastShapeRecognition(Parent[l])
                                usr_parent = UltrafastShapeRecognition(Parent[l])
                                usr_parent_vector = usr_parent.all_in_one()
				#rmsd_parent=CA_rmsd(com_pose,Parent[l])
				#rmsd_child=CA_rmsd(com_pose,Child[l])
				#usr2.set_protein(Parent[l])
				#usr_parent = usr(usr_vector, usr1.all_in_one())
				usr2.set_protein(Child[l])
				usr_child = usr(usr_parent_vector, usr2.all_in_one())
				#switch.apply(Parent[l])
				#switch.apply(Child[l])
                                #current_score[l] = CA_rmsd(Parent[l],Child[i])
				current_score[l]=usr_child
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
                                usr1.set_protein(Parent[i])
                                usr1_vector = usr1.all_in_one()
                                sys.stdout=fsock
                                sys.stdout.write(str(usr(usr_vector, usr1_vector)))
                                Xaxes.append(usr(usr_vector, usr1_vector))
                                sys.stdout.write("\t")
                                CA_rmsd_score = CA_rmsd(com_pose,Parent[i])
                                sys.stdout.write(str(CA_rmsd_score))
                                Yaxes.append(CA_rmsd_score)
                                sys.stdout.write("\t")
                                energy_score = scorefxn(Parent[i])
                                if(CA_rmsd_score<4):
                                        output_name="CA_rmsd"+str(CA_rmsd_score)+"+"+str(energy_score)+".pdb"
                                        Parent[i].dump_pdb(output_name)
                                
                                usr_parent = usr(usr_vector, usr1_vector)	
                                if(usr_parent > 0.6):
                                        output_name="USR"+str(usr_parent)+"+"+str(energy_score)+".pdb"
                                        Parent[i].dump_pdb(output_name)
                                sys.stdout.write(str(scorefxn(Parent[i])))
                                sys.stdout.write("\n")
                                sys.stdout=saveout
                                #fsock.close()			
                                #recoverer.apply(Parent[i])             #输出最终种群
		time1_end = time.clock()
		print "The %d times' all 300-generation used %f time to compleste"%(generation,time1_end - time1_start)
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
	time_all_end = time.clock()
	print "This program used total time is:%f"%(time_all_end - time_all_start)
	#figure(figsize = (8,6),dpi=80)
finally:
	fsock.close()
        	
			
			
			
		
	
	



	
		
