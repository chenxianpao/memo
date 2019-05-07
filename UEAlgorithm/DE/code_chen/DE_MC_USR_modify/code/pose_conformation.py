#!/usr/bin/env python
# -*- coding: utf-8 -*-
from rosetta import *
import sys
import os
reload(sys)
from USR import *
import random
from cleanfile import *


class PoseConformation:
	"""这个类用于保存在迭代的过程中产生的pose（蛋白质）对象和其他信息

        #Author: QCQ
        
	"""
	def __init__(self, score, protein, usr, CA_rmsd):
                """参数说明：

                score:该蛋白质的能量得分
                usr:该蛋白质的简化USR
                """
		self.score = score
		self.protein = Pose()
		self.protein.assign(protein)
		self.usr = usr
		self.CA_rmsd = CA_rmsd

	def set_all(self, score, protein ,usr):
                self.score = score
                self.protein.assign(protein)
                self.usr = usr
                self.CA_rmsd = CA_rmsd
	
	
	def set_protein(self,protein):
		self.protein.assign(protein)
	
	def set_energy(self,score):
		self.score = score

	def set_usr(self,usr):
                self.usr = usr

        def set_CA_rmsd(self, CA_rmsd):
                self.CA_rmsd = CA_rmsd

	def get_protein(self):
                return self.protein
		
	def get_energy(self):
                return self.score

        def get_usr(self):
                return self.usr

        def get_CA_rmsd(self):
                return self.CA_rmsd

        def  __str__(self):
                return "energy："+str(self.score)+" USR:("+str(self.usr[0])+","+str(self.usr[1])+","+str(self.usr[2])+")"
		
		
	
		
	


