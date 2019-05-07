#!/usr/bin/env python
# -*- coding: utf-8 -*-
from rosetta import *
init()
import sys
import os
reload(sys)
from USR import *
import random
import math
from cleanfile import *
from pose_conformation import *
import time
def index_to_energy_level(index):
	"""this function used for 用于列表索引表示的能量级和列表内容能量级之间的转换

	#Author: QCQ
	
	"""
	return (index + 1) * 2

def energy_real_level_to_energy_level(real_level):
	"""这个函数用于把求的的蛋白质能量值进行归一化处理（离散化）
	"""
	return int(math.fabs(real_level) / 2) + 1

def energy_real_level_to_index(real_level):
	"""把求得的蛋白质能量转换成索引相对应的能量级
	"""
	return int(math.fabs(real_level) / 2)
	
def energy_level_to_index(level):
	"""这个函数和第一个函数功能相反
	"""
	return level / 2 -1
	
def usr_to_pose_key(usr):
	"""把USR转换成字符串，可以作为字典的键用来索引
	把小房子进行归一化处理，采用5个单位为一个小房子
	"""
	return str(int(usr[0]/5))+str(int(usr[1]/5))+str(int(usr[2]/5))
	
def usr_to_cell(usr):
	"""把USR转换成一个整数，利于判断是不是在同一个房间
	"""
	return int(str(int(usr[0]/5))+str(int(usr[1]/5))+str(int(usr[2]/5)))


	
	
	

