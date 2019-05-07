#!/usr/bin/env python
# -*- coding: utf-8 -*-
from rosetta import *
init()
import sys
import os
import shutil
reload(sys)
from USR import *
import random
import math
from pylab import *  #画图所需
from cleanfile import *
from pose_conformation import *
from AllFunction import *
import time
import pdb
import optparse  #for sorting options
from read_fragment import *


###The data defination region
"""
computer_memory = 4  # This variable maybe used in following program
short_inserts = 1 #define times applied 3-mer fragment assembly
cycles = 10  ##定义进行多少次MC
iteration = 1000 #define the num of the conformations
file_directory_linux='/home/qcq/tree_based_memory/'#linux绝对目录
file_directory_windows = ''#windows相对目录（当前目录）
file_directory = file_directory_windows  #use windows
protein_directory='protein/1GYZ/1GYZ.clean1.pdb'
protein_name='1GYZ'
protein_sequence='WIARINAAVRAYGLNYSTFINGLKKAGIELDRKILADMAVRDPQAFEQVVNKVKEALQVQ'
fragment_file='protein/1GYZ/aat000_03_05.200_v1_3'

#Author: QCQ
modified in 2014/3/12   add temperature layer 
modified in 2014/3/13   add the function can run this program many times once
modified in 2014/3/14   fix a bug which can not deal with temperature_list
                        add change size to FA and MC
"""
def fold_protein(run_time = 1,KT_list = [1.0],file_directory = '',KT=1.0,protein_name='1GYZ',protein_sequence='WIARINAAVRAYGLNYSTFINGLKKAGIELDRKILADMAVRDPQAFEQVVNKVKEALQVQ',\
                 fragment_file='protein/1GYZ/aat000_03_05.200_v1_3',protein_directory='protein/1GYZ/1GYZ.clean1.pdb',\
                 short_inserts = 1,cycles = 10,iteration = 100,PDB_out = False, figure_num = 1):
        '''This function's difference from FELTR_new_FA is used usual way to initial the protein
        
        what's difference between this and FELTR_new is add temperature layer
        
        parameters:
        file_directory: stand for the folder's place which contain sub-folder protein et.
        KT: MC parameter
        protein_name: the protein's name
        protein_sequence: the given protein's sequence
        fragment_file: the given protein's fragment file path
        protein_directory: the given protein's pdb file 's path
        short_inserts: the FA run time
        cycles: the of MC run times
        iteration: the ensemble number
        PDB_out: whether output pdb file

        #Author: QCQ
        
        modified in 2014/3/12
        add temperature layer
        '''
        #xaxes=[] #存储以后画图所需的数据
        #yaxes=[]
        nsel = 0  #改能量级的房子被选中的频率
        nconfs = 0 #选中的能量级的cell中构象的个数
        counter = 1
        com_pose = pose_from_pdb(file_directory + protein_directory)
        #output all information into a folder called "protein_name"
        
        file_directory_information = file_directory + str(run_time) + '/'
        #if os.path.exists(file_directory_information):   #create the directory to contain pose
            #shutil.rmtree(file_directory_information)
        
        if not os.path.exists(file_directory_information):
            os.mkdir(file_directory_information)
        
        #add
        if os.path.exists(file_directory_information + protein_name + '/'):
            shutil.rmtree(file_directory_information + protein_name + '/')
        os.mkdir(file_directory_information + protein_name + '/')
        file_directory_information = file_directory_information + protein_name + '/'

        pose = pose_from_sequence(protein_sequence,"fa_standard") 
        #temp_pose =Pose()
        to_centroid = SwitchResidueTypeSetMover("centroid")   #质心switch
        to_fullatom = SwitchResidueTypeSetMover("fa_standard")   #fullatom switch
        scorefxn = create_score_function("score3")      #score3设置
        fragset = ConstantLengthFragSet(3)
        fragset.read_fragment_file(file_directory + fragment_file)

        to_centroid.apply(pose)  #all-atoms-->backbone
        to_centroid.apply(com_pose)  #all-atoms-->backbone

        #KT=1.0
        #mc = MonteCarlo(pose, scorefxn, KT)    #MC
        
        #define the temperature list(add)
        #KT_list = [i * 0.2 for i in range(1,2)]
        #KT_number = len(KT_list)
        length_of_protein = len(protein_sequence)
        #define a list of MC object for different layers in ensemble(add)
        mc_list = [MonteCarlo(com_pose, scorefxn, i) for i in KT_list]
        
        #mc = MonteCarlo(pose, scorefxn, KT)
        movemap = MoveMap()
        movemap.set_bb(True)
        #fragset = ConstantLengthFragSet(3)
        #fragset.read_fragment_file(file_directory + fragment_file)
        mover_3mer = ClassicFragmentMover(fragset, movemap)   #片段库
        '''
        insert_short_frag = RepeatMover(mover_3mer, short_inserts)   #apply short_inserts times 3-mer fragment assmebly
        folding_mover = SequenceMover()  #apply a sequence perturbation to protin conformation
        folding_mover.add_mover(insert_short_frag)
    
        folding_list = []
        for index in mc_list:   
            trial = TrialMover(folding_mover, index)
            #folding_list.append(RepeatMover(trial, cycles))  #定义进行多少次MC
                                        #这里有一个奇怪现象，如果把RepeatMover整体包装进列表
                                      #里面，下边执行apply方法时会提示出错，只能封装到TrivalMover
                                      #这一层。
            
            folding_list.append(trial)
            
        '''    
        #folding = RepeatMover(trial, cycles)  #定义进行多少次MC
        #folding_mover.add_mover    #this object can add other mover
        
        #mc = MonteCarlo(pose, scorefxn, 1.0)
        #trial = TrialMover(folding_mover, mc)
        #folding = RepeatMover(trial, cycles)  #定义进行多少次MC
        #random.seed(1)

        score_max = scorefxn(pose)   #获取一个可能最高的能量值，用于给定数据结构的上限
        
        #define a list to get temperature
        temperature_list = [score_max / len(KT_list) * 1.0 * i for i in range(1,len(KT_list) + 1)]
        temperature_list.insert(0,0.0)
        energy_num = int(math.fabs(score_max)/2) + 1
        for i in range(1,pose.total_residue()+1):
                pose.set_phi(i,120)    #设置蛋白质构象的初始角度
                pose.set_psi(i,-120)
                pose.set_omega(i,180)
        #数据结构定义：用一个列表嵌套所有的代表特定能量水平的字典数据结构，在字典数据结构中-
        #-用其键值代表其cell类型（整型值），字典的值部分是一个列表，把属于这个cell的构象放入此列表-
        #-列表的第一个元素用于记载此cell被引用的次数
        pose_list = [{} for i in range(0,energy_num)]   #用于存储数据的数据结构//增加能量水平

        pose_score = scorefxn(pose)  #测定pose的能量函数
        #yaxes.append(pose_score)
        #xaxes.append(CA_rmsd(com_pose,pose))
        #pose_score = 101.9  #这是一个测试使用的能量值
        usr = UltrafastShapeRecognition(pose)    #创建求取USR的对象
        vector_temp = usr.low_dimensional_geometric_projection()   #求取三维的USR向量
        for i in vector_temp:   #测试使用，验证是否求得对应的USR值
            print i,

        current_temp_pose = Pose()
        current_temp_pose.assign(pose)
        #构建一个包含能量值，蛋白质构象和usr值的对象
        pose_conformation_current = PoseConformation(pose_score, current_temp_pose, vector_temp, CA_rmsd(com_pose,pose))
        
        #acquire the temperature of the current protein
        for i,j in zip(range(len(temperature_list) - 1), KT_list):
            if pose_score > temperature_list[i] and pose_score <= temperature_list[i + 1]:
                temperature = j
                break
                    
        pose_conformation_current.set_temperature(temperature)
        #将当前构建的对象放入上边定义的数据结构的对应位置
        pose_list[energy_real_level_to_index(pose_score)][usr_to_cell(vector_temp)] = [0]#引用次数设置为0
        pose_list[energy_real_level_to_index(pose_score)][usr_to_cell(vector_temp)].append(pose_conformation_current)


        #开始迭代
        time_all_start = time.clock()  #用于计算整个迭代过程所需要的时间
        try:
            test_file = open(file_directory_information + "test.txt","w")   #output test file
            
            temperature_file = open(file_directory_information + "temperature.txt","w")
                
            for iter_time in range(1,iteration):
                time_start = time.clock()  #用于计算一次迭代所需要的时间
                #pdb.set_trace()
                energy_level_list = []  #根据相应的权重函数构建一个概率列表
                energy_level_list_weight=[]
                temp1 = []
                energy_level_chance = []
                #只是用当前能量级里面cell内有构象的能量级组成一个概率列表，防止随机数产生在根本就没有构象的能量级上
                for j,x in zip(pose_list,range(0,len(pose_list))):
                    if j:   #判断某一个能量级是不是有cell：也就是判断是不是能量级所在的字典为空
                        #print "yes",
                        energy_level_list.append(index_to_energy_level(x))
                    else:
                        #print "no",
                        pass
                            
                #print energy_level_list
                """概率列表的构建具体方法是：
                
                                            根据现在数据结构里面有的能量级加上其权重构建一个概率列表
                                            类似[0.0,……,1.0]这样的概率列表，同时用产生的随机数看其落在那个相应的区间位置，即为当前选定的能量水平
                """
                energy_level_list_weight = [y*y for y in energy_level_list]  #权重计算
                energy_level_list_weight.reverse()
                for x,y in zip(energy_level_list_weight, range(0,len(energy_level_list_weight))):
                    temp1.append(x * 1.0 / sum(energy_level_list_weight))
                    energy_level_chance.append(sum(temp1[0:y+1]))
                energy_level_chance.insert(0,0.0)
                #选定能量水平
                level = 0
                test_file.write(str(iter_time)+" ")
                rand = random.random()
                for index1,index2 in zip(range(0,len(energy_level_chance)-1),energy_level_list):
                    if rand > energy_level_chance[index1] and rand <= energy_level_chance[index1+1]:
                        level = index2
                        break
                
                #print level
                test_file.write("随机数："+str(rand)+" "+"energy:"+str(level)+" ")
                #选定cell的权重和概率分配
                cell_list = []
                cell_list_weight = []
                temp2 = []
                cell_list_chance = []
                index_of_level = energy_level_to_index(level)
                nsel = 0
                for pose_list_keys in pose_list[index_of_level].keys():  #求的所有cell被索引的次数
                    nsel = nsel + pose_list[index_of_level][pose_list_keys][0]
                        
                for pose_list_keys in pose_list[index_of_level].keys():
                    #print pose_list_keys
                    cell_list.append(pose_list_keys)
                    nsel_current = pose_list[index_of_level][pose_list_keys][0]  #当前cell的索引次数
                    nconfs = len(pose_list[index_of_level][pose_list_keys]) - 1 #当前cell里面构象的数目
                    #print nconfs
                    #print nconfs
                    if(nsel == 0):
                        nsel = 1
                    cell_list_weight.append(1.0/((1.0 + nsel_current * 1.0 / nsel) * nconfs))
                        
                       
                for x,y in zip(cell_list_weight, range(0,len(cell_list_weight))):
                    temp2.append(x * 1.0 / sum(cell_list_weight))
                    cell_list_chance.append(sum(temp2[0:y+1]))
                cell_list_chance.insert(0,0.0)
                #选定cell
                cell = 0
                rand1 = random.random()
                
                for index3,index4 in zip(range(0,len(cell_list_chance) - 1),cell_list):
                    if rand1 > cell_list_chance[index3] and rand1 <= cell_list_chance[index3+1]:
                        cell = index4
                        break
                #print cell
                test_file.write("随机数："+str(rand1)+" "+"cell:"+str(cell)+" ")
                #选定特定level里面特定cell里面的conformation
                rand2 = random.randint(1,len(pose_list[index_of_level][cell]) - 1 )
                #rand2 = 1
                test_file.write("随机数："+str(rand2)+" "+str(len(pose_list[index_of_level][cell]) - 1 )+'\n')
                pose_conformation_current = pose_list[index_of_level][cell][rand2]  #取得这个要作变化的构象
                #print pose_conformation_current
                #print pose_conformation_current.get_protein()
                pose = pose_conformation_current.get_protein()
                
                #get the protein's temperature(add)
                before_temperature = pose_conformation_current.get_temperature()                    
                temperature_file.write('Before the change, the temperature is: ' + str(before_temperature) + " ")
                
                pose2 = Pose()    #由于是指针传递，防止改变原来存储起来的对象，进行新建对象
                pose2.assign(pose)
                #pose_conformation_current.get_energy()
                #pose_conformation_current.get_energy()
                #according the temperature to get the index of the corresponding MC object and TrialMover below
                mc_list[KT_list.index(before_temperature)].reset(pose2)  #reset the MonteCarlo object
                
                #added by qcq
                short_inserts_add = short_inserts
                insert_short_frag = RepeatMover(mover_3mer, short_inserts_add)   #apply short_inserts times 3-mer fragment assmebly
                folding_mover = SequenceMover()  #apply a sequence perturbation to protin conformation
                folding_mover.add_mover(insert_short_frag) 
                trial = TrialMover(folding_mover, mc_list[KT_list.index(before_temperature)])
                cycles_add = cycles   
                RepeatMover(trial,cycles_add).apply(pose2)
                #recover the lowest scoring decoy structure to pose2
                mc_list[KT_list.index(before_temperature)].recover_low(pose2)
                #get current energy
                score_current = scorefxn(pose2) 
                '''
                if (pose_conformation_current.get_energy() != score_current ):#and score_current <= 199.9):  #采用判断当前的能量值和上一个是否相同来判断是否MC接受了当前的蛋白质构象
                    counter = counter + 1    #用于计数被接受了的蛋白质构象数目，相当于计算当前数据结构里面共有多少构象的计数
                    pose_list[index_of_level][cell][0] = pose_list[index_of_level][cell][0] + 1  #增加cell的引用次数
                    #pose1 = Pose()
                    #pose1.assign(pose2)
                    usr.set_protein(pose2)
                    #yaxes.append(score_current)
                    #xaxes.append(CA_rmsd(com_pose,pose2))
                    #usr = UltrafastShapeRecognition(pose2)
                    usr_vector = usr.low_dimensional_geometric_projection()
                    
                    #calculate the new protein conformation's temperature
                    pose_conformation_current = PoseConformation(score_current, pose2, usr_vector, CA_rmsd(com_pose, pose2))
                    for i,j in zip(range(len(temperature_list) - 1), KT_list):
                        if score_current > temperature_list[i] and score_current <= temperature_list[i + 1]:
                            after_temperature = j
                            break
                        
                    pose_conformation_current.set_temperature(after_temperature)
                    
                    temperature_file.write('After the change, the temperature is: ' + str(after_temperature) + "\n")
                    
                    #judge the temperature changed
                    if before_temperature != after_temperature:
                        pose_conformation_temp = pose_list[index_of_level][cell][rand2]
                        pose_list[index_of_level][cell][rand2] = pose_conformation_current
                        if not pose_list[energy_real_level_to_index(score_current)].has_key(usr_to_cell(usr_vector)):  #如果当前的数据结构里面没有这个cell类型就新建一个此种类型的键值对
                            pose_list[energy_real_level_to_index(score_current)][usr_to_cell(usr_vector)] = [0]  #这个cell的引用次数
                            pose_list[energy_real_level_to_index(score_current)][usr_to_cell(usr_vector)].append(pose_conformation_temp)
                        else:
                            pose_list[energy_real_level_to_index(score_current)][usr_to_cell(usr_vector)].append(pose_conformation_temp)
                    else:
                    
                        #pose_conformation_current = PoseConformation(score_current, pose2, usr_vector, CA_rmsd(com_pose,pose2))
                        #pose_conformation_current.set_all(score_current, pose2, usr_vector)
                        print pose_conformation_current.get_usr()[0],pose_conformation_current.get_usr()[1],pose_conformation_current.get_usr()[2]
                        if not pose_list[energy_real_level_to_index(score_current)].has_key(usr_to_cell(usr_vector)):  #如果当前的数据结构里面没有这个cell类型就新建一个此种类型的键值对
                            pose_list[energy_real_level_to_index(score_current)][usr_to_cell(usr_vector)] = [0]  #这个cell的引用次数
                            pose_list[energy_real_level_to_index(score_current)][usr_to_cell(usr_vector)].append(pose_conformation_current)
                        else:
                            pose_list[energy_real_level_to_index(score_current)][usr_to_cell(usr_vector)].append(pose_conformation_current)
                    
                '''    
                counter = counter + 1    #用于计数被接受了的蛋白质构象数目，相当于计算当前数据结构里面共有多少构象的计数
                pose_list[index_of_level][cell][0] = pose_list[index_of_level][cell][0] + 1  #增加cell的引用次数
                #pose1 = Pose()
                #pose1.assign(pose2)
                usr.set_protein(pose2)
                #yaxes.append(score_current)
                #xaxes.append(CA_rmsd(com_pose,pose2))
                #usr = UltrafastShapeRecognition(pose2)
                usr_vector = usr.low_dimensional_geometric_projection()
                
                #calculate the new protein conformation's temperature
                pose_conformation_current = PoseConformation(score_current, pose2, usr_vector, CA_rmsd(com_pose, pose2))
                for i,j in zip(range(len(temperature_list) - 1), KT_list):
                    if score_current > temperature_list[i] and score_current <= temperature_list[i + 1]:
                        after_temperature = j
                        break
                    
                pose_conformation_current.set_temperature(after_temperature)
                
                temperature_file.write('After the change, the temperature is: ' + str(after_temperature) + "\n")
                
                #judge the temperature changed
                '''
                if before_temperature != after_temperature:
                    pose_conformation_temp = pose_list[index_of_level][cell][rand2]
                    pose_list[index_of_level][cell][rand2] = pose_conformation_current
                    if not pose_list[energy_real_level_to_index(score_current)].has_key(usr_to_cell(usr_vector)):  #如果当前的数据结构里面没有这个cell类型就新建一个此种类型的键值对
                        pose_list[energy_real_level_to_index(score_current)][usr_to_cell(usr_vector)] = [0]  #这个cell的引用次数
                        pose_list[energy_real_level_to_index(score_current)][usr_to_cell(usr_vector)].append(pose_conformation_temp)
                    else:
                        pose_list[energy_real_level_to_index(score_current)][usr_to_cell(usr_vector)].append(pose_conformation_temp)
                else:    
                '''    
                     
                #pose_conformation_current.set_all(score_current, pose2, usr_vector)
                print pose_conformation_current.get_usr()[0],pose_conformation_current.get_usr()[1],pose_conformation_current.get_usr()[2],pose_conformation_current.get_energy(),rand,rand1,rand2
                if not pose_list[energy_real_level_to_index(score_current)].has_key(usr_to_cell(usr_vector)):  #如果当前的数据结构里面没有这个cell类型就新建一个此种类型的键值对
                    pose_list[energy_real_level_to_index(score_current)][usr_to_cell(usr_vector)] = [0]  #这个cell的引用次数
                    pose_list[energy_real_level_to_index(score_current)][usr_to_cell(usr_vector)].append(pose_conformation_current)
                else:
                    pose_list[energy_real_level_to_index(score_current)][usr_to_cell(usr_vector)].append(pose_conformation_current)
                  
                time_end = time.clock()
                print "The %d iteration of %s takes %f seconds"%(iter_time,protein_name,time_end - time_start)
        finally:
            test_file.close()
            #(add)
            temperature_file.close()
            
        time_all_end = time.clock()
        print "The %s function takes %f seconds"%(protein_name,time_all_end - time_all_start)

        figure(num = figure_num)
        #xaxes=[]
        #yaxes=[]
        #将信息输出到information.txt文件，同时作图
        try:
            f = open(file_directory_information + "information.txt","w")
            fdata = open(file_directory_information + "data.txt","w")
            #fdata1 = open(file_directory_information + "data1.txt","w")
            f.write('\t\t%s\n'%(protein_name))
            fdata.write('\t\t%s\n'%(protein_name))
            pdb_directory_information = file_directory_information +'pdb/'
            for i,j in zip(pose_list,range(0, len(pose_list))):
                flag1 = False
                if not flag1:
                    f.write("the current energy level is: %d."%(index_to_energy_level(j)))
                    f.write("\n")
                    flag1 = True
                    
                
                for x,y in i.items():
                    flag = False
                    if not flag:
                        f.write("\tthis cell is :%d ,referenced %d times."%(x,y[0]))
                        f.write("\n")
                        flag = True
                    for index in range(1,len(y)):   
                        f.write("\t\tThe USR is:(%f,%f,%f) score:%f"%(y[index].get_usr()[0],y[index].get_usr()[1],y[index].get_usr()[2],y[index].get_energy()))
                        f.write("\n")
                        #xaxes.append(CA_rmsd(com_pose,y[index].get_protein()))
                        #yaxes.append(y[index].get_energy())
                        scatter(y[index].get_CA_rmsd(), y[index].get_energy(), 1, color ='red')  #作点图
                        #output indormation
                        fdata.write(str(y[index].get_CA_rmsd()))
                        fdata.write('\t')
                        fdata.write(str(y[index].get_energy()))
                        fdata.write('\n')
                
                        if PDB_out:
                                
                                if not os.path.exists(pdb_directory_information):
                                        os.mkdir(pdb_directory_information)
                                if float(y[index].get_CA_rmsd()) < 4:
                                        pdb_file_name = pdb_directory_information + \
                                                        str(y[index].get_CA_rmsd()) + '_' + \
                                                        str(y[index].get_usr()[0]) + '_' + \
                                                        str(y[index].get_usr()[1]) + '_' + \
                                                        str(y[index].get_usr()[2]) + '_' + \
                                                        str(y[index].get_energy()) + '.pdb'
                                                
                                        y[index].get_protein().dump_pdb(pdb_file_name)
                        #annotate(str(index+1),xy=(CA_rmsd(com_pose,y[index].get_protein()), y[index].get_energy()), xycoords='data',size='xx-small') #注解
            f.write("The %s function takes %f seconds and have %d conforamtions"%(protein_name,time_all_end - time_all_start, counter))
            '''
            for index_x,index_y in zip(xaxes,yaxes):
                fdata1.write(str(index_x))
                fdata1.write('\t')
                fdata1.write(str(index_y))
                fdata1.write('\n')
            '''
                
            ax = gca()
            ax.spines['right'].set_color('none')
            ax.spines['top'].set_color('none')
            ax.xaxis.set_ticks_position('bottom')
            ax.spines['bottom'].set_position(('data',0))
            ax.yaxis.set_ticks_position('left')
            ax.spines['left'].set_position(('data',0))
            xticks([dexi*1 for dexi in range(0,20)])
            yticks([dexi*1 for dexi in range(0,energy_num * 2,25)])
            xlabel('CA_rmsd')
            #yab=str(energy_num)
            ylabel('energy')
            title('The comparison between %s ensemble and standard protein'%(protein_name),fontsize=15,color='red')
            grid(True)
            xlim(0.0,20.0)
            ylim(0.0,energy_num * 2)
            savefig(file_directory_information + 'hello1.png',dpi=90)        
                    
                                                   
        finally:
            fdata.close()
            f.close()
            #break


parser = optparse.OptionParser()
parser.add_option('--file_directory', dest = 'file_directory',
    default = '',    # default example PDB
    help = 'The project\'s path')
parser.add_option('--kT', dest='kT',
    default = '1.0',
    help = 'the \"temperature\" of the sample refinement protocol')
parser.add_option( '--protein_name', dest='protein_name',
    default = '1GYZ',
    help = 'the name of the protein' )
parser.add_option('--protein_sequence', dest='protein_sequence',
    default = 'WIARINAAVRAYGLNYSTFINGLKKAGIELDRKILADMAVRDPQAFEQVVNKVKEALQVQ',
    help = 'The sequence of the protein' )
parser.add_option( '--fragment_file', dest='fragment_file',
    default = 'protein/1GYZ/aat000_03_05.200_v1_3',
    help = 'the fragment file')
parser.add_option('--cycles', dest='cycles',
    default = '10',
    help = 'the number of refinement rounds (small, shear, min, pack) in\
        the sample refinement protocol')
parser.add_option('--protein_directory', dest='protein_directory',
    default = 'protein/1GYZ/1GYZ.clean1.pdb',    # default to single trajectory for speed
    help = 'the place of the protein')

parser.add_option('--short_inserts', dest = 'short_inserts',
    default = '1',    # if a specific output name is desired
    help = 'the fragment assmbly applied times')
parser.add_option('--iteration', dest = 'iteration',
    default = '100',    # if a specific output name is desired
    help = 'the conformation generation nums')
parser.add_option('--PDB_out', dest = 'PDB_out',
    default = 'False',    # if a specific output name is desired
    help = 'if output pdb')
(options,args) = parser.parse_args()

KT = float(options.kT)

short_inserts = int(options.short_inserts) #define times applied 3-mer fragment assembly
cycles = int(options.cycles)  ##定义进行多少次MC
iteration = int(options.iteration) #define the num of the conformations
file_directory_linux='/home/qcq/tree_based_memory/'#linux绝对目录
file_directory_windows = ''#windows相对目录（当前目录）
file_directory = file_directory_windows  #use windows
protein_directory=options.protein_directory
protein_name = options.protein_name
protein_sequence = options.protein_sequence
fragment_file = options.fragment_file
PDB_out = bool(options.PDB_out)

"""
fold_protein(file_directory,KT,protein_name,protein_sequence,fragment_file,protein_directory,\
             short_inserts,cycles,iteration,PDB_out)

"""

















    
    
