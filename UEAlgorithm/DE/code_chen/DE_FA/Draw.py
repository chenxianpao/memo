#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The code in this  file used to draw the picture
"""
from rosetta import *
from pylab import *  #画图所需
import time
def drawing_picture(file_in,file_out,figure_num,max_rmsd,max_energy,protein_pdb_file_path):
    figure(num = figure_num)
    with open(file_in[0],'r') as f0:    
        for line in f0:
            if line[0] == 'T' or line =='':
                continue
            linedata = line.split('\t')
            scatter(linedata[0],linedata[1],1,color='black')
    
    with open(file_in[1],'r') as f1:    
        for line in f1:
            linedata = line.split('\t')
            scatter(linedata[0],linedata[1],1,color='red')
    with open(file_in[2],'r') as f2:
        for line in f2:
            linedata = line.split('\t')
            scatter(linedata[0],linedata[1],1,color='blue')
            

    #xticks([dexi*1 for dexi in range(0,25)])
    ax = gca()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.spines['bottom'].set_position(('data',0))
    ax.yaxis.set_ticks_position('left')
    ax.spines['left'].set_position(('data',0))
    if max_rmsd:
        xticks([dexi*1 for dexi in range(0,max_rmsd)])
        xlabel('CA_rmsd')
        #xticks([dexi*0.1 for dexi in range(0,1)])
    else:
        xticks([dexi*0.1 for dexi in range(0,11)])
        xlabel('usr')
    yticks([dexi*1 for dexi in range(0,max_energy,50)])
    
    #yab=str(energy_num)
    ylabel('energy')
    title('The comparison between %s ensemble and standard protein'%(protein_pdb_file_path.split(".")[0].split('/')[-1]),fontsize=15,color='red')
    grid(True)
    savefig(file_out,dpi=90)

'''
file_in=[]
file_directory_information = '1GYZ/'
max_rmsd = 20
max_energy = 1000
protein_1GYZ_clean='protein/1GYZ/1GYZ.clean1.pdb'
protein_pdb_file_path = protein_1GYZ_clean
file_in.append(file_directory_information + "test.txt")
file_in.append(file_directory_information + "child_data.txt")
file_in.append(file_directory_information + "populations_data.txt")
drawing_picture(file_in,"hello1",1,max_rmsd,max_energy,protein_pdb_file_path)
'''

     
        
