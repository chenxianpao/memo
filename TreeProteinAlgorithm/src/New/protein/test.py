#!/usr/bin/env python
# -*- coding: utf-8 -*-
import random
from rosetta import *
init()
import os
from cleanfile import *
protein_1ENH_output='1ENH/1ENH.clean.pdb'
protein_1ENH_clean='1ENH/1ENH.clean1.pdb'
protien_1ENH_input = '1ENH/1ENH.pdb'
sequence_1ENH = 'RPRTAFSSEQLARLKREFNENRYLTERRRQQLSSELGLNEAQIKIWFQNKRAKI'
fragment_file_1ENH='1ENH/aat000_03_05.200_v1_3'
pose_1ENH = pose_from_sequence(sequence_1ENH,"fa_standard")
pose_standard_1ENH = pose_from_pdb(protein_1ENH_clean)

protein_1GB1_output='1GB1/1GB1.clean.pdb'
protein_1GB1_clean='1GB1/1GB1.clean1.pdb'
protien_1GB1_input = '1GB1/1GB1.pdb'
sequence_1GB1 = 'MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE'
fragment_file_1GB1='1GB1/aat000_03_05.200_v1_3'
pose_1GB1 = pose_from_sequence(sequence_1GB1,"fa_standard")
pose_standard_1GB1 = pose_from_pdb(protein_1GB1_clean)

protein_1GYZ_output='1GYZ/1GYZ.clean.pdb'
protein_1GYZ_clean='1GYZ/1GYZ.clean1.pdb'
protien_1GYZ_input = '1GYZ/1GYZ.pdb'
sequence_1GYZ = 'WIARINAAVRAYGLNYSTFINGLKKAGIELDRKILADMAVRDPQAFEQVVNKVKEALQVQ'
fragment_file_1GYZ = '1GYZ/aat000_03_05.200_v1_3'
pose_1GYZ = pose_from_sequence(sequence_1GYZ,"fa_standard")
pose_standard_1GYZ = pose_from_pdb(protein_1GYZ_clean)

protein_1I6C_output='1I6C/1I6C.clean.pdb'
protein_1I6C_clean='1I6C/1I6C.clean1.pdb'
protien_1I6C_input = '1I6C/1I6C.pdb'
sequence_1I6C = 'KLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGNSSSG'
fragment_file_1I6C = '1I6C/aat000_03_05.200_v1_3'
pose_1I6C = pose_from_sequence(sequence_1I6C,"fa_standard")
pose_standard_1I6C = pose_from_pdb(protein_1I6C_clean)

protein_1VII_output='1VII/1VII.clean.pdb'
protein_1VII_clean='1VII/1VII.clean1.pdb'
protien_1VII_input = '1VII/1VII.pdb'
sequence_1VII = 'MLSDEDFKAVFGMTRSAFANLPLWKQQNLKKEKGLF'
fragment_file_1VII='1VII/aat000_03_05.200_v1_3'
pose_1VII = pose_from_sequence(sequence_1VII,"fa_standard")
pose_standard_1VII = pose_from_pdb(protein_1VII_clean)

protein_1L2Y_output='1L2Y/1L2Y.clean.pdb'   #too short no fragment library
protein_1L2Y_clean='1L2Y/1L2Y.clean1.pdb'
protien_1L2Y_input = '1L2Y/1L2Y.pdb'
sequence_1L2Y = 'NLYIQWLKDGGPSSGRPPPS'
fragment_file_1L2Y ='1L2Y/aat000_03_05.200_v1_3'
pose_1L2Y = pose_from_sequence(sequence_1L2Y,"fa_standard")
pose_standard_1L2Y = pose_from_pdb(protein_1L2Y_clean)

protein_4ICB_output='4ICB/4ICB.clean.pdb'
protein_4ICB_clean='4ICB/4ICB.clean1.pdb'
protien_4ICB_input = '4ICB/4ICB.pdb'
sequence_4ICB = 'MKSPEELKGIFEKYAAKEGDPNQLSKEELKLLLQTEFPSLLKGPSTLDELFEELDKNGDGEVSFEEFQVLVKKISQ'
fragment_file_4ICB='4ICB/aat000_03_05.200_v1_3'
pose_4ICB = pose_from_sequence(sequence_4ICB,"fa_standard")
pose_standard_4ICB = pose_from_pdb(protein_4ICB_clean)

input_file = []
output_file = []
pose=[]
pose_standard=[]
input_file.append(protien_1ENH_input)
input_file.append(protien_1GB1_input)
input_file.append(protien_1GYZ_input)
input_file.append(protien_1I6C_input)
input_file.append(protien_1VII_input)
#input_file.append(protien_1L2Y_input)
input_file.append(protien_4ICB_input)

output_file.append(protein_1ENH_output)
output_file.append(protein_1GB1_output)
output_file.append(protein_1GYZ_output)
output_file.append(protein_1I6C_output)
output_file.append(protein_1VII_output)
#output_file.append(protein_1L2Y_output)
output_file.append(protein_4ICB_output)

pose.append(pose_1ENH)
pose.append(pose_1GB1)
pose.append(pose_1GYZ)
pose.append(pose_1I6C)
pose.append(pose_1VII)
pose.append(pose_1L2Y)
pose.append(pose_4ICB)

pose_standard.append(pose_standard_1ENH)
pose_standard.append(pose_standard_1GB1)
pose_standard.append(pose_standard_1GYZ)
pose_standard.append(pose_standard_1I6C)
pose_standard.append(pose_standard_1VII)
pose_standard.append(pose_standard_1L2Y)
pose_standard.append(pose_standard_4ICB)

'''
for index_x,index_y in zip(input_file,output_file):
    clean_file(index_x,index_y)
'''    
to_centroid = SwitchResidueTypeSetMover("centroid")   #质心switch
to_fullatom = SwitchResidueTypeSetMover("fa_standard")   #fullatom switch
scorefxn=create_score_function("score3")      #score3设置
'''
for i in range(1,pose.total_residue()+1):
	pose.set_psi(i,120)    #设置蛋白质构象的初始角度
	pose.set_phi(i,-120)
	pose.set_omega(i,180)
'''
for j,standard in zip(pose,pose_standard):
    to_centroid.apply(j)
    score = scorefxn(j)
    to_centroid.apply(standard)
    for i in range(1,j.total_residue()+1):
        j.set_psi(i,120)
        j.set_phi(i,-120)
        j.set_omega(i,180)
    content = str(scorefxn(j))+'\t'+str(scorefxn(standard))
    print score,content
    
