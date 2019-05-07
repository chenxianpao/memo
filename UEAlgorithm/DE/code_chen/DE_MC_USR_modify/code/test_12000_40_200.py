#!/usr/bin/env python
# -*- coding: utf-8 -*-
import random
from rosetta import *
#init()
import os
from FELTR_new import *
protein_1ENH_clean='protein/1ENH/1ENH.clean1.pdb'
protien_1ENH_name = '1ENH'
sequence_1ENH = 'RPRTAFSSEQLARLKREFNENRYLTERRRQQLSSELGLNEAQIKIWFQNKRAKI'
fragment_file_1ENH='protein/1ENH/aat000_03_05.200_v1_3'

protein_1GB1_clean='protein/1GB1/1GB1.clean1.pdb'
protien_1GB1_name = '1GB1'
sequence_1GB1 = 'MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE'
fragment_file_1GB1='protein/1GB1/aat000_03_05.200_v1_3'

protein_1GYZ_clean='protein/1GYZ/1GYZ.clean1.pdb'
protien_1GYZ_name = '1GYZ'
sequence_1GYZ = 'WIARINAAVRAYGLNYSTFINGLKKAGIELDRKILADMAVRDPQAFEQVVNKVKEALQVQ'
fragment_file_1GYZ = 'protein/1GYZ/aat000_03_05.200_v1_3'

protein_1I6C_clean='protein/1I6C/1I6C.clean1.pdb'
protien_1I6C_name = '1I6C'
sequence_1I6C = 'KLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGNSSSG'
fragment_file_1I6C = 'protein/1I6C/aat000_03_05.200_v1_3'

protein_1VII_clean='protein/1VII/1VII.clean1.pdb'
protien_1VII_name = '1VII'
sequence_1VII = 'MLSDEDFKAVFGMTRSAFANLPLWKQQNLKKEKGLF'
fragment_file_1VII='protein/1VII/aat000_03_05.200_v1_3'

protein_1L2Y_output='protein/1L2Y/1L2Y.clean.pdb'   #too short no fragment library
protein_1L2Y_clean='1L2Y/1L2Y.clean1.pdb'
protien_1L2Y_input = '1L2Y/1L2Y.pdb'
sequence_1L2Y = 'NLYIQWLKDGGPSSGRPPPS'
fragment_file_1L2Y ='1L2Y/aat000_03_05.200_v1_3'


protein_4ICB_clean='protein/4ICB/4ICB.clean1.pdb'
protien_4ICB_name = '4ICB'
sequence_4ICB = 'MKSPEELKGIFEKYAAKEGDPNQLSKEELKLLLQTEFPSLLKGPSTLDELFEELDKNGDGEVSFEEFQVLVKKISQ'
fragment_file_4ICB='protein/4ICB/aat000_03_05.200_v1_3'
#cxp's protein start
protein_1BBO_clean='protein/1BBO/1BBO.clean1.pdb'
protein_1BBO_name = '1BBO'
sequence_1BBO = 'KYICEECGIRAKKPSMLKKHIRTHTDVRPYHCTYCNFSFKTKGNLTKHMKSKAHSKK'
fragment_file_1BBO='protein/1BBO/aat000_03_05.200_v1_3'

protein_1MN3_clean='protein/1MN3/1MN3.clean1.pdb'
protein_1MN3_name = '1MN3'
sequence_1MN3 = 'SSLIKKIEENERKDTLNTLQNMFPDMDPSLIEDVCIAKKSRIEPCVDALLSLSE'
fragment_file_1MN3='protein/1MN3/aat000_03_05.200_v1_3'

protein_2IMU_clean='protein/2IMU/2IMU.clean1.pdb'
protein_2IMU_name = '2IMU'
sequence_2IMU = 'FGFKDIIRAIRRIAVPVVSTLFPPAAPLAHAIGEGVDYLLGDEAQA'
fragment_file_2IMU='protein/2IMU/aat000_03_05.200_v1_3'

protein_2JUJ_clean='protein/2JUJ/2JUJ.clean1.pdb'
protein_2JUJ_name = '2JUJ'
sequence_2JUJ = 'ATASPQLSSEIENLMSQGYSYQDIQKALVIAQNNIEMAKNILREFVSISSPAHVAT'
fragment_file_2JUJ='protein/2JUJ/aat000_03_05.200_v1_3'

protein_1FD4_clean='protein/1FD4/1FD4.clean1.pdb'
protein_1FD4_name = '1FD4'
sequence_1FD4 = 'GIGDPVTCLKSGAICHPVFCPRRYKQIGTCGLPGTKCCKKP'
fragment_file_1FD4='protein/1FD4/aat000_03_05.200_v1_3'
#end cxp
#long protein
protein_3GWL_clean='protein/3GWL/3GWL.clean1.pdb'
protein_3GWL_name = '3GWL'
sequence_3GWL = 'GSHMLHWGPKYWRSLHLYAIFFSDAPSWKEKYEAIQWILNFIESLPCTRCQHHAFSYLTKNPLTLNNSEDFQYWTFAFHNNVNNRLNKKIISWSEYKNIYEQSILK'
fragment_file_3GWL='protein/3GWL/aat000_03_05.200_v1_3'


protein_clean_file_name = []
protein_clean_file_name.append(protein_1ENH_clean)
protein_clean_file_name.append(protein_1GB1_clean)
protein_clean_file_name.append(protein_1GYZ_clean)
protein_clean_file_name.append(protein_1I6C_clean)
protein_clean_file_name.append(protein_1VII_clean)
protein_clean_file_name.append(protein_4ICB_clean)
protein_clean_file_name.append(protein_1BBO_clean)
protein_clean_file_name.append(protein_1MN3_clean)
protein_clean_file_name.append(protein_2IMU_clean)
protein_clean_file_name.append(protein_2JUJ_clean)
protein_clean_file_name.append(protein_3GWL_clean)
protein_clean_file_name.append(protein_1FD4_clean)

protein_name = []
protein_name.append(protien_1ENH_name)
protein_name.append(protien_1GB1_name)
protein_name.append(protien_1GYZ_name)
protein_name.append(protien_1I6C_name)
protein_name.append(protien_1VII_name)
protein_name.append(protien_4ICB_name)
protein_name.append(protein_1BBO_name)
protein_name.append(protein_1MN3_name)
protein_name.append(protein_2IMU_name)
protein_name.append(protein_2JUJ_name)
protein_name.append(protein_3GWL_name)
protein_name.append(protein_1FD4_name)


sequence_name = []
sequence_name.append(sequence_1ENH)
sequence_name.append(sequence_1GB1)
sequence_name.append(sequence_1GYZ)
sequence_name.append(sequence_1I6C)
sequence_name.append(sequence_1VII)
sequence_name.append(sequence_4ICB)
sequence_name.append(sequence_1BBO)
sequence_name.append(sequence_1MN3)
sequence_name.append(sequence_2IMU)
sequence_name.append(sequence_2JUJ)
sequence_name.append(sequence_3GWL)
sequence_name.append(sequence_1FD4)

fragment_name=[]
fragment_name.append(fragment_file_1ENH)
fragment_name.append(fragment_file_1GB1)
fragment_name.append(fragment_file_1GYZ)
fragment_name.append(fragment_file_1I6C)
fragment_name.append(fragment_file_1VII)
fragment_name.append(fragment_file_4ICB)
fragment_name.append(fragment_file_1BBO)
fragment_name.append(fragment_file_1MN3)
fragment_name.append(fragment_file_2IMU)
fragment_name.append(fragment_file_2JUJ)
fragment_name.append(fragment_file_3GWL)
fragment_name.append(fragment_file_1FD4)

file_directory_path = '/home/qcq/tree_based_memory/'
iterate = 2500



want_to_run = ['1ENH','1GB1','1GYZ','1I6C','1VII','4ICB',
               '1BBO','1MN3','2IMU','2JUJ'#,'1FD4','3GWL'\
               ]

want_to_run1 = ['2JUJ']

for i,j,k,z,v in zip(protein_name,sequence_name,fragment_name,\
                 protein_clean_file_name, range(1,len(want_to_run) + 1)):
    if i in want_to_run:
	mc_number = 50 #MC times is sequence - 2
	short_insert = 10
	print len(j),mc_number,short_insert
        fold_protein(file_directory='',protein_name = i, protein_sequence = j, fragment_file = k,\
                     iteration = iterate,protein_directory= z,short_inserts=short_insert\
                     ,cycles = mc_number,PDB_out = True,figure_num = v)
'''
file_directory_path1 = '/home/qcq/tree_based_memory/'
iterate1 = 12000
short_insert1 = 20
mc_number1 = 100	

want_to_run1 = [#'1ENH','1GB1','1GYZ','1I6C','1VII','4ICB',
               #'1BBO','1MN3','2IMU','2JUJ',
'3GWL'
]

for i,j,k,z in zip(protein_name,sequence_name,fragment_name,\
                 protein_clean_file_name):
    if i in want_to_run1:
        fold_protein(file_directory=file_directory_path1,protein_name = i, protein_sequence = j, fragment_file = k,\
                     iteration = iterate1,protein_directory= z,short_inserts=short_insert1\
                     ,cycles = mc_number1)
    
'''


    


