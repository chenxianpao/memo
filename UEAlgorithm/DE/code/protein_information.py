#!/usr/bin/env python
# encoding: utf-8
import random
#from rosetta import *
#init()
import os
'''
strictly speaking, this is not a real program.
What contains in this file is many types of string information\
about protein on which to predict its structure.

created on 2014/3/18
Author:QCQ
'''

protein_1ENH_clean='protein/1ENH/1ENH.clean1.pdb'
protein_1ENH_name = '1ENH'
sequence_1ENH = 'RPRTAFSSEQLARLKREFNENRYLTERRRQQLSSELGLNEAQIKIWFQNKRAKI'
fragment_file_1ENH='protein/1ENH/aat000_03_05.200_v1_3'

protein_1GB1_clean='protein/1GB1/1GB1.clean1.pdb'
protein_1GB1_name = '1GB1'
sequence_1GB1 = 'MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE'
fragment_file_1GB1='protein/1GB1/aat000_03_05.200_v1_3'

protein_1GYZ_clean='protein/1GYZ/1GYZ.clean1.pdb'
protein_1GYZ_name = '1GYZ'
sequence_1GYZ = 'WIARINAAVRAYGLNYSTFINGLKKAGIELDRKILADMAVRDPQAFEQVVNKVKEALQVQ'
fragment_file_1GYZ = 'protein/1GYZ/aat000_03_05.200_v1_3'

protein_1I6C_clean='protein/1I6C/1I6C.clean1.pdb'
protein_1I6C_name = '1I6C'
sequence_1I6C = 'KLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGNSSSG'
fragment_file_1I6C = 'protein/1I6C/aat000_03_05.200_v1_3'

protein_1VII_clean='protein/1VII/1VII.clean1.pdb'
protein_1VII_name = '1VII'
sequence_1VII = 'MLSDEDFKAVFGMTRSAFANLPLWKQQNLKKEKGLF'
fragment_file_1VII='protein/1VII/aat000_03_05.200_v1_3'

protein_1L2Y_output='protein/1L2Y/1L2Y.clean.pdb'   #too short no fragment library
protein_1L2Y_clean='1L2Y/1L2Y.clean1.pdb'
protein_1L2Y_input = '1L2Y/1L2Y.pdb'
sequence_1L2Y = 'NLYIQWLKDGGPSSGRPPPS'
fragment_file_1L2Y ='1L2Y/aat000_03_05.200_v1_3'


protein_4ICB_clean='protein/4ICB/4ICB.clean1.pdb'
protein_4ICB_name = '4ICB'
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

#added by tiger 
protein_1AIL_clean='protein/1AIL/1AIL.clean1.pdb'
protein_1AIL_name = '1AIL'
sequence_1AIL = 'MDSNTVSSFQVDCFLWHVRKQVVDQELGDAPFLDRLRRDQKSLRGRGSTLGLNIEAATHVGKQIVEKILKEED'
fragment_file_1AIL = 'protein/1AIL/aat000_03_05.200_v1_3'

protein_1AOY_clean='protein/1AOY/1AOY.clean1.pdb'
protein_1AOY_name = '1AOY'
sequence_1AOY = 'MRSSAKQEELVKAFKALLKEEKFSSQGEIVAALQEQGFDNINQSKVSRMLTKFGAVRTRNAKMEMVYCLPAELGVPTT'
fragment_file_1AOY = 'protein/1AOY/aat000_03_05.200_v1_3'

protein_1CC5_clean='protein/1CC5/1CC5.clean1.pdb'
protein_1CC5_name = '1CC5'
sequence_1CC5 = 'GGGARSGDDVVAKYCNACHGTGLLNAPKVGDSAAWKTRADAKGGLDGLLAQSLSGLNAMPPKGTCADCSDDELKAAIGKMSGL'
fragment_file_1CC5 = 'protein/1CC5/aat000_03_05.200_v1_3'

protein_1FWP_clean='protein/1FWP/1FWP.clean1.pdb'
protein_1FWP_name = '1FWP'
sequence_1FWP = 'RQLALEAKGETPSAVTRLSVVAKSEPQDEQSRSQSPRRIILSRLKAGEVDLLEEELGHLTTLTDVVKGADSLSAILPGDIAEDDITAVLCFVIEADQITFETVEVSPKISTPPVLKLAAEQAPTGRVEREKTTRIKLGT'
fragment_file_1FWP = 'protein/1FWP/aat000_03_05.200_v1_3'



protein_1SAP_clean='protein/1SAP/1SAP.clean1.pdb'
protein_1SAP_name = '1SAP'
sequence_1SAP = 'MVKVKFKYKGEEKEVDTSKIKKVWRVGKMVSFTYDDNGKTGRGAVSEKDAPKELLDMLARAEREKK'
fragment_file_1SAP = 'protein/1SAP/aat000_03_05.200_v1_3'

protein_2EZK_clean='protein/2EZK/2EZK.clean1.pdb'
protein_2EZK_name = '2EZK'
sequence_2EZK = 'MIARPTLEAHDYDREALWSKWDNASDSQRRLAEKWLPAVQAADEMLNQGISTKTAFATVAGHYQVSASTLRDKYYQVQKFAKPDWAAALVDGRGASRRN'
fragment_file_2EZK = 'protein/2EZK/aat000_03_05.200_v1_3'

protein_2H5N_clean='protein/2H5N/2H5N.clean1.pdb'
protein_2H5N_name = '2H5N'
sequence_2H5N = 'MGLGRQSLNIMTFSGQELTAIIKMAKSMVMADGKIKPAEIAVMTREFMRFGILQDQVDLLLKASDSIEASQAVALIARMDEERKKYVASYLGVIMASDGDIDDNELALWTLISTLCGLPTMTVMEAINNMKNL'
fragment_file_2H5N = 'protein/2H5N/aat000_03_05.200_v1_3'


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
#added by tiger
protein_clean_file_name.append(protein_1AIL_clean)
protein_clean_file_name.append(protein_1AOY_clean)
protein_clean_file_name.append(protein_1CC5_clean)
protein_clean_file_name.append(protein_1FWP_clean)
protein_clean_file_name.append(protein_1SAP_clean)
protein_clean_file_name.append(protein_2EZK_clean)
protein_clean_file_name.append(protein_2H5N_clean)

protein_name = []
protein_name.append(protein_1ENH_name)
protein_name.append(protein_1GB1_name)
protein_name.append(protein_1GYZ_name)
protein_name.append(protein_1I6C_name)
protein_name.append(protein_1VII_name)
protein_name.append(protein_4ICB_name)
protein_name.append(protein_1BBO_name)
protein_name.append(protein_1MN3_name)
protein_name.append(protein_2IMU_name)
protein_name.append(protein_2JUJ_name)
protein_name.append(protein_3GWL_name)
protein_name.append(protein_1FD4_name)

#added by tiger
protein_name.append(protein_1AIL_name)
protein_name.append(protein_1AOY_name)
protein_name.append(protein_1CC5_name)
protein_name.append(protein_1FWP_name)
protein_name.append(protein_1SAP_name)
protein_name.append(protein_2EZK_name)
protein_name.append(protein_2H5N_name)


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

#added bu tiger
sequence_name.append(sequence_1AIL)
sequence_name.append(sequence_1AOY)
sequence_name.append(sequence_1CC5)
sequence_name.append(sequence_1FWP)
sequence_name.append(sequence_1SAP)
sequence_name.append(sequence_2EZK)
sequence_name.append(sequence_2H5N)


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

#added by tiger
fragment_name.append(fragment_file_1AIL)
fragment_name.append(fragment_file_1AOY)
fragment_name.append(fragment_file_1CC5)
fragment_name.append(fragment_file_1FWP)
fragment_name.append(fragment_file_1SAP)
fragment_name.append(fragment_file_2EZK)
fragment_name.append(fragment_file_2H5N)


file_directory_path = '/home/qcq/tree_based_memory/'
iterate = 1



want_to_run = ['1ENH','1GB1','1GYZ','1I6C','1VII','4ICB',
               '1BBO','1MN3','2IMU','2JUJ','1FD4','3GWL'\
               ,'1AIL','1AOY','1CC5','1FWP','1SAP','2EZK'\
               ,'2H5N'
               ]

KT_list_input = []
KT_list_input.append([32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66])    # this temperture list is the idealest
KT_list_input.append([32,34,36,38,40,41,42,43,44,45,52,54,56,58,60,62,64,66])    # this temperture list is the second 
KT_list_input.append([36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70])    # this temperture list is the third
KT_list_input.append([32,34,36,38,40,41,42,43,48,50,52,54,56,58,60,62,64,66])    # this temperture list is the last
KT_list_input.append([i*2+32 for i in range(1,19)])

FA_list = [1,4,8]  #times of 
#FA_list.append(1)

MC_list = [10,40,80]  # Mc times' list


