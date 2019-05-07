# -*- coding: utf-8 -*-
from rosetta import *
init()
import random

def read_mer3(file_data, last_fragment = False):
    first_line_flag = True
    read_one_of_three = 0
    fragment_torsion_angle = []
    return_angle = []
    for i,j in zip(file_data, range(0, len(file_data))):#range(802 * index , 802 * (index + 1))):
            
        if first_line_flag:
            first_line_flag = False    
            continue
        if i == '\n':    
            continue

        read_one_of_three = read_one_of_three + 1
        if read_one_of_three <= 3:
            if not last_fragment and read_one_of_three == 1:
                #temp_angle = []
                #print j,i[18:18+9], i[27:27+9], i[36:36+9]
                fragment_torsion_angle.append([i[18:18+9], i[27:27+9], i[36:36+9]])
            if last_fragment :
                fragment_torsion_angle.append([i[18:18+9], i[27:27+9], i[36:36+9]])
                #fragment_torsion_angle[j].append(i[18:18+9])
                #fragment_torsion_angle[j].append(i[27:27+9])
                #fragment_torsion_angle[j].append(i[36:36+9])
                pass
            if read_one_of_three == 3:
                read_one_of_three = 0     
            continue
    if last_fragment:
        rand_number = random.randint(0,199)
         
        return_angle.append(fragment_torsion_angle[rand_number * 3])
        return_angle.append(fragment_torsion_angle[rand_number * 3 + 1])
        return_angle.append(fragment_torsion_angle[rand_number * 3 + 2])
    else:
        rand_number = random.randint(0,199)
        return_angle = fragment_torsion_angle[rand_number]
    #print return_angle
    return return_angle
       
    

def read_fragment(protein, fragment_file_path):
    torsion_angle = []
    with open(fragment_file_path, "r") as f:
        file_data = f.readlines()
        sequence_length = len(file_data)/802 + 2
    for index in range(0, sequence_length - 2):
        if index == sequence_length - 3:   #for the last fragment to deal with three-set angle
            temp_list = read_mer3(file_data[802 * index : 802 * (index + 1)],last_fragment = True)
            for index_i in range(0, 3):
                torsion_angle.append(temp_list[index_i])
        #for not the last fragment to deal with only one set angle
        else:
            torsion_angle.append(read_mer3(file_data[802 * index : 802 * (index + 1)]))
    #print len(torsion_angle)
    #print torsion_angle[0]
    
    for protein_length_index in range(1, sequence_length + 1):
        protein.set_phi(protein_length_index, float(torsion_angle[protein_length_index - 1][0].strip()))
        protein.set_psi(protein_length_index, float(torsion_angle[protein_length_index - 1][1].strip()))
        protein.set_omega(protein_length_index, float(torsion_angle[protein_length_index - 1][2].strip()))
    

        
    #print len(torsion_angle)
'''
            
protein_sequence = "WIARINAAVRAYGLNYSTFINGLKKAGIELDRKILADMAVRDPQAFEQVVNKVKEALQVQ"
scorefxn=create_score_function("score3")      #score3设置       
Parent = pose_from_sequence(protein_sequence, "fa_standard")
to_centroid = SwitchResidueTypeSetMover("centroid")   #centroid->fullatom convert
to_centroid.apply(Parent)
print scorefxn(Parent)
fragment_file_path = 'aat000_03_05.200_v1_3'
read_fragment(Parent, fragment_file_path)
print scorefxn(Parent)
'''
