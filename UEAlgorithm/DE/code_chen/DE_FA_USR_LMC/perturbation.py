#!/usr/bin/env python
# -*- coding: utf-8 -*-
import math
from rosetta import *
def fragment_localMC(pose,scorefxn,rejects_times,mover_3mer,fragment_times=1):
    insert_short_frag = RepeatMover(mover_3mer,fragment_times)   #apply short_inserts times 3-mer fragment assmebly
    folding_mover = SequenceMover()  #apply a sequence perturbation to protin conformation
    folding_mover.add_mover(insert_short_frag)  #
    initial_score = scorefxn(pose) #initial score of the protein object
    intermediate_pose = Pose()  #intermediary pose object in case change the original object
    accepted_pose = Pose()
    accepted_pose.assign(pose)
    intermediate_pose.assign(pose)
    reject_counter = 0 #plus 1 when rejected
    while True:
        folding_mover.apply(intermediate_pose)
        current_score = scorefxn(intermediate_pose)#get current score of the protein
        if(current_score <= initial_score):
            #pass
            reject_counter = 0
            initial_score = current_score
            accepted_pose.assign(intermediate_pose)
        else:
            intermediate_pose.assign(accepted_pose)
            reject_counter = reject_counter + 1
            if(reject_counter >= rejects_times):
                pose.assign(intermediate_pose)
                break
                
def fragment_localMC_simple(pose,intermediate_pose,accepted_pose,scorefxn,rejects_times,folding_mover):
    initial_score = scorefxn(pose) #initial score of the protein object
    accepted_pose.assign(pose)
    intermediate_pose.assign(pose)
    reject_counter = 0 #plus 1 when rejected
    while True:
        folding_mover.apply(intermediate_pose)
        current_score = scorefxn(intermediate_pose)#get current score of the protein
        if(current_score <= initial_score):
            #pass
            reject_counter = 0
            initial_score = current_score
            accepted_pose.assign(intermediate_pose)
        else:
            intermediate_pose.assign(accepted_pose)
            reject_counter = reject_counter + 1
            if(reject_counter >= rejects_times):
                pose.assign(intermediate_pose)
                break
        
