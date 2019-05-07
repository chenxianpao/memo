#!/usr/local/bin/python2.7
# encoding: utf-8
'''
This program used to convert the unix format txt to Dos txt file and vice versa.
This program used the linux shell command, so should be used on the linux platform



Author:QCQ
'''

import sys
import os
#from protein_information import *

def unix2dos(file_directory):
    '''
    parameter:
    file_directory: the unix file which want to convert
    '''
    with open(file_directory,'r') as f:
        os.popen('unix2dos ' + f.name)
        
def dos2unix(file_directory):
    '''
    parameter:
    file_directory: the Dos file which want to convert
    '''
    with open(file_directory,'r') as f:
        os.popen('dos2unix ' + f.name)
    
#unix2dos('test.txt')




