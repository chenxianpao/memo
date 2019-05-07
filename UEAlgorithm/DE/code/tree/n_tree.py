#!/usr/bin/python  
# -*- coding: utf-8 -*- 
'''
Created on 2014-3-26

@author: qcq
'''
from Tree import *

tree = Tree(None,99)
tree_1 = Tree(tree,88)
tree_2 = Tree(tree_1,89)
tree_3 = Tree(tree_1,238)
tree_4 = Tree(tree_1,242)
tree_5 = Tree(tree_3,2342)
tree_root = Tree(None,78)

print tree_4.fullPath()
node = NodeId(tree_4.fullPath())
abc = IterTree()
abc.iter_tree(tree)
print abc.child_nodes()
print tree_size(tree)

print '*' * 18
print tree_depth(tree)
print tree_size(tree)

print tree_contains(tree,tree_root)