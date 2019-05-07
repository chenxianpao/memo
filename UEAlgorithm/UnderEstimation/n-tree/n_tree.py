#!/usr/bin/python  
# -*- coding: utf-8 -*- 
'''
Created on 2014-3-26

@author: qcq
'''
from Tree import *
from TreeIter import *

tree = Tree(None,99)
tree_1 = Tree(tree,88)
tree_2 = Tree(tree_1,89)
tree_3 = Tree(tree_1,238)
tree_4 = Tree(tree_1,242)
tree_5 = Tree(tree_3,2342)
tree_6 = Tree(tree_3, 566)
tree_7 = Tree(tree_1, 3242)
tree_8 = Tree(tree, 983)
tree_root = Tree(None,566)

'''

print tree_4.fullPath()
node = NodeId(tree_4.fullPath())

'''
abc = IterTree()
bcd = IterTreeNonRecursive()
'''
print '8'*18
bcd.iter_tree_depth(tree)
#print bcd.nodes_depth
for i in bcd.nodes_depth:
    print i
print '9'*18
'''
'''
abc.iter_tree(tree)
abc.iter_tree_depth_first(tree)
for i in abc.child_nodes():
    print i
print '*'*11
for i in abc.nodes_width:
    print i
print tree_size(tree)
'''
print '*' * 18
print tree_depth(tree)
print tree_size(tree)

print tree_contains(tree,tree_root)
print
print tree_containsNonRecursive(tree, tree_root)
