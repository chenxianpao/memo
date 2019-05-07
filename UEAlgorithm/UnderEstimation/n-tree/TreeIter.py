#!/usr/bin/python  
# -*- coding: utf-8 -*- 
'''
Created on 2014-4-24

@author: qcq
'''
from collections import deque
 #==============================================================================
 # TreeIter.py:  Iter the tree node for Python
 #==============================================================================
#--
# $Revision: 1.7 $
# $Date: 2000/03/29 22:36:12 $
# modified by qcq in 2014/3/27
# modified by qcq in 2014/4/23 增加了遍历树的非递归版本的广度优先和深度优先
#--

#================================================================
# Contents
#----------------------------------------------------------------
# class IterTree:  Generic n-ary tree node Iter object using recursion
# class IterTreeNonRecursive:  Generic n-ary tree node Iter object without recursion
#----------------------------------------------------------------

class IterTree:
    '''
    used for operation of tree iteration with recursive procedure 
    '''
    def __init__(self):
        
        self.nodes = []
        self.nodes_width = []
        self.contain_flag = False
        
        
        

    def iter_tree(self, root):
        '''Iterating the whole tree from the root(One kind of depth first search)
        '''
        #print root.value
        if len(root.childList) != 0:
            for i in range(len(root.childList)):
                #print root.nthChild(i).value,
                self.nodes.append(root.nthChild(i))
                self.iter_tree(root.nthChild(i))
        else:
            return
        
    def iter_tree_depth_first(self, root):
        '''
        depth-first search the whole tree
        '''
        if root is None:
            return 
        if root.nChildren() == 0:
            self.nodes_width.append(root)
            return 
        #count = 0
        self.nodes_width.append(root)
        self.iter_tree_depth_first(root.nthChild(0))
        if root.nChildren() == 1:
            #self.nodes_width.append(root.nthChild(0).value)
            return 
        else:
            for i in range(1,root.nChildren()):
                #count = count + 1
                #self.nodes_width.append(root.nthChild(i).value)
                self.iter_tree_depth_first(root.nthChild(i))
        return 
                
        
    def child_nodes(self):
        '''
        return the list which contains all tree node obiect
        '''
        return self.nodes
        
    def tree_size(self):
        return len(self.nodes)
        
    def contains(self,root, node):
        '''
        judge if the node in the tree rooted in root
        '''
        if node.value in self.nodes:
            return 1
        return 0

    
class IterTreeNonRecursive:
    '''
    used for operation of tree iteration without recursive procedure 
    '''
    def __init__(self):
        self.nodes_width = []
        self.nodes_depth = []
        self.depth = 0
    '''
    width first with queue
    '''    
    def iter_tree_width(self, root):
        queue = deque([])
        queue.append(root)
        
        while queue:
            #print 'yes'
            temp = queue.popleft()
            print temp
            self.nodes_width.append(temp)
            for i in temp.childList:
                queue.append(i)
                
    def iter_tree_depth(self, root):
        '''depth first with stack
        '''
        stack = []
        #self.depth = 1
        stack.append(root)
        while stack:
            temp = stack.pop()
            #print temp
            self.nodes_depth.append(temp)
            for i in temp.childList[::-1]:
                
                stack.append(i)
            #self.depth = self.depth + 1 
            
    def tree_size(self):
        return len(self.nodes_width)
            
        
        
#===============================================================================
# tree_size  return the size of the tree, This code referenced Tman's code style
#===============================================================================
def tree_size(root):
    if root is None:
        return 0
    if root.nChildren() == 0:
        return 0
    count = 0
    count = 1 + tree_size(root.nthChild(0))
    if root.nChildren() == 1:
        return count
    else:
        for i in range(1,root.nChildren()):
            count = count + 1
            count = count + tree_size(root.nthChild(i))
    return count

#===============================================================================
# tree_depth return the depth of the tree
#===============================================================================
def tree_depth(root):
    if root is None:
        return 0
    if root.nChildren() == 0:
        return 0
    depth = 0
    maxdepth = 0
    maxdepth = 1 + tree_depth(root.nthChild(0))
    if root.nChildren() == 1:
        return maxdepth
    else:
        for i in range(1, root.nChildren()):
            depth = 1 + tree_depth(root.nthChild(i))
            maxdepth = max(maxdepth, depth)
            
    return maxdepth
        
 #==============================================================================
 # tree_conains judge is one node exists in the tree rooted in root with Recursive
 #==============================================================================
def tree_contains(root, node):
    if root is None:
        #print 'False1'
        return False
    if root.nChildren() == 0:
        #print str(root.equals(node)) + '2'
        return root.equals(node)
    flag = False
    
    #self.nodes_width.append(root)
    if root.equals(node):
        #print 'True' + '3'
        return True
    flag = tree_contains(root.nthChild(0), node)
    if flag:
        return flag
    if root.nChildren() == 1:
        #self.nodes_width.append(root.nthChild(0).value)
        #print str(flag) + '4'
        return flag
    else:
        for i in range(1,root.nChildren()):
            #count = count + 1
            #self.nodes_width.append(root.nthChild(i).value)
            flag = tree_contains(root.nthChild(i), node)
            if flag :
                return flag
    #print str(flag) + '5'
    return flag

def tree_containsNonRecursive(root, node):
    '''depth first with stack
        '''
    stack = []
    #self.depth = 1
    stack.append(root)
    while stack:
        temp = stack.pop()
        #print temp
        #self.nodes_depth.append(temp)
        if temp.equals(node):
            return True
        for i in temp.childList[::-1]:
            
            stack.append(i)
        #self.depth = self.depth + 1 
    return False