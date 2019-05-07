#!/usr/bin/python  
# -*- coding: utf-8 -*- 
'''
Created on 2014-3-26

@author: qcq
'''

 #==============================================================================
 # tree.py:  Generic tree node object for Python
 #==============================================================================
#--
# $Revision: 1.7 $
# $Date: 2000/03/29 22:36:12 $
#modified by qcq in 2014/3/27
#--

#================================================================
# Contents
#----------------------------------------------------------------
# class Tree:  Generic n-ary tree node object
# class NodeId:  Represents the path to a node in an n-ary tree
#----------------------------------------------------------------

import string

class Tree:
    """ Generic n-ary tree node object

        Children are additive; no provision for deleting them.
        The birth order of children is recorded: 0 for the first
        child added, 1 for the second, and so on.
        
        modified by qcq in 2014/3/27. add the function for deleting one node in the tree structure

        Exports:
            Tree(parent, value=None)    Constructor
            .parent         If this is the root node, None, otherwise
                            the parent's Tree object.
            .childList      List of children, zero or more Tree objects.
            .value          Value passed to constructor; can be any type.
            .birthOrder     If this is the root node, 0, otherwise the
                            index of this child in the parent's .childList
            .nChildren()    Returns the number of self's children.
            .nthChild(n)    Returns the nth child; raises IndexError if
                            n is not a valid child number.
            .fullPath():    Returns path to self as a list of child numbers.
            .nodeId():      Returns path to self as a NodeId.
    """


# - - -   T r e e . _ _ i n i t _ _   - - -

    def __init__ ( self, parent, value=None ):
        """ Constructor for a Tree object.

            [ if (parent is None or a Tree object) ->
                if (parent is None) ->
                  return a new Tree object with no parent, no children,
                  and value (value)
                else ->
                  return a new Tree object added as the next new child
                  of parent (parent) and no children, and value (value)
            ]
        """
        #-- 1 --
        self.parent     =  parent
        self.value      =  value
        self.childList  =  []

        #-- 2 --
        #-[ if (parent is None) ->
        #     self.birthOrder  :=  0
        #   else ->
        #     parent           :=  parent with self added as its next child
        #     self.birthOrder  :=  size of parent's .childList
        #-]
        if  parent is None:
            self.birthOrder  =  0
        else:
            self.birthOrder  =  len(parent.childList)
            parent.childList.append ( self )


# - - -   T r e e . n C h i l d r e n   - - -

    def nChildren ( self ):
        """ [ return the number of children of self
            ]
        """
        return len(self.childList)


# - - -   T r e e . n t h C h i l d   - - -

    def nthChild ( self, n ):
        """ [ if (n is an integer) ->
                if (0 <= n < (number of self's children) ->
                  return self's (n)th child, counting from 0
                else ->
                  raise IndexError
            ]
        """
        return self.childList[n]


# - - -   T r e e . f u l l P a t h   - - -

    def fullPath ( self ):
        """Returns a list of child numbers from root to self.

          [ return a sequence [c0, c1, ..., ci] such that self is
            root.nthChild(c0).nthChild(c1). ... .nthChild(ci), or
            an empty list if self is the root ]
        """

        #-- 1 --
        result  =  []
        parent  =  self.parent
        kid     =  self

        #-- 2 --
        # [ result  +:=  child numbers from parent to root, in
        #                reverse order ]
        while  parent:
            result.insert ( 0, kid.birthOrder )
            parent, kid  =  parent.parent, parent

        #-- 3 --
        return result


# - - -   T r e e . n o d e I d   - - -

    def nodeId ( self ):
        """Returns the path to self in the tree as a NodeId.
        """
        #-- 1 --
        # [ fullPath  :=  sequence [c0, c1, ..., ci] such that self is
        #   root.nthChild(c0).nthChild(c1). ... .nthChild(ci), or
        #   an empty list if self is the root ]
        fullPath  =  self.fullPath()

        #-- 2 --
        return NodeId ( fullPath )
    
    def equals(self, node):
        '''judge if the two tree object is equal
        '''
        return self.value == node.value
    
    #===========================================================================
    # delete the node from the tree
    #===========================================================================
    def delete(self):
        if self.parent is None:
            return
        else:
            #temp = self.birthOrder
            self.parent.childList.remove(self.parent.childList[self.birthOrder])
            for i,j in zip(range(self.birthOrder + 1, self.parent.nChildren()), self.parent.childList[self.birthOrder + 1:]):
                j.birthOrder = j.birthOrder - 1
                
            
    



# - - - - -   c l a s s   N o d e I d   - - - - -

class NodeId:
    """Represents the location of a node in a tree as a path from the root.

      Exports:
        NodeId(path):
          [ if path is a list of zero or more nonnegative integers ->
              return a new NodeId object representing that node-id ]
        .path:      [ as passed to constructor ]
        .__str__():     [ return self as a string ]
        .find(root):
          [ if root is a Tree object ->
              if self describes a path to a node that exists in the
              tree rooted at (root) ->
                return the .value of that node
              else ->
                return None ]
        .isOnPath(node):
          [ if node is a Tree object ->
              if the path from the root to node is a prefix of self ->
                return 1
              else ->
                return 0 ]
    """

# - - -   N o d e I d . _ _ i n i t _ _   - - -

    def __init__ ( self, path ):
        """Constructor for the NodeId object
        """
        self.path  =  path


# - - -   N o d e I d . _ _ s t r _ _   - - -

    def __str__ ( self ):
        """Return self in displayable form
        """

        #-- 1 --
        # [ L  :=  a list of the elements of self.path converted to strings ]
        L  =  map ( str, self.path )

        #-- 2 --
        # [ return the elements of L concatenated and separated by "/" ]
        return string.join ( L, "/" )


# - - -   N o d e I d . f i n d   - - -

    def find ( self, node ):
        """Locate the tree node described by self and return its value
        """
        return self.__reFind ( node, 0 )


# - - -   N o d e I d . _ _ r e F i n d   - - -

    def __reFind ( self, node, i ):
        """Recursive node finding routine.  Starts at self.path[i:].

          [ if (node is a Tree object)
            and (0 <= i <= len(self.path)) ->
              if  i == len(self.path) ->
                return node's value
              else if self.path[i:] describes a path from node
              to some tree object T ->
                return T
              else ->
                return None ]   
        """

        #-- 1 --
        if  i >= len(self.path):
            return node.value       # We're there!
        else:
            childNo  =  self.path[i]

        #-- 2 --
        # [ if node has a child of number childNo ->
        #     child  :=  that child node
        #   else ->
        #     return None ]
        try:
            child  =  node.nthChild ( childNo )
        except IndexError:
            return None

        #-- 3 --
        # [ if (i+1) == len(self.path) ->
        #     return child
        #   else if self.path[i+1:] describes a path from node to
        #   some tree object T ->
        #     return T
        #   else ->
        #     return None ]
        return self.__reFind ( child, i+1 )


# - - -   N o d e I d . i s O n P a t h   - - -

    def isOnPath ( self, node ):
        """Is self's path to or through the given node?
        """

        #-- 1 --
        # [ nodePath  :=  path list for node ]
        nodePath  =  node.fullPath()

        #-- 2 --
        # [ if nodePath is a prefix of, or identical to self.path ->
        #     return 1
        #   else ->
        #     return 0 ]
        if  len(nodePath) > len(self.path):
            return 0        # Node is deeper than self.path

        for  i in range(len(nodePath)):
            if  nodePath[i] != self.path[i]:
                return 0    # Node is a different route than self.path

        return 1
    
class IterTree:
    '''
    used for operation of tree iteration
    '''
    def __init__(self):
        
        self.nodes = []
        
        
        

    def iter_tree(self, root):
        '''Iterating whole tree from the root
        '''
        #print root.value
        if len(root.childList) != 0:
            for i in range(len(root.childList)):
                #print root.nthChild(i).value,
                self.nodes.append(root.nthChild(i).value)
                self.iter_tree(root.nthChild(i))
        else:
            return
        
    def child_nodes(self):
        return self.nodes
        
    def tree_size(self):
        return len(self.nodes)
        
    def contains(root, node):
        '''
        judge if the node in the tree rooted in root
        '''
        if node.value in self.nodes:
            return 1
        return 0
    

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
 # tree_conains judge is one node exists in the tree rooted in root
 #==============================================================================
def tree_contains(root, node):
    if root is None:
        return False
    if root.equals(node):
        return True
    if root.nChildren() == 0:
        return False
    flag = False
    #print root.nthChild(0).value,
    if root.nthChild(0).equals(node):
        flag = True
        return flag
    flag = tree_contains(root.nthChild(0), node)
    if root.nChildren() == 1:
        return flag
    else:
        for i in range(1,root.nChildren()):
            #print root.nthChild(i).value,
            if root.nthChild(i).equals(node):
                flag = True
                return flag
            flag = tree_contains(root.nthChild(i), node)
    return flag
    
            
    

        
    