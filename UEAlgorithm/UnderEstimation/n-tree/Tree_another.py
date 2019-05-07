#!/usr/bin/python  
# -*- coding: utf-8 -*- 
'''
Created on 2014-3-27

@author: qcq
'''
class treeModel:

    '''tree view'''

    def __init__(self,Id,value,fatherId):

        self.Id=Id

        self.value=value

        self.fatherId=fatherId

    def show(self):

        return self.value

# 树的遍历和展示

class treeShow:

    '''tree show'''

    logList = [treeModel(0,'addTree',0)]  #记录已经遍历过的节点

    writtenList = [treeModel(0,'addTree',0)]  #记录已经打印出的节点

    def __init__(self,rootId,list):

        self.rootId = rootId

        self.list=list

    #通过Id获取节点 

    def getModelById(self,Id):

        for t in self.list:

            if t.Id == Id:

                return t

        return None

    #判断是否有子节点

    def haveChild(self,t):

        for t1 in self.list:

            if t1.fatherId == t.Id and not self.IsInLogList(t1):

                return True

        return False

    #获取第一个没有遍历的子节点

    def getFirstChild(self,t):

        for t1 in self.list:

            if t1.fatherId == t.Id and not self.IsInLogList(t1):

                return t1

        return None

    #判断某节点是否已经被遍历

    def IsInLogList(self,t):

        for t1 in self.logList:

            if t1.Id == t.Id:

                return True

        return False

    #判断某节点是否已经打印

    def IsInWrittenList(self,t):

        for t1 in self.writtenList:

            if t1.Id == t.Id:

                return True

        return False

    #获取父节点

    def getFatherTree(self,t):

        for t1 in self.list:

            if t1.Id == t.fatherId:

                return t1

        return None

    #遍历打印

    def show(self):

        currentTree = self.getModelById(self.rootId)

        s = '  '

        strNum = 1

        while(True):

            if self.haveChild(currentTree):

                if not self.IsInWrittenList(currentTree):

                    print s*strNum,currentTree.show()

                    self.writtenList.append(currentTree)

                currentTree = self.getFirstChild(currentTree)

                strNum += 1

                continue

            else:

                if(currentTree.Id == self.rootId):

                    break

                else:

                    if not self.IsInWrittenList(currentTree):

                        print s*strNum,currentTree.show()

                    self.logList.append(currentTree)

                    currentTree = self.getFatherTree(currentTree)

                    strNum -= 1

                    continue