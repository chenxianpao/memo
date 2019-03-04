from numpy import *
import operator
import matplotlib
import matplotlib.pyplot as plt

def createDataSet():
    group = array([[1.0,1.1], [1.0,1.0], [0,0], [0,0.1]])
    labels = ['A', 'A', 'B', 'B']
    return group, labels

def classify0(inX, dataSet, labels, k):
    print dataSet
    dataSetSize = dataSet.shape[0]
    # print dataSetSize
    # calc distance. first calc the diff, second get sq, third get sum, finally extract a root
    diffMat = tile(inX, (dataSetSize, 1)) - dataSet
    print diffMat
    sqDiffMat = diffMat ** 2
    print sqDiffMat
    sqDistances = sqDiffMat.sum(axis=1)
    print sqDistances
    distances = sqDistances ** 0.5
    print distances
    sortedDistIndicies = distances.argsort()
    print sortedDistIndicies
    classCount = {}
    # calc least distance point
    for i in range(k):
        voteIlabel = labels[sortedDistIndicies[i]]
        classCount[voteIlabel] = classCount.get(voteIlabel, 0) + 1
    print classCount
    sortedClassCount = sorted(classCount.iteritems(),
                              key=operator.itemgetter(1), reverse=True)
    print sortedClassCount
    return sortedClassCount[0][0]



group, labels = createDataSet()
# print classify0([0, 0], group, labels, 3)

def file2matrix(filename):
    fr = open(filename)
    arrayOLines = fr.readlines()
    numberOfLines = len(arrayOLines)
    returnMat = zeros((numberOfLines, 3))
    classLabelVector = []
    index = 0
    for line in arrayOLines:
        line = line.strip()
        listFromLine = line.split('\t')
        returnMat[index, :] = listFromLine[0:3]
        classLabelVector.append(listFromLine[-1])
        index += 1
    return returnMat, classLabelVector

datingDataMat, datingLabels = file2matrix("./datingTestSet.txt")
# print datingLabels

def autoNorm(dataSet):
    minVals = dataSet.min(0)
    maxVals = dataSet.max(0)
    ranges = maxVals - minVals
    normDataSet = zeros(shape(dataSet))
    m = dataSet.shape[0]
    normDataSet = dataSet - tile(minVals, (m, 1))
    normDataSet = normDataSet / tile(ranges, (m, 1))
    return normDataSet, ranges, minVals


# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.scatter(datingDataMat[:,1], datingDataMat[:,2])
#
# plt.show()
print autoNorm(datingDataMat)