import numpy as np

a=np.array([[1,2,3,4],[5,6,7,8],[9,10,11,12],[13,14,15,16]])
b=a[1:2,]
#print b
c=np.array([[1,2,3,4]])
d=np.array([[5,6,7,8]])
e=c<d
#print c.shape[0]
#print c[0,]
usr=np.array([[1,2,3,4]])
#print usr[0,1]
#print usr
temp=np.array([[0,0,0,0]])
xcurrent_energy_score = [np.array([[0,0,0,0]]) for i in range(0, 10)]  
#print xcurrent_energy_score
f=np.zeros((4, 4))
#print f
unit_matrix = np.zeros((4,4))
for i in range(0, 4):
    unit_matrix[i,i] = 1
#print unit_matrix
a= []
a.append(1)
a.append(2)
a.append(3)
a.append(4)
a.insert(4,5)
#print a[-2:]
#print a
b=[]
b.append(a.index(1))
b.append(a.index(2))
b.append(a.index(3))
#print len(b)
'''
usr_bound = ([0,0,0,0],[100,100,100,100])
bound_sum=0.0
for i in usr_bound:
    bound_sum = bound_sum + (i[1] - i[0])  
    print i[1][1]
    print i[1] 
'''
usr_bound = ([0,100],[0,100],[0,100],[0,100])
bound_sum=0.0
for i in range(0,4):
    bound_sum=bound_sum+(usr_bound[i][1]-usr_bound[i][0])
for bound in usr_bound:
    print bound[0]
convertedusr_vector=np.array([[1,2,3,4]])   
print  convertedusr_vector[:,:-1].sum(axis=1)[0,0]
#convertedusr_vector[0, len(usr)] = 1 - convertedusr_vector[:,:-1].sum(axis = 1)[0, 0]