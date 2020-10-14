import pandas as pd 
import numpy as np
import os
import glob
os.chdir('/home/fjy/16种癌症/PCCn')
path = '/home/fjy/16种癌症/PCCn'
os.listdir(path)
file = glob.glob(os.path.join(path, "*.npy"))

from math import sqrt
def multipl(a,b):
	sumofab=0.0
	for i in range(len(a)):
		temp=a[i]*b[i]
		sumofab+=temp
	return sumofab
 
def corrcoef(x,y):
	n=len(x)
	sum1=sum(x)
	sum2=sum(y)
	sumofxy=multipl(x,y)
	sumofx2 = sum([pow(i,2) for i in x])
	sumofy2 = sum([pow(j,2) for j in y])
	num=sumofxy-(float(sum1)*float(sum2)/n)
	den=sqrt((sumofx2-float(sum1**2)/n)*(sumofy2-float(sum2**2)/n))
	if den == 0:
		return 0
	else:
		return num/den
PCCn = []
for f in file:
    PCCn.append(np.load(f))

os.chdir('/home/fjy/16种癌症')
PCC = np.load('PCC.npy')


alldeltPCC = []
for i in range(len(PCC[0])):
	for j in range(len(PCCpairs[i][0])):
		alldeltPCC.append(PCCpairs[i][:,j] - PCC[:,i])

PCC = pd.DataFrame(PCC)
alldeltPCC = []
for i in range(len(PCC.iloc[0])):
	for j in range(len(PCCpairs16_list[i].iloc[0])):
		alldeltPCC.append(PCCpairs16_list[i].iloc[:,j] - PCC[i])


deltPCC = np.array(alldeltPCC).T
edges = pd.DataFrame(np.loadtxt("edges.csv",dtype=str))

deltPCC_all = pd.DataFrame(np.hstack((edges,deltPCC)))#151892*6549

os.chdir('/home/fjy')
reactom_pathways = pd.read_csv('reactom.pathways.csv',delimiter=',')

os.chdir('/home/fjy/16种癌症')
test = pd.merge(reactom_pathways,deltPCC_all,left_on='edges',right_on= 0 ,how = 'inner')
orderpathways = test.sort_values(by = 'V2',axis = 0,ascending = True)

temp = orderpathways.iloc[:,1].tolist()
from collections import Counter
import collections
def counter(temp):
    return Counter(temp)

count = counter(temp)
pathways_number = collections.Counter(count).most_common()
delet = []
for i in range(len(pathways_number)):
	number = pathways_number[i][1]
	if number <= 3:
		delet.append(i)

del pathways_number[216:224]

pathways_number = pd.read_csv('pathways_numbernew.csv',delimiter='\t', header = None)

temppathway = (np.unique(orderpathways.iloc[:,1].tolist())).tolist()
index = list()
filter_pathways = list()
for i in temppathway:
	for j in pathways_number.iloc[:,0]:
		if i == j:
			index.append(pathways_number.iloc[(pathways_number.iloc[:,0].tolist()).index(j),1])
			filter_pathways.append(pathways_number.iloc[(pathways_number.iloc[:,0].tolist()).index(j),0])

orderdeltPCC = []
for i in filter_pathways:
	orderdeltPCC.append(np.array(orderpathways)[np.array(orderpathways)[:,1] == i,:])

np.savetxt('index.csv',index,delimiter=',')


score = []
for i in range(len(orderdeltPCC)):
	score.append(sum(abs(orderdeltPCC[i][:,6:]))/index[i])

allscore = pd.DataFrame(score)
normal16score = pd.DataFrame(score)

np.savetxt('score16.csv',allscore,delimiter=',')
np.savetxt('normalscore16.csv',normal16score,delimiter=',')

