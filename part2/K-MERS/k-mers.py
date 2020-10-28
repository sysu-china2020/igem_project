import os
import pandas as pd
import numpy as np
from collections import Counter
import matplotlib.pyplot as plt
import math


#（0，1）标准化
def standFunc(rawData,k):
    Min = rawData.iloc[:,6].min()
    Max = rawData.iloc[:,6].max()
    norFuc = lambda x : (x - Min) / (Max - Min);
    rawData['standardization'] = rawData.iloc[:,6].apply(norFuc)
    Affinity = pd.DataFrame(sorted(rawData.iloc[:,9]))
    Affinity.plot.density(figsize = (16,8),xlim = [-0.1,1.1])
    plt.savefig(str(k)+'-mers'+'standardization.png',dpi =1200,\
                bbbox_inches = 'tight')
    Affinity.describe()
    positive = rawData[rawData['sigmoid']>=0.002]
    negative = rawData[rawData['sigmoid']<=0.002]
    return positive,negative

#Sigmoid函数标准化,个人觉得标准化的效果更好,后续步骤基于标准化
def sigFunc(rawData,k):
    sigFuc = lambda x : 1.0 / (1 + np.exp(-float(x)))
    rawData['sigmoid'] = rawData.iloc[:,6].apply(sigFuc)
    Affinity = pd.DataFrame(sorted(rawData.iloc[:,9]))
    Affinity.plot.density(figsize = (16,8),xlim = [0.4,1.1])
    plt.savefig(str(k)+'-mers'+'SigmoidFunc.png',dpi =1200,\
                bbbox_inches = 'tight')
    plt.clf()
    Affinity.describe()
    positive = rawData[rawData['sigmoid']>=0.62]
    negative = rawData[rawData['sigmoid']<=0.58]
    return positive,negative


def cutIntoMers(dataSet,k):
    rows = dataSet.shape[0]
    allMers = []
    for i in range(rows):
        sequence = dataSet.iloc[i,7]
        sequence = sequence.replace(" ","")
        length = len(sequence)
        for j in range(length-k+1):
            mers = sequence[j:j+k]
            allMers.append(mers)
    result = Counter(allMers)
    return result

def richness(k):
    positiveSet,negativeSet = sigFunc(myData,k)
    positivedic = cutIntoMers(positiveSet,k)
    negativedic = cutIntoMers(negativeSet,k)
    positiveMers = pd.DataFrame(pd.Series(positivedic),columns  = ['Merspos'])
    negativeMers = pd.DataFrame(pd.Series(negativedic),columns  = ['Mersneg'])
    Mers = positiveMers.join(negativeMers)
    Mers['richness'] = list(map(lambda x,y:math.log((x + 1)/(y + 1)),\
        Mers.iloc[:,0],Mers.iloc[:,1]))
    return Mers

def graph(k):
    Mers = richness(k)
    Mers = Mers.dropna()
    richnessValue = Mers.iloc[:,2]
    mers = pd.Series(richnessValue.index)
    positive10 = richnessValue.sort_values().tail(10)
    negative10 = richnessValue.sort_values().head(10)
    aPositive10 = pd.Series(positive10.index)
    aNegative10 = pd.Series(negative10.index)
    richnessValue.describe()
    fig, ax = plt.subplots(figsize=(18,10),dpi =1000)
    plt.xticks([])
    ax.scatter(mers, richnessValue,s=2,color = 'blueviolet')
    ax.scatter(aPositive10, positive10,s=10,color = 'coral')
    for i in range(len(positive10)):
        plt.annotate(aPositive10[i],(aPositive10[i],positive10[i]),\
                     xytext = (2,2), textcoords = 'offset points',rotation=25)
    ax.scatter(aNegative10, negative10,s=10,color = 'turquoise')
    for j in range(len(negative10)):
        plt.annotate(aNegative10[j],(aNegative10[j],negative10[j]),\
                     xytext = (2,2),textcoords = 'offset points',rotation=25)
    ax.set_xlabel(str(k)+'-mers',fontsize = 20)
    ax.set_ylabel('richness',fontsize = 20)
    plt.savefig(str(k)+'-mers'+'Richness.png')
    plt.clf()


myPath = 'D:\\HFJ data\\Programming\\python\\iGEM'
os.chdir(myPath)
myData = pd.read_csv('output_all_ADAR1_Endogenous_DNA.csv',sep = ',')
graph(5)
graph(6)


