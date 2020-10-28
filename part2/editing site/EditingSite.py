# -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 09:02:41 2020

@author: ASUS
"""


import os,csv,time
import pandas as pd

def r_basic(): #分别读入编辑位点和dsRNA数据集
    global data,string
    read_path = 'D:\\HFJ data\\Programming\\python\\iGEM\\editing site'
    os.chdir(read_path)
    site = 'human_all_merged_editing-sites_add-strand'
    data = pd.read_csv(site,header = None, sep = '\t')
    
    #myString = 'output_all_ADAR1_Endogenous_DNA.csv'
    myString = 'output_all_ADAR1_Exogenous_DNA.csv'
    string = pd.read_csv(myString,sep = ',')
    
def readsite(k): #读入k号染色体上的dsRNA序列和被编辑位点
    global site1,string1,stringRow,siteRow
    site1 = data[data.iloc[:,0]==('chr'+k)] #判断编辑位点中的k号染色体
    string1 = string[string.iloc[:,2]== k] #判断自提底物中的k号染色体
    siteRow = site1.shape[0] #.shape[0]获取行数
    stringRow = string1.shape[0]

def w_basic(): #结果写入Excel
    global f,writer
    #f = open('output_all_ADAR1_Endogenous_DNA_with_ESS.csv','w',newline = '')
    f = open('output_all_ADAR1_Exogenous_DNA_with_ESS.csv','w',newline = '')
    writer = csv.writer(f)
    writer.writerow(["Substrate_id","Arm","Chromosome","Strand","Start","End",
                     "Affinity","Seq","Seq_Len","Editing_site","ES_Seq"])    
    
def jw(): #找出ESS，写入Excel
    for i in range(stringRow):
        strandString = string1.iloc[i,3]
        start = int(string1.iloc[i,4])
        end = int(string1.iloc[i,5])
        sequence = string1.iloc[i,7]
        sequence = '***' + sequence + '***'
        for j in range(siteRow):
            editSite = int(site1.iloc[j,1])
            strandSite = site1.iloc[j,2]
            if (start <= editSite <= end) & (strandString == strandSite) :
                ess = sequence[editSite-start:editSite-start+7]
                writer.writerow([string1.iloc[i,0],string1.iloc[i,1],
                                 string1.iloc[i,2],string1.iloc[i,3],
                                 string1.iloc[i,4],string1.iloc[i,5],
                                 string1.iloc[i,6],string1.iloc[i,7],
                                 string1.iloc[i,8],editSite,ess])

def main():
    startTime = time.perf_counter()
    r_basic()
    w_basic()
    chromosome = ["1","2","3","4","5","6","7","8","9","10","11","12","13",
                  "14","15","16","17","18","19","20","21","22","X","Y","M"]
    for k in range(25):
        readsite(chromosome[k])
        jw()
    endTime = time.perf_counter()
    print(endTime - startTime)
    
main()
f.close()
