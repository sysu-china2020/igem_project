# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 16:29:14 2020

@author: mi
"""

import os            #用于改变工作路径
import pandas as pd  #用于处理dataframe
import time  #用于计算程序运行用时

class ReadSeq:
    
    def __init__(self,path_csv,path_fasta1,path_fasta2,chromo,write_fasta,write_csv):
                                         #内部定义五个参数
        self.csv = path_csv         #包含底物信息的文件
        self.fasta1 = path_fasta1   #基因组文件名称前一半
        self.fasta2 = path_fasta2   #后一半。分开是因为染色体的编号要改
        self.chromo = chromo             #染色体编号
        self.write1 = write_fasta     #待写入的文件名
        self.write2 = write_csv
    
    def open_csv(self):                  #打开文件读取底物信息
        data = pd.read_csv(self.csv,sep = ',').iloc[:,[0,1,2,3,4,5,10]]
        substrate = data.loc[data['Chromosome']==self.chromo]
        return substrate
    
    def open_fasta(self):                #读取染色体fasta文件
        this_path = self.fasta1 + self.chromo + self.fasta2
        fin = open(this_path)
        line = fin.readlines()           #逐行读入
        sline = line[1:]                 #除去fasta格式文件第一行的解释信息 
        rownum = len(sline)
        for i in range(0,rownum):
            sline[i] = sline[i].strip();
        fin.close()
        return sline                     #返回该条染色体的所有序列，列表形式
    
    def write_seq(self):                 
        #根据csv中的起始终止位点信息，再fasta文件中定位序列，写入新的fasta文件
        substrate = self.open_csv()
        start = substrate.iloc[:,4] #获取起始位点,str,不包括表头
        end = substrate.iloc[:,5]   #获取终止位点
        subRow = substrate.shape[0]           #获取行数
        chromo_seq = self.open_fasta()
        col = len(chromo_seq[0]);           #列数/每行字符串长度 减掉末尾空格！
        substrate['Sequence'] = 0
        substrate['SequenceLength'] = 0
        seqDict = dict();                   #构建一个字典
        for j in range(0,subRow):           #准备提取序列
            seqStart = int(start.iloc[j]);         #将浮点数转换为整型再运算
            seqEnd = int(end.iloc[j]);
            startRow,startCol = divmod(seqStart,col) #记录初始行列
            endRow,endCol = divmod(seqEnd,col)       #记录末尾行列
            if startCol == 0:                 #起始位点编号刚好整除列数的情况
                startCol = col
                startRow = startRow - 1
            if endCol == 0:                          #结尾位点同理
                endRow = endRow - 1
            length = seqEnd-seqStart                 #记录序列长度-1
            mySeq = ''                               #用以拼接字符串
            while startRow <= endRow:
                mySeq = mySeq + chromo_seq[startRow]
                startRow += 1;
            mySeq = mySeq[(startCol-1):(startCol+length)] + ' ' 
                                                     #防止后面位置不够
            if substrate.iloc[j,3]=="-":                  #负链上碱基的转换
                mySeq = mySeq.replace("A","B").replace("T","A").replace("B","T")
                mySeq = mySeq.replace("C","B").replace("G","C").replace("B","G")                           
            seqDict[str(j)] = mySeq                        #k字典存储了所有序列
            substrate.iloc[j,7] = mySeq
            lengthSeq = len(mySeq)
            substrate.iloc[j,8] = lengthSeq
        
        
        substrate = substrate[substrate.iloc[:,8]>=15]
        subRow2 = substrate.shape[0]
        with open(self.write1 ,'a') as fout:
            for i in range(0,subRow2):
                info = '>'+substrate.iloc[i,0]+'_'+substrate.iloc[i,1]+\
                    '_chromosome'+str(substrate.iloc[i,2]) +\
                    '_'+substrate.iloc[i,3]+' '+ str(substrate.iloc[i,4])+\
                    '-' + str(substrate.iloc[i,5]) + " Affinity: " + \
                    str(substrate.iloc[i,6])
                fout.write(info+"\n"+seqDict[str(i)]+"\n\n")
            fout.close()
                
        if self.chromo == '1':
            substrate.to_csv(self.write2, mode = 'a', header=True, index=False)
        else:
            substrate.to_csv(self.write2, mode = 'a', header=False, index=False)
        

def main():    
    startTime = time.perf_counter()
    read_path = 'D:\\HFJ data\\4大二下\\python\\iGEM'
    chromosome = ["1","2","3","4","5","6","7","8","9","10","11","12","13",
                  "14","15","16","17","18","19","20","21","22","X","Y","MT"]    
    csv_ADAR1_Exogenous = 'ADAR1 with affinity.csv'
    csv_ADAR1_Endogenous = 'Endogenous ADAR1 (HEK293) with affinity.csv'
    fasta1 = 'Human_Chromosome_all/chromosome_'
    fasta2 = '_GRCh37.fasta'
    myFasta1 = 'output_all_ADAR1_Exogenous_DNA.fasta'
    myFasta2 = 'output_all_ADAR1_Endogenous_DNA.fasta'
    myCSV1 = 'output_all_ADAR1_Exogenous_DNA.csv'
    myCSV2 = 'output_all_ADAR1_Endogenous_DNA.csv'
    os.chdir(read_path)
    for k in range(25):
        chromo = chromosome[k] #每次循环，只考虑一种染色体，避免多次打开同个文件
        for_exo = ReadSeq(csv_ADAR1_Exogenous,fasta1,fasta2,chromo,myFasta1,myCSV1)
        for_exo.write_seq()    #处理外源底物
        for_endo = ReadSeq(csv_ADAR1_Endogenous,fasta1,fasta2,chromo,myFasta2,myCSV2)
        for_endo.write_seq()   #处理内源底物
    endTime = time.perf_counter()
    print(endTime - startTime) #运用Time模块，记录一下程序运行的时间
    
if __name__ == '__main__':
    main()


