# DataMining
数据挖掘课程实验

## 实验一：频繁模式挖掘---DNA序列上的Motif发现
###[【要求】](https://github.com/suliangxd/DataMining/tree/master/doc)
1.导入Sequences.txt（/1.txt）文件中的DNA序列数据，共n=30条序列，每条长度m=281个字符

2.分别采用贪心算法、Gibbs抽样算法搜索长度w=30的频繁序列模式作为Motif，要求输出Motif在每个序列中的
  起始位置、Consensus Motif序列、得分score值以及运行时间；

3.比较贪心算法和Gibbs抽样算法的效果及效率。

###【Gibbs算法流程】
#####具体可参考[PDF](https://github.com/suliangxd/DataMining/tree/master/%E9%A2%91%E7%B9%81%E6%A8%A1%E5%BC%8F%E6%8C%96%E6%8E%98/doc)
1. 随机起点S[s1,...sn],(n行)  
2. 随机挑选一行,记为x   
3. 计算除了第X行以外的其余行的概率分布矩阵P[4][W]    
4. 计算所有剩余行（n-1行）的非匹配串的概率列向量Q[4]    
5. 遍历第x行，得到第x行所有起点关于P矩阵的概率P(a|P)    
6. 遍历第x行，得到第x行所有起点关于Q向量的概率P(a|Q)      
7. 5/6，得到关于x行所有起点的概率分布     
8. 挑选出最大概率点作为起点Sx     
9. 重复2-8，直到score没有增长     
10. 输出结果      
   
