//
//  main.cpp
//  DNAMotif
//
//  Created by 苏亮 on 15/7/6.
//  Copyright (c) 2015年 苏亮. All rights reserved.
//  title:频繁模式挖掘---DNA序列上的Motif发现
/*
 1、导入Sequences.txt(1.txt)文件中的DNA序列数据，共n=30条序列，每条长度m=280个字符；
 2、分别采用贪心算法、Gibbs抽样算法搜索长度w=30的频繁序列模式作为Motif，
 要求输出Motif在每个序列中的起始位置、Consensus Motif序列、得分score值以及运行时间；
 3、比较贪心算法和Gibbs抽样算法的效果及效率。
 */
//  算法：贪心 and Gibbs

#include <iostream>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>

using namespace std;

const int N = 30;
const int M = 280;
const int W = 30;

class Profiles{
public:
    /*
     *输入原始矩阵
     */
    void init_input();
    
    /*
     *随机起始点数组或者随机行，只能一种情况
     */
    void randStartPos(int randRow = false);
    
    /*
     *获取统计矩阵, 并得到ProfileMax_array[W+1]
     */
    void getProfile(bool IsUseGibbs=false);
    
    /*
     *获取Consensus模板串
     */
    void getConsensus();
    
    /*
     *获取概率矩阵P
     */
    void getProbablity();
    
    /*
     *获取Q向量,Gibbs算法
     */
    void getProb_Qarray();
    
    /*
     *随机行所有起点的P概率分布
     */
    void getRandRow_ProbP();
    
    /*
     *随机行所有起点的Q概率分布
     */
    void getRandRow_ProbQ();
    
    /*
     *获取随机行P/Q的概率分布,返回最大概率的起点
     */
    int getProb_PdivQ();
    
    /*
     *获取当前分数
     */
    int  getScore();
    
    /*
     *贪心算法处理
     */
    void SolveByGreedy();
    
    /*
     *Gibbs算法处理
     */
    void SolveByGibbs();
    
    /*
     *输出结果
     */
    void output(int IsUseGibbs = false);
    
private:
    /*
     *原始矩阵 30*281
     */
    char Alignment_matrix[N+1][M+1];
    
    /*
     * 概率矩阵 4*30(W)
     */
    double Probability_matrix[4][W+1];
    
    /*
     *概率向量Q[4]，Gibbs算法使用
     */
    double Prob_Qarray[4];
    
    /*
     *随机行的关于P概率分布
     */
    double Prob_rowOfP[M+1];
    
    /*
     *随机行的关于Q概率分布
     */
    double Prob_rowOfQ[M+1];
    
    /*
     *随机行的P/Q概率分布
     */
    double Prob_PdivQ[M+1];
    
    /*
     *模式串
     */
    char  Consensus_array[W+1];
    
    /*
     *随机行，Gibbs算法使用
     */
    int   Random_row;
    
    /*
     *每行起始点数组 [30]
     */
    int   StartPos_array[N+1];
    
    /*
     *统计矩阵 4*30(W) A T G C
     */
    int   Profile_matrix[4][W+1];
    
    /*
     *每列最大的基因出现次数
     */
    int   ProfileMax_array[W+1];
    
    /*
     *最优分数
     */
    int   bestScore;
};

void Profiles::init_input()
{
    freopen("1.txt","r",stdin);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            cin >> this->Alignment_matrix[i][j];
        }
    }
}
void Profiles::randStartPos(int randRow)
{
    //随机任意行
    if (randRow == true)
    {
        this->Random_row = rand()%N;
        return ;
    }
    //随机起始位置
    for (int i = 0; i < N; i++)
    {
        int tmp = rand()%N;
        while(tmp+W>N)
            tmp = rand()%N;
        this->StartPos_array[i] = rand()%N;
    }
}
void Profiles::getProfile(bool IsUseGibbs)
{
    //使用Gibbs算法
    if (IsUseGibbs == true)
    {
        int deleteRow = this->Random_row;//被删除的行
        //列
        for (int i = 0; i < W; i++)
        {
            int A,T,G,C;
            A = T = G = C = 0;
            //行
            for (int j = 0; j < N; j++)
            {
                if (j == deleteRow)
                    continue;
                switch (this->Alignment_matrix[j][i+this->StartPos_array[j]])
                {
                    case 'A':
                        A++;
                        break;
                    case 'T':
                        T++;
                        break;
                    case 'G':
                        G++;
                        break;
                    case 'C':
                        C++;
                        break;
                    default:
                        break;
                }
            }
            this->Profile_matrix[0][i] = A;
            this->Profile_matrix[1][i] = T;
            this->Profile_matrix[2][i] = G;
            this->Profile_matrix[3][i] = C;
            int tmpMax;
            tmpMax = max(max(max(T,G),C),A);
            this->ProfileMax_array[i] = tmpMax;
        }
    }
    //不使用Gibbs算法
    else
    {
        //列
        for (int i = 0; i < W; i++)
        {
            int A,T,G,C;
            A = T = G = C = 0;
            //列
            for (int j = 0; j < N; j++)
            {
                switch (this->Alignment_matrix[j][i+this->StartPos_array[j]])
                {
                    case 'A':
                        A++;
                        break;
                    case 'T':
                        T++;
                        break;
                    case 'G':
                        G++;
                        break;
                    case 'C':
                        C++;
                        break;
                    default:
                        break;
                }
            }
            this->Profile_matrix[0][i] = A;
            this->Profile_matrix[1][i] = T;
            this->Profile_matrix[2][i] = G;
            this->Profile_matrix[3][i] = C;
            int tmpMax;
            tmpMax = max(max(max(T,G),C),A);
            this->ProfileMax_array[i] = tmpMax;
        }
    }
}

void Profiles::getProb_Qarray()
{
    int A,T,G,C;
    A = T = G = C = 0;
    for (int i = 0; i < N; i++)
    {
        int deleteRow = this->Random_row;
        if (i == deleteRow)
            continue;
        for (int j = 0; j < M; j++)
        {
            if (j == this->StartPos_array[i])
            {
                j += (W - 1);
                continue;
            }
            else
            {
                switch (this->Alignment_matrix[i][j])
                {
                    case 'A':
                        A++;
                        break;
                    case 'T':
                        T++;
                        break;
                    case 'G':
                        G++;
                        break;
                    case 'C':
                        C++;
                        break;
                    default:
                        break;
                }
            }
        }
    }
    int sum = A + T + G + C;
    this->Prob_Qarray[0] = A*1.0/sum;
    this->Prob_Qarray[1] = T*1.0/sum;
    this->Prob_Qarray[2] = G*1.0/sum;
    this->Prob_Qarray[3] = C*1.0/sum;
}

void Profiles::getConsensus()
{
    for (int j = 0; j < W; j++)
    {
        int maxrow = 0;
        int maxrowNum = this->Profile_matrix[0][j];
        //行
        for (int i = 1; i < 4; i++)
        {
            if (this->Profile_matrix[i][j] > maxrowNum)
            {
                maxrowNum = this->Profile_matrix[i][j];
                maxrow = i;
            }
        }
        switch (maxrow)
        {
            case 0:
                this->Consensus_array[j] = 'A';
                break;
            case 1:
                this->Consensus_array[j] = 'T';
                break;
            case 2:
                this->Consensus_array[j] = 'G';
                break;
            case 3:
                this->Consensus_array[j] = 'C';
                break;
            default:
                break;
        }
    }
}
void Profiles::getProbablity()
{
    for (int j = 0; j < W; j++)
    {
        int tmpsum = 0;
        for (int i = 0; i < 4; i++)//row
        {
            tmpsum+=this->Profile_matrix[i][j];
        }
        for (int i = 0; i < 4; i++)
        {
            this->Probability_matrix[i][j] = this->Profile_matrix[i][j]*1.0/tmpsum;
        }
    }
}

void Profiles::getRandRow_ProbP()
{
    int deleteRow = this->Random_row;
    //枚举随机行的起点
    for (int i = 0; i+W-1 < M; i++)
    {
        double tmpProb = 1.0;
        for(int j = 0; j < W; j++)
        {
            switch (this->Alignment_matrix[deleteRow][i+j])
            {
                case 'A':
                    tmpProb *= this->Profile_matrix[0][j];
                    break;
                case 'T':
                    tmpProb *= this->Profile_matrix[1][j];
                    break;
                case 'G':
                    tmpProb *= this->Profile_matrix[2][j];
                    break;
                case 'C':
                    tmpProb *= this->Profile_matrix[3][j];
                    break;
                default:
                    break;
            }
        }
        this->Prob_rowOfP[i] = tmpProb;
    }
}

void Profiles::getRandRow_ProbQ()
{
    int deleteRow = this->Random_row;
    //枚举随机行的列
    for (int i = 0; i+W-1 < M; i++)
    {
        double tmpProb = 1.0;
        switch (this->Alignment_matrix[deleteRow][i])
        {
            case 'A':
                tmpProb *= this->Prob_Qarray[0];
                break;
            case 'T':
                tmpProb *= this->Prob_Qarray[1];
                break;
            case 'G':
                tmpProb *= this->Prob_Qarray[2];
                break;
            case 'C':
                tmpProb *= this->Prob_Qarray[3];
                break;
            default:
                break;
        }
        this->Prob_rowOfQ[i] = tmpProb;
    }
}

int Profiles::getProb_PdivQ()
{
    int tmpRow = 0;
    double tmpMaxRowNum = 0.0;
    for (int i = 0; i+W-1< M; i++)
    {
        this->Prob_PdivQ[i] = this->Prob_rowOfP[i] / this->Prob_rowOfQ[i];
        if (this->Prob_PdivQ[i] > tmpMaxRowNum)
        {
            tmpRow = i;
            tmpMaxRowNum = this->Prob_PdivQ[i];
        }
    }
    return tmpRow;
}
int Profiles::getScore(){
    int score = 0;
    for (int i = 0; i < W; i++)
        score += this->ProfileMax_array[i];
    return score;
}
void Profiles::output(int IsUseGibbs){
    if(IsUseGibbs == true)
    {
        cout << "----------  Gibbs  ------------"<<endl<<endl;
    }
    else
    {
        cout << "----------  Greedy ------------"<<endl<<endl;
    }
    cout << "起始位置信息: " <<endl;
    for (int i = 0; i < N; i++)
        cout <<"第"<<i<<"行"<<": "<<this->StartPos_array[i]<<endl;
    cout << "Consensus Motif序列: " <<endl;
    for (int i = 0; i < N; i++)
        printf("%c",this->Consensus_array[i]);
    printf("\n");
    cout << "得分: ";
    cout << this->bestScore <<endl;
}
void Profiles::SolveByGreedy(){
    //cout <<"进入SolveByGreedy()函数"<<endl;
    this->init_input();   //输入原始矩阵
    //cout << "输入原始矩阵完成" <<endl;
    this->randStartPos(); //随机起点
    //cout << "随机起点完成" <<endl;
    this->getProfile();   //统计矩阵
    //cout << "统计矩阵完成" <<endl;
    this->getConsensus(); //得模板串
    //cout << "模板串完成" <<endl;
    this->getProbablity(); //统计概率
    //cout << "统计概率完成" <<endl;
    
    this->bestScore = 0;
    int tmpScore = this->getScore();
    while(tmpScore > this->bestScore)
    {
        this->bestScore = tmpScore; //更新
        for (int i = 0; i < N; i++)
        {
            int tmpcol = 0;
            double tmpcolProba = 0;//最优概率
            
            //枚举起始点，找最优
            for (int j = 0; j+W-1 < M; j++)
            {
                double nowcolProba = 1.0; //当前概率值
                for (int k = 0; k < W; k++){
                    switch (this->Alignment_matrix[i][j+k])
                    {
                        case 'A':
                            nowcolProba *= this->Probability_matrix[0][k];
                            break;
                        case 'T':
                            nowcolProba *= this->Probability_matrix[1][k];
                            break;
                        case 'G':
                            nowcolProba *= this->Probability_matrix[2][k];
                            break;
                        case 'C':
                            nowcolProba *= this->Probability_matrix[3][k];
                            break;
                        default:
                            break;
                    }
                }
                if (tmpcolProba < nowcolProba)
                {
                    tmpcolProba = nowcolProba;
                    tmpcol = j;
                }
            }
            //更新起始点位置,选取最优概率情况下的起始点
            this->StartPos_array[i] = tmpcol;
        }
        //更新起始点位置以后，需要重新计算Profile,Probablity
        this->getProfile();
        this->getConsensus();
        this->getProbablity();
    }
    this->output();
}

void Profiles::SolveByGibbs(){
    this->init_input();     //输入原始矩阵
    this->randStartPos();   //随机起点
    this->randStartPos(true); //随机行
    this->getProfile(true); //统计矩阵
    this->getConsensus();   //模板串
    this->getProbablity();  //概率矩阵
    this->getProb_Qarray(); //得到向量Q
    this->getRandRow_ProbP(); //随机行关于P的概率分布
    this->getRandRow_ProbQ(); //随机行关于Q的概率分布
    //随机行关于P/Q的概率分布,返回最大概率情况的起始点
    int MaxstartPos = this->getProb_PdivQ();
    int tmpScore = this->getScore();
    this->bestScore = 0;
    while(tmpScore > this->bestScore)
    {
        //更新
        this->bestScore = tmpScore;
        this->StartPos_array[this->Random_row] = MaxstartPos;
        
        //重复算法过程
        this->randStartPos(true); //随机行
        this->getProfile(true); //统计矩阵
        this->getConsensus();   //模板串
        this->getProbablity();  //概率矩阵
        this->getProb_Qarray(); //得到向量Q
        this->getRandRow_ProbP(); //随机行关于P的概率分布
        this->getRandRow_ProbQ(); //随机行关于Q的概率分布
        //随机行关于P/Q的概率分布,返回最大概率情况的起始点
        MaxstartPos = this->getProb_PdivQ();
        
        tmpScore = this->getScore();
    }
    this->output(true);
}
int main()
{
    srand(time(NULL));
    Profiles ProfilesTest_Greedy;
    Profiles ProfilesTest_Gibbs;
    ProfilesTest_Greedy.SolveByGreedy();
    ProfilesTest_Gibbs.SolveByGibbs();
    return 0;
}
