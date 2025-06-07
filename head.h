#ifndef RAND_UTILS_H
#define RAND_UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define SAMPLE_SIZE 1000000
#define N 100
#define MAX_GAP 1000

int histogram[N];  // 频率统计
int pos[N];  // 记录上次出现的位置
int gap_hist[MAX_GAP];  // 间隔分布统计
int seq[SAMPLE_SIZE];  // 存储生成的随机数序列

int rand_uniform(int n); //均匀分布的随机数
double M_Pi(); //计算圆周率
double rand_normal(double mean, double stddev); //正态分布的随机数

double Chi_Square(); //卡方检验
double Entropy(); //熵计算
double Autocorrelation(); //自相关函数
int Freedom_Degrees(int n); //计算自由度
double Chi_Square_P_Value(double chi_square_stat, int df); // 卡方分布的p值计算

#endif 
