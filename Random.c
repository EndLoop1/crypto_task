#include"head.c"

int main()
{
    int x,i,gap;

    srand((unsigned int)time(NULL));

    for (i = 0; i < N; i++)  pos[i] = -1;
    for (i = 0; i < SAMPLE_SIZE; i++) 
    {
        //x = rand() % N;
        x = rand_normal(50.0, 10.0);
        //x= rand_uniform(N);
        seq[i] = x;
        histogram[x]++;
        if (pos[x] != -1) 
        {
            gap = i - pos[x];
            if (gap < MAX_GAP) gap_hist[gap]++;
        }
        pos[x] = i;
    }

    double chi_sq = Chi_Square();
    int dof = Freedom_Degrees(N);
    double p_val = Chi_Square_P_Value(chi_sq, dof);
    double ent = Entropy();
    double max_ent = log2(N);
    double r1 = Autocorrelation();

    printf("\n=== 随机数测评结果 ===\n");
    printf("卡方检验(Chi-Square) = %.2f\n", chi_sq);
    printf("自由度 = %d\n", dof);
    printf("p值 = %.6f\n", p_val);
    if (p_val < 0.05) 
    {
        printf("结论：拒绝均匀分布假设（分布显著偏离均匀）\n");
    } else 
    {
        printf("结论：接受均匀分布假设（数据符合均匀分布）\n");
    }
    printf("\n熵(Entropy) = %.4f / %.4f (最大熵)\n", ent, max_ent);
    if (ent > 0.9 * max_ent) 
    {
        printf("结论：熵值与最大熵较接近，分布较均匀\n");
    } else 
    {
        printf("结论：熵值较低，分布不均匀或存在规律\n");
    }
    printf("\n自相关系数R(1) = %.4f\n", r1);
    if (fabs(r1) < 0.1) 
    {
        printf("结论：序列无明显自相关，随机性较好\n");
    } else 
    {
        printf("结论：序列存在自相关，随机性可能不足\n");
    }
    printf("\n");

    // 导出 .csv
    FILE* f1 = fopen("histogram.csv", "w");
    FILE* f2 = fopen("gap_hist.csv", "w");
    FILE *f3 = fopen("data_normal.csv", "w");
    FILE *f4 = fopen("data_uniform.csv", "w");
    fprintf(f1, "index,value\n");
    fprintf(f2, "gap,count\n");
    for (int i = 0; i < N; i++) 
    {
        fprintf(f1, "%d,%d\n", i, histogram[i]);
    }
    for (int i = 1; i < MAX_GAP; i++) 
    {
        if (gap_hist[i] > 0)
            fprintf(f2, "%d,%d\n", i, gap_hist[i]);
    }
    for (int i = 0; i < SAMPLE_SIZE; i++) 
    {
        fprintf(f3, "%.5f\n", rand_normal(50.0, 10.0));
        fprintf(f4, "%d\n", rand_uniform(N));
    }
    fclose(f1);
    fclose(f2);
    fclose(f3);
    fclose(f4);

    int bits = 8;
    export_seq_as_ascii_bits(seq, SAMPLE_SIZE, bits, "random_seq_ascii_bits.txt");

    return 0;
}