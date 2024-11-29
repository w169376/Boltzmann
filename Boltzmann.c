/*****************************************************************************80
    Copyright 2016 Songsong Wang <sswang@hit.edu.cn>
/*----------------------------------------------------------------------------80
    Author: Songsong Wang <sswang@hit.edu.cn>
    Date: 2016.01.04
------------------------------------------------------------------------------*/
#include<stdio.h>
#include <stdlib.h>
#include<math.h>
#include<time.h>
#include <string.h>  // 添加 string.h 头文件
#define SIZE 20
#define INPUT_SIZE 100                                          // 使用更高精度的数据类型（如果编译器支持）
#define ARRAY_SIZE  20
#define BOLTZMANN_CONSTANT 1.380650e-23                         // 波尔兹曼常量

long double Fact(long double n);                                         // 本程序多次用到阶乘，故定义用递归法求阶乘
void Boltzmann_System(int a[], int w[], int l);                 // 玻尔兹曼系统微观状态数求解
void Bose_System(int a[], int w[], int l);                      // 波色系统求解微观状态数
void Fermi_System(int a[], int w[], int l);                     // 费米系统求解微观状态数                                             // 基于波尔兹曼统计的相关求解
void read_int_array(int *array, int size, const char *prompt);
void calculate_common(int a[], int w[], int l, long double (*fact_func)(int), void (*print_func)(double));
void Boltzmann(double k0);


// 计算阶乘
long double Fact(long double n) {
    if (n < 0) return -1;
    if (n == 0) return 1; // 阶乘为0时返回1
    long double result = 1;
    for (int i = 2; i <= n; ++i) {
        result *= i;
    }
    printf("Fact(%d) = %Lf\n", (long double)result); //调试输出  
    return result;
}

void calculate_common(int a[], int w[], int l, long double (*fact_func), void (*print_func)(double)) {
    int n = 0, i = 0;
    int a1 = 0, w1 = 0;
    double N = 1, W = 1.0, A = 1, SUM = 0;

    for (i = 0; i < l; i++) {
        n += a[i];
        a1 = a[i];
        w1 = w[i];
        W *= pow(w1, a1);
        A *= fact_func(a1);

    }

    N = fact_func(n);
    SUM = N * W / A;
    print_func(SUM);
}
// 打印波尔兹曼系统总的微观状态数
void print_boltzmann_system(double sum) {
    printf("波尔兹曼系统总的微观状态数为：%g\n", sum);
}

// 打印波色系统总的微观状态数
void print_bose_system(double sum) {
    printf("波色系统总微观状态数%g\n", sum);
}

// 打印费米系统总的微观状态数
void print_fermi_system(double sum) {
    printf("费米系统总微观状态数%g\n", sum);
}

// 计算波尔兹曼系统微观状态数
void Boltzmann_System(int a[], int w[], int l) {
    calculate_common(a, w, l, Fact, print_boltzmann_system);
}

// 计算波色系统微观状态数
void Bose_System(int a[], int w[], int l) {
    calculate_common(a, w, l, Fact, print_bose_system);
}

// 计算费米系统微观状态数
void Fermi_System(int a[], int w[], int l) {
    calculate_common(a, w, l, Fact, print_fermi_system);
}
// 读取整数数组
void read_int_array(int *array, int size, const char *prompt) {
    int i;
    char input[INPUT_SIZE];
    printf("%s\n", prompt);
    for (i = 0; i < size; i++) {
        while (1) {
            if (fgets(input, sizeof(input), stdin) == NULL) {
                fprintf(stderr, "读取输入失败。\n");
                exit(1);
            }
            if (sscanf(input, "%d", &array[i]) != 1) {
                fprintf(stderr, "输入错误，请输入一个整数。\n");
                continue;
            }
            if (array[i] < 0 || array[i] >= SIZE) {
                fprintf(stderr, "输入错误，请输入一个正整数且不大于%d。\n", SIZE);
                continue;
            }
            break;
        }
    }
}

// 主函数
void main() {
    //const int SIZE = 20;
    const double K0 = 1.380650e-23;

    int a[SIZE];
    int w[SIZE];
    int k, i = 0;
    char input[INPUT_SIZE];

    printf("请输入能级数l\n");
while (1) {
    if (fgets(input, sizeof(input), stdin) == NULL) {
        fprintf(stderr, "读取输入失败。\n");
        return;
    }
    if (sscanf(input, "%d", &k) != 1 || k <= 0 || k > SIZE) {
        fprintf(stderr, "输入错误，请输入一个正整数且不大于%d。\n", SIZE);
        continue;
    }
    break;
}
    read_int_array(a, k, "请依次输入l个各能级上的粒子数");
    read_int_array(w, k, "请依次输入l个各能级上的简并度");

    Boltzmann_System(a, w, k);
    Bose_System(a, w, k);
    Fermi_System(a, w, k);
    Boltzmann(K0);
}


// 波尔兹曼统计的相关求解   
// 功能：基于波尔兹曼统计计算粒子分配函数Z、系统总粒子数N、总内能U、自由能F、定体热容Ct\Cv、Cr等
void Boltzmann(double k0) {
    const int size = 20;
    int wq[size];
    double aq[size];
    long double eq[size];
    int i, r;
    float T, T0 = 20.3, t1, TV;
    double p, z1, Z = 0, U = 0, Ct, Cr, Cv, UV0, C0, Ut, Uv, Ur;
    long double er;

    // 初始化变长数组
    memset(wq, 0, sizeof(wq));
    memset(aq, 0, sizeof(aq));
    memset(eq, 0, sizeof(eq));

    srand((unsigned)time(NULL));
    printf("请输入不小于224的温度（单位为K）T：\n");
    while (1) {
        if (scanf("%f", &T) != 1 || T < 224) {
            fprintf(stderr, "输入错误，请输入一个不小于224的温度。\n");
            continue;
        }
        break;
    }

    p = k0 * T;
    p = 1 / p;
    t1 = T / T0;
    er = 140 * powl(t1, 1.5);

    for (i = 0; i < size; i++) {
        r = (i + 1) * (i + 1);
        eq[i] = -13.6 * 1.60217e-19L / r;
        wq[i] = rand() % 10;
        double q = -p * eq[i];
        z1 = wq[i] * expl(q);
        aq[i] = z1 / er;
        Z += z1;
        U += aq[i] * eq[i];
    }

    printf("粒子分配函数Z为：%g\n", Z);
    double N = Z / er;
    printf("系统总的粒子数为N：%g\n", N);
    printf("系统总内能为U：%g\n", U);

    printf("请输入转动特征温度TV\n");
    while (1) {
        if (scanf("%f", &TV) != 1) {
            fprintf(stderr, "输入错误，请输入一个有效的转动特征温度。\n");
            continue;
        }
        break;
    }

    Cr = k0 * T;
    Ur = Cr * N;
    Ct = 1.5 * k0 * T;
    Ut = Ct * N;
    double m = TV / T;
    double F = -Cr * N * logl(Z);
    printf("系统自由能F为：%g\n", F);

    if (m != 0) {
        Cv = k0 * N * m * m * exp(m) / ((exp(m) - 1) * (exp(m) - 1));
        double n0 = N * k0 * TV;
        Uv = n0 / 2 + n0 / (exp(m) - 1);
    } else {
        Cv = 0;
        Uv = 0;
    }

    printf("平动定体热容Ct为：%g\n", Ct);
    printf("转动定体热容Cv为：%g\n", Cv);
    printf("振动定体热容Cr为：%g\n", Cr);
    printf("平动能量Ut为：%g\n", Ut);
    printf("转动能量Uv为：%g\n", Uv);
    printf("振动能量Ur为：%g\n", Ur);
    UV0 = Ut + Uv + Ur;
    C0 = Ct + Cv + Cr;
    printf("定体热容为：%g\n", C0);
    printf("系统内能为：%g\n", UV0);
}
