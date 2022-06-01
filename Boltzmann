/*****************************************************************************80
    Copyright 2016 Songsong Wang <sswang@hit.edu.cn>
/*----------------------------------------------------------------------------80
    Author: Songsong Wang <sswang@hit.edu.cn>
    Date: 2016.01.04
------------------------------------------------------------------------------*/

#include<stdio.h>
#include<math.h>
#include<time.h>
#define size  20
#define k0   1.380650e-23                                        //波尔兹曼常量
double Fact(int n);                                              //本程序多次用到阶乘，故定义用递归法求阶乘    
void Boltzmann_System(int a[],int w[],int l);                    //玻尔兹曼系统微观状态数求解
void Bose_System(int a[],int w[],int l);                         //波色系统求解微观状态数
void Fermi_System(int a[],int w[],int l);                        //费米系统求解微观状态数
Boltzmann();                                                     //基于波尔兹曼统计的相关求解

void main()
{   	
	int a[size]={0};   //用数组a[]表示al的大小
	int w[size]={0};   //数组w[]表示w的值
	int k,i=0;
	printf("请输入能级数l\n");
	scanf("%d",&k);    //输入能级数
	printf("请依次输入l个各能级上的粒子数\n");                   //若输入的数字大于l个，将自动赋值给w[]
	for(i=0;i<k;i++)   
	scanf("%d",&a[i]);
	printf("请依次输入l个各能级上的简并度\n");
	for(i=0;i<k;i++)
	scanf("%d",&w[i]);
    Boltzmann_System(a,w,k);
	Bose_System(a,w,k);
	Fermi_System(a,w,k);
    Boltzmann();
}

double Fact(int temp)
{
	double n=temp;
	if(n<0)
		return -1;             //处理非法数据 
	else if(n==0 || n==1)      //基线情况，即递归终止条件
		return 1;
	else                       //一般情况
		return (n*Fact(n-1));  //调用递归，利用(n-1)!计算n!
}

void Boltzmann_System(int a[],int w[],int l)                       //玻尔兹曼系统微观状态数求解
{
	
	int n=0,i=0,j=0;
	int a1=0,w1=0;
	long double N=1,W=1.0,A=1,SUM=0;                           //N表示各能级粒子数的阶乘，W表示兼并度的粒子数次方的连乘，A表示粒子数的截成的连乘，SUM表示总的微观状态数
	
	for(i=0;i<l;i++)
	n=n+a[i];

	N=Fact(n);

	for(i=l-1;i>=0;i--)
	{
	    a1=a[i];
       	    w1=w[i];	
	    W=pow(w1,a1)*W;              
	    A*=Fact(a1);
	}	
	SUM=N*W/A;
	printf("波尔兹曼系统总的微观状态数为：%g\n",SUM);
}

void Bose_System(int a[],int w[],int l)
{   
	int n=0,i=0;
	long a1=0,w1=0,aw=0;
	long double N=1,W=1,A=1,WA=1;
	long double SUM=1;
	
	for(i=0;i<l-1;i++)
	{
	    n=n+a[i];
    }
	N=Fact(n);
	for(i=l-1;i>=0;i--)
	{
	    a1=a[i];
            w1=w[i]-1;
            aw=a1+w1;
            W=Fact(w1);
	    WA=Fact(aw);
	    A=Fact(a1);
            SUM=SUM*WA/(A*W);
	}	
	printf("波色系统总微观状态数%g\n",SUM);
}
void Fermi_System(int a[],int w[],int l)
{ 
	int n=0,i=0;
	long a1=0,w1=0,aw=0;
	double N=1,W=1,A=1,WA=1;
	double SUM=1;	
	for(i=0;i<l-1;i++)
	{
	    n=n+a[i];
    }
	while(n>0)
	{
	    N=N*n;
            n=n-1;
	}
	for(i=l-1;i>=0;i--)
	{
	    a1=a[i];
	    w1=w[i];
	    aw=w1-a1;
	    W=Fact(w1);
	    WA=Fact(aw);
	    A=Fact(a1);
            SUM=SUM*W/(A*WA);
	}	
	printf("费米系统总微观状态数%g\n",SUM);
}
                     
Boltzmann()
{
    int wq[size]={0};
    double aq[size]={0};
    long double eq[size]={0};
    int i,w1,r;                        
    float T=0,T0=20.3,t1=0,TV=0;                      //T表示输入的温度，T0为氢气在1pn时的沸点（单位为开尔文）,TV振动特征温度
    double e1,er;                                      //er表示相应T值对应的e的阿尔法次方的值
    long double a1=0,q=0,p=0,z1=0,Z=0,N=0,u1=0,U=0,N2,F=0,UV0=0,C0=0,Ut=0,Uv=0,Ur=0,Cv=0,Ct=0,Cr=0,m=0,n0=0;   //F表示自由能,Ct表示平动定体热容，Cv表示振动定体热容,Cr振动定体热容，Ut表示平动内能，Uv表示转动动能，Ur表示振动内能
	                                                                                     //UV0表示理想双原子气体的内能,C0表示定体热容
    srand((unsigned)time(NULL));                       //播种子
    printf("请输入不小于224的温度（单位为K）T：\n");                           //本程序中T的最小值为224，,否则因数据溢出无法得到结果
    scanf("%f",&T);
    p=k0*T;
    p=1/p;
    t1=T/T0;
    er=140*pow(t1,1.5);
    for(i=0;i<size;i++)
    {
        r=(i+1)*(i+1);                 //计算l的平方
	eq[i]=-13.6*(1.60217e-19)/r;    //氢原子各能级对应能量的值
        wq[i]=rand()%10;                //产生100以内的随机整数
	w1=wq[i];
	e1=eq[i];
	q=-p*e1;
	z1=w1*exp(q);
	aq[i]=z1/er;
	a1=aq[i];
	Z+=z1;
	u1=a1*e1;
	U+=u1;
	}
	printf("粒子分配函数Z为：%g\n",Z);
    N=Z/er;
	printf("系统总的粒子数为N：%g\n",N);              //N为系统粒子总数
	printf("系统总内能为U：%g\n",U);                  //输出总内能
    printf("请输入转动特征温度TV\n");
	scanf("%f",&TV);
	Cr=k0*T;
	Ur=Cr*N;
	Ct=1.5*k0*T;
    Ut=Ct*N;
	m=TV/T;
    F=-Cr*N*log(Z);
	printf("系统自由能F为：%g\n",F);
	if(m!=0)
	{
	    Cv=k0*N*m*m*exp(m)/((exp(m)-1)*(exp(m)-1));
	    n0=N*k0*TV;
	}
	Uv=n0/2+n0/(exp(m)-1);
	printf("平动定体热容Ct为：%g\n",Ct);
	printf("转动定体热容Cv为：%g\n",Cv);
	printf("振动定体热容Cr为：%g\n",Cr);
	printf("平动能量Ut为：%g\n",Ut);
	printf("转动能量Uv为：%g\n",Uv);
	printf("振动能量Ur为：%g\n",Ur);
	UV0=Ut+Uv+Ur;
	C0=Ct+Cv+Cr;
	printf("定体热容为：%g\n",C0);
	printf("系统内能为：%g\n",UV0);
}
