#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include<time.h>  
#include <cstdlib>
#include<conio.h>
#include "十个函数s.h"

/*
测试函数的执行代码在“源.cpp”文件中

*/
using namespace std;
#define POPSIZE 100                //一个种群里的个体数
#define MAXGENS 500                      //程序迭代的次数，或者说种群pop1和种群pop2协同进化的次数                   
#define PXOVER 0.5               
#define PMUTATION 0.045         
#define TRUE 1
#define FALSE 0
#define NVARS  15                 //一个个体有30维度的变量，取其中一半的变量即15个维度来进化，
#define upper 500               //每个变量的取值范围的最大值和最小值
#define lower -500
struct genotype
{
	double gene[NVARS];                        //取维度为NVARS的变量进行优化
	double genefixed[NVARS];                    //另一半维度的变量固定不变
	double fitness;
	
	
};
//pop类实现一个种群的功能
class pop {
public:
	int generation;
	int cur_best;
	int flag;
	struct genotype population[POPSIZE + 1];
	struct genotype newpopulation[POPSIZE + 1];




	void initialize(void)
	{


		//	printf("initialize ......\n");
		int i, j;
	
		//	fscanf(infile, "%lf", &lbound);
		//	fscanf(infile, "%lf", &ubound);
	
		for (j = 0; j < POPSIZE + 1; j++)
		{
			population[j].fitness = 100;
			
			for (i = 0; i < NVARS; i++) {

				
				population[j].gene[i] = randval(lower,
					upper);
				population[j].genefixed[i] = randval(lower,
					upper);

			}

		}

		population[POPSIZE].fitness = 10000000;


	}


	double randval(double low, double high)
	{
		double val;
		val = ((double)(rand() % 1000) / 1000.0)*(high - low) + low;
		return(val);
	}

	//evaluate函数有两个分支，一个flag等于1是计算种群pop1的fitness，flag等于2即计算种群Pop2的fitness
	void evaluate(void)                                                   
	{
		int j;
		int i;
		double rt1[POPSIZE];
		double x[NVARS * 2];
		if (flag == 1) {
			for (j = 0; j < POPSIZE; j++)
			{
				rt1[j] = 0;
				for (i = 0; i < NVARS; i++)
				{
					x[i] = population[j].gene[i];
					x[i + NVARS] = population[j].genefixed[i];
				}
				rt1[j] = f8(x);               //这里调用别的文件里的函数，输入x[i],计算种群Pop1的fitness
			}
		}
		if (flag == 0) {
			for (j = 0; j < POPSIZE; j++)
			{
				rt1[j] = 0;
				for (i = 0; i < NVARS; i++)
				{
					x[i] = population[j].genefixed[i];
					x[i + NVARS] = population[j].gene[i];
				}

				rt1[j] = f8(x);               //这里调用别的文件里的函数，输入x[i],计算种群pop2的fitness
				
			}
		}
		for (j = 0; j < POPSIZE; j++)
			population[j].fitness = rt1[j];
		//printf("%d\t%g\t%g\t%g\t%g\n", j, population[j].gene[0], population[j].genefixed[1], population[j].gene[2], population[j].fitness);

	}



	void keep_the_best()
	{
		int mem;
		int i;
		cur_best = 0;

		for (mem = 0; mem < POPSIZE; mem++)
		{
			if (population[mem].fitness < population[POPSIZE].fitness)
			{
				cur_best = mem;
				population[POPSIZE].fitness = population[mem].fitness;
			}
		}

		for (i = 0; i < NVARS; i++) {
			population[POPSIZE].gene[i] = population[cur_best].gene[i];
		}
	}

	void elitist()                              //遗传变异后找到最好和最坏的fitness
	{
		int i;
		double best, worst;
		int best_mem = 0, worst_mem = 0;

		best = population[0].fitness;
		worst = population[0].fitness;
		for (i = 0; i < POPSIZE - 1; ++i)
		{
			if (population[i].fitness < population[i + 1].fitness)
			{
				if (population[i].fitness <= best)
				{
					best = population[i].fitness;
					best_mem = i;
				}
				if (population[i + 1].fitness >= worst)
				{
					worst = population[i + 1].fitness;
					worst_mem = i + 1;
				}
			}
			else
			{
				if (population[i].fitness >= worst)
				{
					worst = population[i].fitness;
					worst_mem = i;
				}
				if (population[i + 1].fitness <= best)
				{
					best = population[i + 1].fitness;
					best_mem = i + 1;
				}
			}
		}

		if (best <= population[POPSIZE].fitness)
		{
			for (i = 0; i < NVARS; i++)
				population[POPSIZE].gene[i] = population[best_mem].gene[i];
			population[POPSIZE].fitness = population[best_mem].fitness;
		}
		else
		{
			for (i = 0; i < NVARS; i++)
				population[worst_mem].gene[i] = population[POPSIZE].gene[i];
			population[worst_mem].fitness = population[POPSIZE].fitness;
		}

		//printf("best individual %d\t fitness = %g\t%g\n",best_mem,population[POPSIZE].fitness,population[best_mem].fitness);
	}

	void select(void)                    //这里使用的是锦标赛法选择下一代个体            
	{
		int mem, i, j, k;
		struct genotype a[POPSIZE / 10];
		struct genotype b;
		for (k = 0; k<POPSIZE; k++)
		{

			for (i = 0; i<POPSIZE / 10; i++)
			{
				j = ((int)(rand() % (POPSIZE - 0)));
				a[i] = population[j];

			}
			for (i = 0; i<POPSIZE / 10; i++)
			{
				b = a[0];
				if (a[i].fitness<b.fitness)
					b = a[i];
			}
			newpopulation[k] = b;
		}
		for (i = 0; i<POPSIZE; i++)

			population[i] = newpopulation[i];

	}


	void crossover(void)
	{
		int i, mem, one;
		int first = 0;
		double x;

		//	printf("crossover ......\n");

		for (mem = 0; mem < POPSIZE; ++mem)
		{
			x = rand() % 1000 / 1000.0;
			if (x < PXOVER)
			{
				++first;
				if (first % 2 == 0)
					Xover(one, mem);
				else
					one = mem;
			}
			//		printf("%d\t%g\t%g\t%g\n",mem,population[mem].gene[0],population[mem].gene[1],population[mem].gene[2]);
		}
	}

	void Xover(int one, int two)
	{
		int i;
		int point;


		if (NVARS > 1)
		{
			if (NVARS == 2)
				point = 1;
			else
				point = (rand() % (NVARS - 1)) + 1;

			for (i = 0; i < point; i++)
				swap(&population[one].gene[i], &population[two].gene[i]);

		}
	}


	void swap(double *x, double *y)
	{
		double temp;

		temp = *x;
		*x = *y;
		*y = temp;

	}


	void mutate(void)
	{
		int i, j;
		double x;

		//	printf("mutate ......\n");

		for (i = 0; i < POPSIZE; i++) {
			for (j = 0; j < NVARS; j++)
			{
				x = rand() % 1000 / 1000.0;
				if (x < PMUTATION)
				{
					population[i].gene[j] = randval(lower, upper);
				}
			}
			//		printf("%d\t%g\t%g\t%g\n",i,population[i].gene[0],population[i].gene[1],population[i].gene[2]);
		}
	}



	void op()                              //遗传算法优化函数
	{
		generation = 0;
		//initialize();
		evaluate();
		keep_the_best();
		while (generation<MAXGENS)
		{
			generation++;
			select();
			crossover();
			mutate();
			//	report();
			evaluate();
			elitist();
			//	printf("result of generation %d",generation);
			/*	for (i = 0; i < NVARS; i++)
			{
			printf( "var(%d) = %g\t", i, population[POPSIZE].gene[i]);
			}
			*/
			//	printf("best-val = %g\n", population[POPSIZE].fitness);


		}

		//printf( "best fitness = %g\n",  population[POPSIZE].fitness);	
	}
};
struct f_s                                     //这个结构体是最优解合并的f(s),fitness函数计算它的fitness
{
	double s1[NVARS];                          //s1和s2是两个种群的一半维度的解，分别由pop1和pop2两个种群得到
	double s2[NVARS];
	double x[NVARS * 2];
	double fs = 1000000;
	double f;
	void fitness() {
		int i;
		f = 0;
		for (i = 0; i < NVARS; i++) {
			x[i] = s1[i];
			x[i + NVARS] = s2[i];
		}                                         //s1和s2组合成x【i】，调用函数计算fitness
		f = f8(x);
	}
};
struct sort                                     //用冒泡排序计算最大fitness的50个个体
{
	int sequence[POPSIZE];
	double a[POPSIZE];
	double temp;
	int i, j;
	void sorting()                                  //冒泡排序函数，以及用数组a[]给出排序后的个体序号
	{
		for (i = 0; i < POPSIZE; i++)sequence[i] = i;
		for (j = 0; j < POPSIZE - 1; j++)
			for (i = 0; i < POPSIZE - 1 - j; i++)
			{
				if (a[i] < a[i + 1])
				{
					temp = a[i];
					a[i] = a[i + 1];
					a[i + 1] = temp;
					temp = sequence[i];
					sequence[i] = sequence[i + 1];
					sequence[i + 1] = temp;
				}
			}
	}
};
int main()
{
	clock_t startTime, endTime;
	startTime = clock();
	srand(time(NULL));
	int i, j;
	int generation = 0;
	pop pop1;
	pop pop2;
	f_s f_s;
	sort sort;
	pop1.initialize();
	pop2 = pop1;

	pop1.flag = 1;
	pop2.flag = 0;
	while (generation < 50)             //程序迭代次数为50
	{
		pop1.op();
		//printf("pop1:%d,%g\n",generation,pop1.population[POPSIZE].fitness);
		pop2.op();
		//printf("pop2:%d,%g\n", generation, pop2.population[POPSIZE].fitness);
		for (i = 0; i < NVARS; i++) {
			f_s.s1[i] = pop1.population[POPSIZE].gene[i];
			f_s.s2[i] = pop2.population[POPSIZE].gene[i];
		}
		f_s.fitness();
		//找pop1、pop2、f(s)三者最小fittness，并取代f（s）
		if (pop2.population[POPSIZE].fitness > pop1.population[POPSIZE].fitness &&pop1.population[POPSIZE].fitness<f_s.f) 
		{
			for (i = 0; i < NVARS; i++) {
				f_s.s1[i] = pop1.population[POPSIZE].gene[i];
				f_s.s2[i] = pop1.population[POPSIZE].genefixed[i];
			}
			f_s.f = pop1.population[POPSIZE].fitness;
		}
		
		if (pop1.population[POPSIZE].fitness > pop2.population[POPSIZE].fitness &&pop2.population[POPSIZE].fitness<f_s.f) 
		{
			for (i = 0; i < NVARS; i++) {
				f_s.s1[i] = pop2.population[POPSIZE].genefixed[i];
				f_s.s2[i] = pop2.population[POPSIZE].gene[i];
			}
			f_s.f = pop2.population[POPSIZE].fitness;
		}

		if (f_s.f < f_s.fs)
			f_s.fs = f_s.f;
		printf("第%d次循环后f(s)的fitness：%g\n", generation, f_s.fs);
		for (i = 0; i < POPSIZE; i++)                               //用得到的新解fs去固定gene部分的替代
			for (j = 0; j < NVARS; j++)
				pop1.population[i].genefixed[j] = f_s.s2[j];
		for (i = 0; i < POPSIZE; i++)
			for (j = 0; j < NVARS; j++)
				pop2.population[i].genefixed[j] = f_s.s1[j];
		for (i = 0; i < POPSIZE; i++)sort.a[i]=pop1.population[i].fitness ;   //找出pop中最大的50个（暂定为50）个体然后用f_s的s取代
		sort.sorting();
		for (i = 0; i < POPSIZE / 2; i++)
		{int k = sort.sequence[i];
		//printf("k:%d\n", k);
				for (j = 0; j < NVARS; j++)
					pop1.population[k].gene[j] = f_s.s1[j];
			}
		for (i = 0; i < POPSIZE; i++) sort.a[i]= pop2.population[i].fitness ;
		sort.sorting();
		for (i = 0; i < POPSIZE/2; i++)
		{
			int k = sort.sequence[i];
				for (j = 0; j < NVARS; j++)
					pop2.population[i].gene[j] = f_s.s2[j];
			}
		generation++;
	}
	
	printf("最终的finess：%g\n", f_s.fs);
	endTime = clock();
	cout << "Totle Time : " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
	getchar();
	return 0;
}