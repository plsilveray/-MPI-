#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include<time.h>  
#include <cstdlib>
#include<conio.h>
#include"10个函数.h"


#define POPSIZE 100             
#define MAXGENS 1000                                            
#define PXOVER 0.5               
#define PMUTATION 0.045         
#define TRUE 1
#define FALSE 0
struct genotype
{
	double gene[NVARS];                        //取维度为NVARS的变量进行优化
	double genefixed[NVARS];                    //另一半维度的变量固定不变
	double fitness;
	double upper[NVARS];
	double lower[NVARS];
	double rfitness;
	double cfitness;
};
class pop {
public:
	int generation;
	int cur_best;
	int flag;
	struct genotype population[POPSIZE + 1];
	struct genotype newpopulation[POPSIZE + 1];


	/*void initialize(void);
	double randval(double, double);
	void evaluate(void);
	void keep_the_best(void);
	void elitist(void);
	void select(void);
	void crossover(void);
	void Xover(int, int);
	void swap(double *, double *);
	void mutate(void);*/

	void initialize(void)
	{


		//	printf("initialize ......\n");
		int i, j;
		double lbound, ubound;
		//	fscanf(infile, "%lf", &lbound);
		//	fscanf(infile, "%lf", &ubound);
		lbound = -500;                                        	//	修改变量的取值范围。 
		ubound = 500;
		for (j = 0; j < POPSIZE + 1; j++)
		{
			population[j].fitness = 100;
			population[j].rfitness = 0;
			population[j].cfitness = 0;
			for (i = 0; i < NVARS; i++) {

				population[j].lower[i] = lbound;
				population[j].upper[i] = ubound;
				population[j].gene[i] = randval(population[j].lower[i],
					population[j].upper[i]);
				population[j].genefixed[i] = randval(population[j].lower[i],
					population[j].upper[i]);

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


	void evaluate(void)                                                    /*   在这里进行函数修改公式      ,当输入vector时加入函数参数为vector                */
	{
		int j;
		int i;
		double rt1[POPSIZE];
		double x[NVARS * 2];
		/*	if (flag == 1) {
		for (j = 0; j < POPSIZE; j++)
		{
		rt1[j] = 0;
		for (i = 0; i < NVARS; i++)
		{
		x[i] = population[j].gene[i];
		x[i + NVARS] = population[j].genefixed[i];
		}
		rt1[j] = f8(x);               //这里调用别的文件里的函数，输入x[i],计算rt1[j]
		}
		}*/
		//*if (flag == 0) {
		for (j = 0; j < POPSIZE; j++)
		{
			rt1[j] = 0;
			for (i = 0; i < NVARS; i++)
			{
				x[i] = 0;
				x[i + NVARS] = 2;
			}

			rt1[j] = f1(x);               //这里调用别的文件里的函数，输入x[i],计算rt1[j]
		}
		//	}
		for (j = 0; j < POPSIZE; j++)
			population[j].fitness = rt1[j];
		printf("%d\t%g\t%g\t%g\t%g\n", j, population[j].gene[0], population[j].genefixed[1], population[j].gene[2], population[j].fitness);

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
			population[POPSIZE].genefixed[i] = population[cur_best].genefixed[i];
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
			population[POPSIZE].genefixed[i] = population[best_mem].genefixed[i];
			population[POPSIZE].fitness = population[best_mem].fitness;
		}
		else
		{
			for (i = 0; i < NVARS; i++)
				population[worst_mem].gene[i] = population[POPSIZE].gene[i];
			population[POPSIZE].genefixed[i] = population[best_mem].genefixed[i];
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
		double lbound, hbound;
		double x;

		//	printf("mutate ......\n");

		for (i = 0; i < POPSIZE; i++) {
			for (j = 0; j < NVARS; j++)
			{
				x = rand() % 1000 / 1000.0;
				if (x < PMUTATION)
				{
					lbound = population[i].lower[j];
					hbound = population[i].upper[j];
					population[i].gene[j] = randval(lbound, hbound);
				}
			}
			//		printf("%d\t%g\t%g\t%g\n",i,population[i].gene[0],population[i].gene[1],population[i].gene[2]);
		}
	}


	/*void report(void)
	{
	int i;
	double best_val;
	double avg;
	double stddev;
	double sum_square;
	double square_sum;
	double sum;

	sum = 0.0;
	sum_square = 0.0;

	for (i = 0; i < POPSIZE; i++)
	{
	sum += population[i].fitness;
	sum_square += population[i].fitness * population[i].fitness;
	}

	avg = sum / (double)POPSIZE;
	square_sum = avg * avg * POPSIZE;
	stddev = sqrt((sum_square - square_sum) / (POPSIZE - 1));
	best_val = population[POPSIZE].fitness;


	}*/

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
			printf("best-val = %g\n", population[POPSIZE].fitness);


		}

		//printf( "best fitness = %g\n",  population[POPSIZE].fitness);	
	}
};
struct f_s                                     //这个结构体是最优解合并的f(s),fitness函数计算它的fitness
{
	double s1[NVARS];                          //s1和s2分别由pop1和pop2两个种群得到
	double s2[NVARS];
	double x[NVARS * 2];
	double fs;
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
	void sorting()                                  //冒泡排序函数，以及用数组a【】给出排序后的个体序号
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
void main()
{
	srand(time(NULL));
	int i, j;
	int generation = 0;
	pop pop1;
	pop pop2;
	f_s f_s;
	sort sort;
	pop1.initialize();
	pop2.initialize();
	for (int n = 0; n < POPSIZE; n++) {
		for (int m = 0; m < NVARS; m++) {
			pop1.population->gene[m] = pop2.population->gene[m];
			pop1.population->genefixed[m] = pop2.population->genefixed[m];
		}
	}
	pop1.flag = 1;
	pop2.flag = 0;
	while (generation < 1)
	{
		//	pop1.op();
		//printf("haha");
		pop2.op();
		for (i = 0; i < NVARS; i++) {
			f_s.s1[i] = pop1.population[POPSIZE].gene[i];
			f_s.s2[i] = pop2.population[POPSIZE].gene[i];
		}
		f_s.fitness();
		if (f_s.f > pop1.population[POPSIZE].fitness) {                //找pop1、pop2、f(s)三者最小fittness，并取代f（s）
			for (i = 0; i < NVARS; i++) {
				f_s.s1[i] = pop1.population[POPSIZE].gene[i];
				f_s.s2[i] = pop1.population[POPSIZE].gene[i];
			}
		}
		f_s.fitness();
		if (f_s.f > pop2.population[POPSIZE].fitness) {
			for (i = 0; i < NVARS; i++) {
				f_s.s1[i] = pop2.population[POPSIZE].gene[i];
				f_s.s2[i] = pop2.population[POPSIZE].gene[i];
			}
		}
		f_s.fitness();
		if (pop1.population[POPSIZE].fitness > pop2.population[POPSIZE].fitness)
			f_s.fs = pop2.population[POPSIZE].fitness;
		else
			f_s.fs = pop1.population[POPSIZE].fitness;
		if (f_s.f < f_s.fs)
			f_s.fs = f_s.f;
		printf("第%d次循环后f(s)的fitness：%g\n", generation, f_s.fs);
		for (i = 0; i < POPSIZE; i++)                               //固定gene部分的替代
			for (j = 0; j < NVARS; j++)
				pop1.population[i].genefixed[j] = f_s.s2[j];
		for (i = 0; i < POPSIZE; i++)
			for (j = 0; j < NVARS; j++)
				pop2.population[i].genefixed[j] = f_s.s1[j];
		for (i = 0; i < POPSIZE; i++)pop1.population[i].fitness = sort.a[i];   //找出pop中最大的50个（暂定为50）个体然后用f_s的s取代
		sort.sorting();
		for (i = 0; i < POPSIZE; i++)
			if (i = sort.sequence[i])
			{
				for (j = 0; j < NVARS; j++)
					pop1.population[i].gene[j] = f_s.s1[j];
			}
		for (i = 0; i < POPSIZE; i++)pop2.population[i].fitness = sort.a[i];
		sort.sorting();
		for (i = 0; i < POPSIZE; i++)
			if (i = sort.sequence[i])
			{
				for (j = 0; j < NVARS; j++)
					pop2.population[i].gene[j] = f_s.s2[j];
			}
		generation++;
	}
	for (i = 0; i < NVARS; i++) {
		f_s.s1[i] = pop1.population[POPSIZE].gene[i];
		f_s.s2[i] = pop2.population[POPSIZE].gene[i];
	}
	f_s.fitness();
	printf("最终的finess：%g", f_s.fs);
	getchar();
}