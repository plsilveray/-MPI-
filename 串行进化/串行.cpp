#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include<time.h>  
#include <cstdlib>
#include<conio.h>
#include "ʮ������s.h"

/*
���Ժ�����ִ�д����ڡ�Դ.cpp���ļ���

*/
using namespace std;
#define POPSIZE 100                //һ����Ⱥ��ĸ�����
#define MAXGENS 500                      //��������Ĵ���������˵��Ⱥpop1����Ⱥpop2Эͬ�����Ĵ���                   
#define PXOVER 0.5               
#define PMUTATION 0.045         
#define TRUE 1
#define FALSE 0
#define NVARS  15                 //һ��������30ά�ȵı�����ȡ����һ��ı�����15��ά����������
#define upper 500               //ÿ��������ȡֵ��Χ�����ֵ����Сֵ
#define lower -500
struct genotype
{
	double gene[NVARS];                        //ȡά��ΪNVARS�ı��������Ż�
	double genefixed[NVARS];                    //��һ��ά�ȵı����̶�����
	double fitness;
	
	
};
//pop��ʵ��һ����Ⱥ�Ĺ���
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

	//evaluate������������֧��һ��flag����1�Ǽ�����Ⱥpop1��fitness��flag����2��������ȺPop2��fitness
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
				rt1[j] = f8(x);               //������ñ���ļ���ĺ���������x[i],������ȺPop1��fitness
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

				rt1[j] = f8(x);               //������ñ���ļ���ĺ���������x[i],������Ⱥpop2��fitness
				
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

	void elitist()                              //�Ŵ�������ҵ���ú����fitness
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

	void select(void)                    //����ʹ�õ��ǽ�������ѡ����һ������            
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



	void op()                              //�Ŵ��㷨�Ż�����
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
struct f_s                                     //����ṹ�������Ž�ϲ���f(s),fitness������������fitness
{
	double s1[NVARS];                          //s1��s2��������Ⱥ��һ��ά�ȵĽ⣬�ֱ���pop1��pop2������Ⱥ�õ�
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
		}                                         //s1��s2��ϳ�x��i�������ú�������fitness
		f = f8(x);
	}
};
struct sort                                     //��ð������������fitness��50������
{
	int sequence[POPSIZE];
	double a[POPSIZE];
	double temp;
	int i, j;
	void sorting()                                  //ð�����������Լ�������a[]���������ĸ������
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
	while (generation < 50)             //�����������Ϊ50
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
		//��pop1��pop2��f(s)������Сfittness����ȡ��f��s��
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
		printf("��%d��ѭ����f(s)��fitness��%g\n", generation, f_s.fs);
		for (i = 0; i < POPSIZE; i++)                               //�õõ����½�fsȥ�̶�gene���ֵ����
			for (j = 0; j < NVARS; j++)
				pop1.population[i].genefixed[j] = f_s.s2[j];
		for (i = 0; i < POPSIZE; i++)
			for (j = 0; j < NVARS; j++)
				pop2.population[i].genefixed[j] = f_s.s1[j];
		for (i = 0; i < POPSIZE; i++)sort.a[i]=pop1.population[i].fitness ;   //�ҳ�pop������50�����ݶ�Ϊ50������Ȼ����f_s��sȡ��
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
	
	printf("���յ�finess��%g\n", f_s.fs);
	endTime = clock();
	cout << "Totle Time : " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
	getchar();
	return 0;
}