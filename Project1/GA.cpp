#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include<time.h>  
#include <cstdlib>
#include<conio.h>
#include "10个函数的计算公式.h"

#define POPSIZE 100             
#define MAXGENS 100000            
                                   /*   在这里进行函数参数个数的修改                  */
#define PXOVER 0.5               
#define PMUTATION 0.045         
#define TRUE 1
#define FALSE 0

int generation;
int cur_best;

struct genotype
{
	double gene[NVARS];
	double fitness;
	double upper[NVARS];
	double lower[NVARS];
	double rfitness;
	double cfitness;
};

struct genotype population[POPSIZE + 1];
struct genotype newpopulation[POPSIZE + 1];


void initialize(void);
double randval(double, double);
void evaluate(void);
void keep_the_best(void);
void elitist(void);
void select(void);
void crossover(void);
void Xover(int, int);
void swap(double *, double *);
void mutate(void);
//void report(void);


void initialize(void)
{
	

	//	printf("initialize ......\n");
	int i, j;
	double lbound, ubound;
	//	fscanf(infile, "%lf", &lbound);
	//	fscanf(infile, "%lf", &ubound);
	lbound = -100;                                        	//	修改变量的取值范围。 
	ubound = 100;
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


void evaluate(void)                                        /*   在这里进行函数修改                  */
{
	int j;
	int i;
	double rt1[POPSIZE];
	double x[NVARS * 2];

	for (j = 0; j < POPSIZE; j++)
	{
		rt1[j] = 0;
		for (i = 0; i < NVARS; i++)
		{
			x[i] = population[j].gene[i];
			
		}
		rt1[j] = f4(x);               //这里调用别的文件里的函数，输入x[i],计算rt1[j]
	}
	for (j = 0; j < POPSIZE; j++)
		population[j].fitness = rt1[j];
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

	for (i = 0; i < NVARS; i++)
		population[POPSIZE].gene[i] = population[cur_best].gene[i];
}

void elitist()
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

	//	printf("best individual %d\t fitness = %g\t%g\n",best_mem,population[POPSIZE].fitness,population[best_mem].fitness);
}

void select(void)
{
	int mem, i, j, k;
	struct genotype a[10];
	struct genotype b;
	for (k = 0; k<100; k++)
	{

		for (i = 0; i<10; i++)
		{
			j = ((int)(rand() % (POPSIZE - 0)));
			a[i] = population[j];

		}
		for (i = 0; i<10; i++)
		{
			b = a[0];
			if (a[i].fitness<b.fitness)
				b = a[i];
		}
		newpopulation[k] = b;
	}
	for (i = 0; i<100; i++)

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

int main()
{
	srand(time(NULL));
	int i, j, k;
	double a[5], sum, avg;
	generation = 0;
	initialize();
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
		printf( "best-val = %g\n", population[POPSIZE].fitness);
		*/

	}

	printf("best fitness = %g\n", population[POPSIZE].fitness);
	getchar();
	return 0;
}
