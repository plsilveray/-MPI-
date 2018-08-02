#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include<time.h>  
#include <cstdlib>
#include<conio.h>
#include"10个函数的计算.h"
#include <mpi.h>
#pragma comment(lib,"msmpi.lib")

using namespace std;
#define POPSIZE 100             
#define MAXGENS 500                                         
#define PXOVER 0.5               
#define PMUTATION 0.045         
#define TRUE 1
#define FALSE 0
#define NVARS  15
#define lbound  -500                                       	//	修改变量的取值范围。 
#define ubound  500

struct genotype
{
	//另一半维度的变量固定不变
	double fitness;
	double gene[NVARS];                //取维度为NVARS的变量进行优化
	double genefixed[NVARS];

};
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
		for (j = 0; j < POPSIZE + 1; j++)
		{
			population[j].fitness = 100;
			for (i = 0; i < NVARS; i++) {


				population[j].gene[i] = randval(lbound,
					ubound);
				population[j].genefixed[i] = randval(lbound,
					ubound);

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
		if (flag == 1) {
			for (j = 0; j < POPSIZE; j++)
			{
				rt1[j] = 0;
				for (i = 0; i < NVARS; i++)
				{
					//	printf("\nevaluate里gene的值为:%g,%g\n", population[j].gene[i], population[j].genefixed[i]);
					x[i] = population[j].gene[i];
					x[i + NVARS] = population[j].genefixed[i];
					//全是0	printf("\nx的值为:%g,%g\n", x[i], x[i+NVARS]);
				}
				rt1[j] = f8(x);               //这里调用别的文件里的函数，输入x[i],计算rt1[j]
											  
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

				rt1[j] = f8(x);               //这里调用别的文件里的函数，输入x[i],计算rt1[j]

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

					population[i].gene[j] = randval(lbound, ubound);
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
		//printf("keepbest,%d,best fitness = %g\n", flag, population[POPSIZE].fitness);
		while (generation<MAXGENS)
		{
			generation++;
			select();
			//	printf("select,%d,best fitness = %g\n", flag, population[POPSIZE].fitness);
			crossover();
			//	printf("crossover,%d,best fitness = %g\n", flag, population[POPSIZE].fitness);
			mutate();
			//	printf("mutate,%d,best fitness = %g\n", flag, population[POPSIZE].fitness);
			evaluate();
			//printf("evaluate,%d,best fitness = %g\n", flag, population[POPSIZE].fitness);
			elitist();
			//printf("elitist,%d,best fitness = %g\n", flag, population[POPSIZE].fitness);
			//	printf("result of generation %d",generation);
			//	printf("best-val = %g\n", population[POPSIZE].fitness);
		}
		//printf("%d,best fitness = %g\n", flag,population[POPSIZE].fitness);
	}
};
struct f_s                                     //这个结构体是最优解合并的f(s),fitness函数计算它的fitness
{
	double s1[NVARS];                          //s1和s2分别由pop1和pop2两个种群得到
	double s2[NVARS];
	double fs = 1000000;
	double f;
	double x[NVARS * 2];
	void fitness() {
		int i;
		f = 0;
		for (i = 0; i < NVARS; i++) {
			x[i] = s1[i];
			x[i + NVARS] = s2[i];
		}                                         //s1和s2组合成x【i】，调用函数计算fitness
		f = f8(x);
		//printf("f:%g\n", f);
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
//int main(int argc, char* argv[])
//{
//	clock_t startTime, endTime;
//	startTime = clock();
//	srand(time(NULL));
//	int i, j;
//	//int generation = 0;
//	int rank, size;
//	pop pop0;
//	pop pop1;
//	pop pop2;
//	f_s f_s;
//	sort sort;
//	MPI_Init(&argc, &argv);
//	MPI_Datatype myvar;
//	//指定原来旧的数据类型
//	MPI_Datatype old_types[3] = { MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE };
//
//	//指定每个块中的变量个数,这里只有2个块，其中包含4个MPI_INT,6个MPI_FLOAT
//	int blockLength[] = { 1, NVARS,NVARS };
//	//addressOffsets是用来描述CTestSendRecv中各个成员离该CTestSendRecv类的首地址的偏移量，这里涉及到编译器的字节对齐的知识
//	MPI_Aint indices[3] = { 0,sizeof(double),(NVARS + 1) * sizeof(double) };
//	MPI_Type_create_struct(3, blockLength, indices, old_types, &myvar);
//	MPI_Type_commit(&myvar);
//	//////////////////////////////////////
//
//	MPI_Datatype myvar2;
//	//指定原来旧的数据类型
//	MPI_Datatype old_types2[4] = { MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE };
//
//	//指定每个块中的变量个数,这里只有2个块，其中包含4个MPI_INT,6个MPI_FLOAT
//	int blockLength2[] = { NVARS,NVARS,1,1 };
//	//addressOffsets是用来描述CTestSendRecv中各个成员离该CTestSendRecv类的首地址的偏移量，这里涉及到编译器的字节对齐的知识
//	MPI_Aint indices2[4] = { 0,(NVARS + 1) * sizeof(double),(2 * NVARS + 1) * sizeof(double), (2 * NVARS + 2) * sizeof(double) };
//	MPI_Type_create_struct(4, blockLength2, indices2, old_types2, &myvar2);
//	MPI_Type_commit(&myvar2);
//	///////////////////////////////////////
//	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//	MPI_Comm_size(MPI_COMM_WORLD, &size);
//	if (rank == 0) {
//		pop0.initialize();
//		pop1 = pop0;
//		pop2 = pop0;
//
//
//		MPI_Send(&pop1.population, POPSIZE + 1, myvar, 1, 0, MPI_COMM_WORLD);
//		MPI_Send(&pop2.population, POPSIZE + 1, myvar, 2, 1, MPI_COMM_WORLD);
//	}
//	if (rank % 2 == 1) {
//		if (rank == 1)//奇数进程pop1
//			MPI_Recv(&pop1.population, POPSIZE + 1, myvar, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//
//		if (rank>2)
//			MPI_Recv(&pop1.population, POPSIZE + 1, myvar, rank - 2, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//		pop1.flag = 0;
//		pop1.op();
//		//printf("pop1:%d,%g\n", rank, pop1.population[POPSIZE].fitness);
//
//		//发pop1给偶数进程计算fs
//		MPI_Send(&pop1.population, POPSIZE + 1, myvar, rank + 1, 2, MPI_COMM_WORLD);
//		//接受fs来形成新的pop1
//
//		MPI_Recv(&f_s, 1, myvar2, rank + 1, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);   //flag=4 !!1
//
//
//		for (i = 0; i < POPSIZE; i++)                               //固定gene部分的替代
//			for (j = 0; j < NVARS; j++)
//				pop1.population[i].genefixed[j] = f_s.s2[j];
//		for (i = 0; i < POPSIZE; i++)sort.a[i] = pop1.population[i].fitness;   //找出pop中最大的50个（暂定为50）个体然后用f_s的s取代
//		sort.sorting();
//		for (i = 0; i < POPSIZE / 2; i++)
//		{
//			int k = sort.sequence[i];
//			//printf("k:%d\n", k);
//			for (j = 0; j < NVARS; j++)
//				pop1.population[k].gene[j] = f_s.s1[j];
//		}
//		if (rank <= 97)
//			MPI_Send(&pop1.population, POPSIZE + 1, myvar, rank + 2, 3, MPI_COMM_WORLD);    //通信flag=3!！！！
//																							//得到新的pop1后发送给rank+2的进程开始下一次进化
//	}
//
//
//	if (rank % 2 == 0 && rank>0) {
//		if (rank == 2)//偶数2+
//			MPI_Recv(&pop2.population, POPSIZE + 1, myvar, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//
//		if (rank>3)
//			MPI_Recv(&pop2.population, POPSIZE + 1, myvar, rank - 2, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);   //rank-2存在吗！！！
//		pop2.flag = 1;
//		pop2.op();
//	//	printf("pop2:%d,%g\n", rank-1, pop2.population[POPSIZE].fitness);
//		MPI_Recv(&pop1.population, POPSIZE + 1, myvar, rank - 1, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//		for (i = 0; i < NVARS; i++) {
//			f_s.s1[i] = pop1.population[POPSIZE].gene[i];
//			f_s.s2[i] = pop2.population[POPSIZE].gene[i];
//			//printf("gene:%g,%g\n", pop1.population[POPSIZE].gene[i], pop2.population[POPSIZE].gene[i]);
//		}
//		f_s.fitness();
//		//	printf("f_s.f:%g,pop1andpop2:%g,%g\n", f_s.f, pop1.population[POPSIZE].fitness, pop2.population[POPSIZE].fitness);
//		f_s.fs = f_s.f;
//		//printf("第一个fs：%g", f_s.fs);
//		if (pop1.population[POPSIZE].fitness < pop2.population[POPSIZE].fitness)
//		{
//			if (pop1.population[POPSIZE].fitness < f_s.f)
//			{
//				f_s.fs = pop1.population[POPSIZE].fitness;
//				for (i = 0; i < NVARS; i++)
//				{
//					f_s.s1[i] = pop1.population[POPSIZE].gene[i];
//					f_s.s2[i] = pop1.population[POPSIZE].gene[i];
//				}
//			}
//		}
//		else
//		{
//			if (pop2.population[POPSIZE].fitness < f_s.f)
//			{
//				f_s.fs = pop2.population[POPSIZE].fitness;
//				for (i = 0; i < NVARS; i++)
//				{
//					f_s.s1[i] = pop2.population[POPSIZE].gene[i];
//					f_s.s2[i] = pop2.population[POPSIZE].gene[i];
//				}
//			}
//		}
//		if (rank == 100)
//			{
//        printf("结果：%g\n", f_s.fs);
//endTime = clock();
//                   cout << "Totle Time : " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
//            exit(0); }
//			//printf("第%d次循环后f(s)的fitness：%g\n", rank, f_s.fs);
//		MPI_Send(&f_s, 1, myvar2, rank - 1, 4, MPI_COMM_WORLD);
//		
//		for (i = 0; i < POPSIZE; i++)
//			for (j = 0; j < NVARS; j++)
//				pop2.population[i].genefixed[j] = f_s.s1[j];
//		for (i = 0; i < POPSIZE; i++)sort.a[i] = pop2.population[i].fitness;   //找出pop中最大的50个（暂定为50）个体然后用f_s的s取代
//		sort.sorting();
//		for (i = 0; i < POPSIZE / 2; i++)
//		{
//			int k = sort.sequence[i];
//			//printf("k:%d\n", k);
//			for (j = 0; j < NVARS; j++)
//				pop2.population[k].gene[j] = f_s.s2[j];
//		}
//		//计算完fs ,发送pop1给奇数进程来进化
//		if (rank <= 98)
//			MPI_Send(&pop2.population, POPSIZE + 1, myvar, rank + 2, 5, MPI_COMM_WORLD);
//	}
//	MPI_Finalize();
//	return 0;
//}

int main(int argc, char* argv[])
{
	clock_t startTime, endTime;
	startTime = clock();
	srand(time(NULL));
	int i, j;
	//int generation = 0;
	int rank, size;
	pop pop0;
	pop pop1;
	pop pop2;
	f_s f_s;
	sort sort;
	MPI_Init(&argc, &argv);
	int k = 1;
	MPI_Datatype myvar;
	//指定原来旧的数据类型
	MPI_Datatype old_types[3] = { MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE };

	//指定每个块中的变量个数,这里只有2个块，其中包含4个MPI_INT,6个MPI_FLOAT
	int blockLength[] = { 1, NVARS,NVARS };
	//addressOffsets是用来描述CTestSendRecv中各个成员离该CTestSendRecv类的首地址的偏移量，这里涉及到编译器的字节对齐的知识
	MPI_Aint indices[3] = { 0,sizeof(double),(NVARS + 1) * sizeof(double) };
	MPI_Type_create_struct(3, blockLength, indices, old_types, &myvar);
	MPI_Type_commit(&myvar);
	//////////////////////////////////////

	MPI_Datatype myvar2;
	//指定原来旧的数据类型
	MPI_Datatype old_types2[4] = { MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE };

	//指定每个块中的变量个数,这里只有2个块，其中包含4个MPI_INT,6个MPI_FLOAT
	int blockLength2[] = { NVARS,NVARS,1,1 };
	//addressOffsets是用来描述CTestSendRecv中各个成员离该CTestSendRecv类的首地址的偏移量，这里涉及到编译器的字节对齐的知识
	MPI_Aint indices2[4] = { 0,(NVARS + 1) * sizeof(double),(2 * NVARS + 1) * sizeof(double), (2 * NVARS + 2) * sizeof(double) };
	MPI_Type_create_struct(4, blockLength2, indices2, old_types2, &myvar2);
	MPI_Type_commit(&myvar2);
	
	///////////////////////////////////////
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	//线程0号，用于收集来自种群pop1和pop2的解得到新解fs并计算结果再传出去
	if (rank == 0) {
		pop0.initialize();
		pop1 = pop0;
		pop2 = pop0;
		//cout << "pop0的fitness" << pop0.population[POPSIZE].fitness << endl;
		//cout << "pop1的fitness" << pop1.population[POPSIZE].fitness << endl;
		MPI_Send(&pop1.population, POPSIZE + 1, myvar, 1, 0, MPI_COMM_WORLD);
		MPI_Send(&pop2.population, POPSIZE + 1, myvar, 2, 0, MPI_COMM_WORLD);
		while (k<51 ) {
			MPI_Recv(&pop1.population, POPSIZE + 1, myvar, 1,k, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(&pop2.population, POPSIZE + 1, myvar, 2, k, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			//printf("fitness:%g,%g\n", pop1.population[POPSIZE].fitness, pop2.population[POPSIZE].fitness);
			for (i = 0; i < NVARS; i++) {
				f_s.s1[i] = pop1.population[POPSIZE].gene[i];
				f_s.s2[i] = pop2.population[POPSIZE].gene[i];
				//printf("gene:%g,%g\n", pop1.population[POPSIZE].gene[i], pop2.population[POPSIZE].gene[i]);
			}
			f_s.fitness();
			//	printf("f_s.f:%g,pop1andpop2:%g,%g\n", f_s.f, pop1.population[POPSIZE].fitness, pop2.population[POPSIZE].fitness);
			f_s.fs = f_s.f;
			//printf("第一个fs：%g", f_s.fs);
			if (pop1.population[POPSIZE].fitness < pop2.population[POPSIZE].fitness)
			{
				if (pop1.population[POPSIZE].fitness < f_s.f)
				{
					f_s.fs = pop1.population[POPSIZE].fitness;
					for (i = 0; i < NVARS; i++)
					{
						f_s.s1[i] = pop1.population[POPSIZE].gene[i];
						f_s.s2[i] = pop1.population[POPSIZE].gene[i];
					}
				}
			}
			else
			{
				if (pop2.population[POPSIZE].fitness < f_s.f)
				{
					f_s.fs = pop2.population[POPSIZE].fitness;
					for (i = 0; i < NVARS; i++)
					{
						f_s.s1[i] = pop2.population[POPSIZE].gene[i];
						f_s.s2[i] = pop2.population[POPSIZE].gene[i];
					}
				}
			}
			printf("第%d次：%g\n",k, f_s.fs);
			if (k== 50)
			{
				printf("结果：%g\n", f_s.fs);
				endTime = clock();
				cout << "Totle Time : " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
				exit(0);
			}
			MPI_Send(&f_s, 1, myvar2,1, k, MPI_COMM_WORLD);
			MPI_Send(&f_s, 1, myvar2,2, k, MPI_COMM_WORLD);
			k++;
		}
	}
	//线程1号，用于处理种群pop1的所有动作
	if (rank  == 1) {
		
			MPI_Recv(&pop1.population, POPSIZE + 1, myvar, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		//printf("fitness:%g\n",  pop1.population[POPSIZE].fitness);
			//cout << "pop1其他的finess" << pop1.population[POPSIZE - 20].fitness << pop1.population[POPSIZE - 1].fitness << endl;
		pop1.flag = 0;
		while (k < 51) {
			pop1.op();
		//	printf("%d,fitness:%g\n",rank, pop1.population[POPSIZE].fitness);
			//printf("pop1:%d,%g\n", rank, pop1.population[POPSIZE].fitness);

			//发pop1给0进程计算fs
			MPI_Send(&pop1.population, POPSIZE + 1, myvar, 0, k, MPI_COMM_WORLD);
			//接受fs来形成新的pop1

			MPI_Recv(&f_s, 1, myvar2, 0, k, MPI_COMM_WORLD, MPI_STATUS_IGNORE);   


			for (i = 0; i < POPSIZE; i++)                               //固定gene部分的替代
				for (j = 0; j < NVARS; j++)
					pop1.population[i].genefixed[j] = f_s.s2[j];
			for (i = 0; i < POPSIZE; i++)sort.a[i] = pop1.population[i].fitness;   //找出pop中最大的50个（暂定为50）个体然后用f_s的s取代
			sort.sorting();
			for (i = 0; i < POPSIZE / 2; i++)
			{
				int k = sort.sequence[i];
				//printf("k:%d\n", k);
				for (j = 0; j < NVARS; j++)
					pop1.population[k].gene[j] = f_s.s1[j];
			}
			k++;
		}
	}
//线程2号，用于处理种群pop2的所有动作
	if (rank ==2) {

		MPI_Recv(&pop2.population, POPSIZE + 1, myvar, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		pop2.flag = 1;

		while (k < 51) {
			pop2.op();
			//printf("%d,fitness:%g\n", rank, pop2.population[POPSIZE].fitness);
			MPI_Send(&pop2.population, POPSIZE + 1, myvar, 0, k, MPI_COMM_WORLD);
			//接受fs来形成新的pop1

			MPI_Recv(&f_s, 1, myvar2, 0, k, MPI_COMM_WORLD, MPI_STATUS_IGNORE);   //flag=4 !!1
			for (i = 0; i < POPSIZE; i++)
				for (j = 0; j < NVARS; j++)
					pop2.population[i].genefixed[j] = f_s.s1[j];
			for (i = 0; i < POPSIZE; i++)sort.a[i] = pop2.population[i].fitness;   //找出pop中最大的50个（暂定为50）个体然后用f_s的s取代
			sort.sorting();
			for (i = 0; i < POPSIZE / 2; i++)
			{
				int k = sort.sequence[i];
				//printf("k:%d\n", k);
				for (j = 0; j < NVARS; j++)
					pop2.population[k].gene[j] = f_s.s2[j];
			}
			k++;
		}
	}
	MPI_Finalize();
	return 0;
}


