
#define _CRT_SECURE_NO_WARNINGS

#ifdef _WIN32
#pragma comment(lib, "pthreadVC2.lib")
#endif

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
//#include <iostream>
#include <vector>
using namespace std;

//#include "visual.h"

static const int LATTICE_SIZE = 1000;
static const int N = LATTICE_SIZE*LATTICE_SIZE;
#define NSAMPLE 10	 

#define AVERAGE_DEGREE 5
#define NTHREAD 20
#define MINL 2
#define LOAD_PER_THREAD 50 

#define MIN(a,b) ((a)<(b)?(a):(b))

//stack stucture
/*********************************************************************/
#define	STACK_INCREASEMENT 50
#define STACK_INITIAL 100
typedef struct stack
{
	int *head, *tail;
	int size;
	int empty;
}stack;

stack* stack_init()
{
	stack *ps = (stack *)malloc(sizeof(stack));
	if (ps == NULL) {
		printf("memory error");
		exit(0);
	}
	ps->empty = STACK_INITIAL;
	ps->size=STACK_INITIAL;
	ps->head=(int*)malloc(sizeof(int)*STACK_INITIAL);
	ps->tail=ps->head-1;//+s.len;
	return ps;
}

void stack_push(stack *ps, int pt)
{
	if (ps->empty<1) {
		ps->head = (int *)realloc(ps->head, sizeof(int) * (ps->size+STACK_INCREASEMENT));
		if (ps->head == NULL) {
			printf("memory error");
			exit(0);
		}
		//memcpy(tmp, ps->head, ps->size);
		//ps->head = tmp;
		ps->tail = ps->head + ps->size - 1;
		ps->size += STACK_INCREASEMENT;
		ps->empty = STACK_INCREASEMENT;

	}
	ps->tail++;
	*(ps->tail) = pt;
	//ps->tail->y = pt.y;
	ps->empty--;
}

int stack_pop(stack *ps)
{
	int pt = *(ps->tail);
	ps->tail--;
	ps->empty++;
	return pt;
}

void stack_release(stack *ps)
{
	if (ps->head != NULL) {
		free(ps->head);
		ps->head = NULL;
	}
	free(ps);
}

int stack_len(stack *ps)
{
	return ps->tail - ps->head+1;
}

void stack_clear(stack *ps)
{
	ps->empty = ps->size;
	ps->tail = ps->head - 1;   
}

//nodes, graph
/******************************************************************/

class vertex
{
public:
	vertex();
	~vertex();
	//point id;
	int interdependent;		 //interdepent neighbor
	stack *interconnect;	 //interconnect neighbors
	bool alive;
	int cluster;
	
	void init(int id);
};

class network
{
public:
	vertex *G;
	stack* sps;

	int maxcc_size;
	int secondcc_size; 

	network(int network_size);
	~network();

	void clear();

	int dfs(int pt, int label);

	void lattice();
	void ER_length1(double p, int L);	
};

class couplednetwork
{
public :
	couplednetwork();
	~couplednetwork();
	network *A;
	network *B;

	void gaint_component(network *G1, network *G2);
	void init_attack(network *G1, network *G2, double p);
};


//save network to file
/******************************************************************/
//void save_network(node *G)
//{
//	FILE *fp = fopen("B.csv", "w");
//	fprintf(fp, "node,pathid,type,x,y\n");
//	int pathid = 0;
//	for (int k = 0; k < N; k++) {
//		int fx = k / LATTICE_SIZE;
//		int fy = k % LATTICE_SIZE; 		
//		for (int i = 0; i<stack_len(G[k].base); i++){  
//			pathid++;
//			int neighbor = G[k].base->head[i];
//			int nx = neighbor / LATTICE_SIZE;
//			int ny = neighbor % LATTICE_SIZE;
//			fprintf(fp, "%d,%d,%d,%d,%d\n", k, pathid, 1, fx, fy);
//			fprintf(fp, "%d,%d,%d,%d,%d\n", neighbor, pathid, 2, nx, ny);
//		}
//	}
//	fclose(fp);
//}
//multithread function arguments structure
/******************************************************************/
typedef struct thread_args
{
	int minl;
	int maxl;
}thread_args;




