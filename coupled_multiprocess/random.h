
#define _CRT_SECURE_NO_WARNINGS

#ifdef _WIN32
#pragma comment(lib, "pthreadVC2.lib")
#include <Windows.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <algorithm>
//#include <iostream>
#include <vector>
//#include <hash_set>
#include <unordered_map>
#include <thread>
using namespace std;

//#include "visual.h"

static const int LATTICE_SIZE = 100;
static const int N = LATTICE_SIZE*LATTICE_SIZE;
#define NSAMPLE 100 

#define ER_AVERAGE_DEGREE 5
/*#define SF_GAMMA 3.0
#define SF_GAMMA_NORM_FACTOR 1.202

3.0	1.202
2.9	1.223
2.8 1.247
2.7 1.274
2.6 1.305
2.5 1.341
*/
#define NTHREAD 2
#define MINL 2
#define LOAD_PER_THREAD 3

#define MIN(a,b) ((a)<(b)?(a):(b))

//stack stucture
/*********************************************************************/
class stack
{
public:
	int *head, *tail;
	int size;
	int empty;
	int	STACK_INCREASEMENT;
	int STACK_INITIAL;
	
	void stack_init(int stack_initial_size=20, int stack_incrament_size=20)
	{
		STACK_INCREASEMENT = stack_incrament_size;
		STACK_INITIAL = stack_initial_size;
		empty = STACK_INITIAL;
		size = STACK_INITIAL;
		head = (int*)malloc(sizeof(int)*STACK_INITIAL);
		if (head == NULL) {
			printf("memory error");
			exit(0);
		}
		tail = head - 1;//+s.len;
	}

	void stack_push(int pt)
	{
		if (empty < 1) {
			head = (int *)realloc(head, sizeof(int)* (size + STACK_INCREASEMENT));
			if (head == NULL) {
				printf("memory error");
				exit(0);
			}
			//memcpy(tmp, ps->head, ps->size);
			//ps->head = tmp;
			tail = head + size - 1;
			size += STACK_INCREASEMENT;
			empty = STACK_INCREASEMENT;

		}
		tail++;
		*(tail) = pt;
		//ps->tail->y = pt.y;
		empty--;
	}

	int stack_pop()
	{
		int pt = *(tail);
		tail--;
		empty++;
		return pt;
	}

	void stack_release()
	{
		if (head != NULL) {
			free(head);
			head = NULL;
		}
	}

	int stack_len()
	{
		return tail - head + 1;
	}

	void stack_clear()
	{
		empty = size;
		tail = head - 1;
	}
};

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
	//unordered_map<int, bool> hm;
	//hash_set<int> hs;
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

	void init();

	int dfs(int pt, int lable);

	void lattice();
	void ER_length1(double p, int L);	
	vector<int> node_sequence_scale_free(double gamma);
	void configuration_model(double gamma, int L);
	vector<int> node_degree_scale_free(double gamma);
	void random_geometric_graph(double p, int L);
	//void scale_free_length(int L);
	//void store_nodes_within_distance(stack *pool, bool exist[], int node, int L);
};

class couplednetwork
{
public :
	couplednetwork();
	~couplednetwork();
	network *A;
	network *B;

	//stack *scps;

	void random_couple(void);
	void gaint_component(network *G1, network *G2);
	void init_attack(network *G1, network *G2, double p);
	bool attack(network *G1, network *G2);
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




