
#include "random.h"
//#include "visual.h"
//#include "save.h"
#include <pthread.h>

node *A, *B;//[N][N], B[N][N];
//stack *sps;

//clear nodes' status, and set it back to initial state
void clear(node *G)
{
	for (int k = 0; k < N; k++) {
		G[k].inter = k;
		//G[k].inter = -1;			//not defaultly set
		G[k].alive = true;
		G[k].cluster = 0;
		//G[k].type = MONOMER;
		//stack_clear(G[k].base);
	}
}

//allocate memorys for nodes in graph, no neighbors
node* init()
{
	node *G = (node*)malloc(sizeof(node)* N);
	if (G == NULL) {
		printf("memory error");
		exit(0);
	}
	for (int k = 0; k < N; k++)
		G[k].base = stack_init();
	clear(G);
	return G;
}

void release(node *G)
{
	for (int k = 0; k < N; k++)
		stack_release(G[k].base);
	free(G);
}	


// lattice network 
void lattice(node *G)
{
	int m = sqrt((double)N);
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < m; j++) {
			stack_push(G[i*m+j].base, ((i - 1 + m) % m) * m + j);//up
			stack_push(G[i*m+j].base, ((i + 1 + m) % m) * m + j);//down
			stack_push(G[i*m+j].base, i * m + (j - 1 + m) % m);//left
			stack_push(G[i*m+j].base, i * m + (j + 1 + m) % m);//right
		}									 		
	}
}

//ER graph with N nodes and connection probablity p
//void ER(node *G, double p)
//{
//	double r = 0.0;
//	for (int i = 0; i < N; i++) {
//		//time_t tr = 0, tp = 0;
//		for (int j = i + 1; j < N; j++) {
//			//time_t t3, t2, t1 = clock();
//			r = (double)rand() / RAND_MAX;
//			//t2 = clock();
//			//tr += t2 - t1;
//			//rand();
//			//printf("%f\n", r);
//			if (r < p) {
//				stack_push(G[i].base, j);
//				stack_push(G[j].base, i);
//			}
///*			t3 = clock();
//			tp += t3 - t2;   */			
//		}
//		if (i%1000==0)
//			printf("%d\n", i);
//	}
//}
//FILE *ffp = fopen("data/avg_dis.dat", "w");
//ER graph with N nodes, connecting nodes with probablity p and with distance l
//implemented with true random method
void ER_length(node *G, double p, int l)
{
	for (int i = 0; i < N; i++) {
		stack_clear(G[i].base);
	}
	double r = 0.0;
	//periodic boundry
	int m = sqrt((double)N);
	for (int i = 0; i < N; i++) {
		for (int j = i+1; j < N; j++) {	
			int dx=0, dy=0, distance = 0;
			dx = abs(i / m - j / m);
			dy = abs(i % m - j % m);
			distance = MIN(dx, m-dx) + MIN(dy, m-dy);
			
			if (distance > l)
				continue;
			r = (double)rand() / RAND_MAX;
			if (r < p) {
#ifdef _DEBUG
				printf("yes\t%d\t%d\n", i, j);
#endif
				//printf("%d\t%d\t%f\t%d\n", i, j, p, distance);
				//fprintf(ffp, "%f\t%d\n", p, distance);
				stack_push(G[i].base, j);
				stack_push(G[j].base, i);
			} 
		}
		//printf("%d\n", i);
		if (i % 1000 == 0)
			printf("%d\n", i);
	}
}
//ER graph with N nodes, connecting nodes with probablity p and with distance L
//replace probability with portion, 
//speed up network generating process with time costing of O(N^2)*T(rand()) to O(pN^2)*T(rand()), in which p is very small usually.
void ER_length1(node *G, double p, int L)
{
	for (int i = 0; i < N; i++) {
		stack_clear(G[i].base);
	}
	
	int m = sqrt((double)N);		//square length
	int n = (int)(p*N*(N - 1) / 2);	//number of edges to be connected

	int r1, r2, x1, y1, minx, miny, maxx, maxy, x2, y2, dx, dy, distance;
	for (int k = 0; k < n;k++)	{		
		r1 = rand() % N;
		x1 = r1 / m;
		y1 = r1 % m;
		minx = (x1 - L) > 0 ? (x1 - L) : 0;
		maxx = (x1 + L) > m ? m : (x1 + L);
		miny = (y1 - L) > 0 ? (y1 - L) : 0;
		maxy = (y1 + L) > m ? m : (y1 + L);		
		x2 = minx + rand() % (maxx - minx);
		y2 = miny + rand() % (maxy - miny);
		r2 = x2 * m + y2;

		dx = abs(r1 / m - r2 / m);
		dy = abs(r1 % m - r2 % m);
		//distance = MIN(dx, m - dx) + MIN(dy, m - dy);	   //periodic boundry, 
		distance = dx + dy;

		if (distance <= L) {			//connection distance shorter than L	
			//connect r2 as r1's neighbr.
			/*
			NOTICE: if r2 is already r1's neighbor, this will add r2 again, 
			but it does not matter because the probability of geting two identical edge
			is about n^2/N^4=p^2, usually p is very small.  
			*/
			stack_push(G[r1].base, r2);	
			stack_push(G[r2].base, r1);
		} 
		else k--;			

	//	if (k % 1000 == 0)
	//		printf("%d\t%d\n", n, k);
	}
}

void ER_length2(node *G, double p, double L)
{
	for (int i = 0; i < N; i++) {
		stack_clear(G[i].base);
	}

	int m = sqrt((double)N);		//square length
	int n = (int)(p*N*(N - 1) / 2);	//number of edges to be connected
	for (int k = 0; k < n; k++)	{
		int r1 = rand() % N;
		int r2 = rand() % N;
		int dx = 0, dy = 0, distance = 0;
		dx = abs(r1 / m - r2 / m);
		dy = abs(r1 % m - r2 % m);
		distance = MIN(dx, m - dx) + MIN(dy, m - dy);	   //periodic boundry, 
		//distance = dx + dy;

		if (distance <= L) {			//connection distance shorter than L	
			//connect r2 as r1's neighbr.
			/*
			NOTICE: if r2 is already r1's neighbor, this will add r2 again,
			but it does not matter because the probability of geting two identical edge
			is about n^2/N^4=p^2, usually p is very small.
			*/
			stack_push(G[r1].base, r2);
			stack_push(G[r2].base, r1);
		}
		else k--;

		//	if (k % 1000 == 0)
		//		printf("%d\t%d\n", n, k);
	}
}

void release(int *list)
{
	if (list != NULL) {
		free(list);
		list = NULL;
	}
}

gcsize get_size(node *G)
{
	gcsize s;
	memset(&s, 0, sizeof(s));

	for (int j = 0; j < N; j++) {
		if (G[j].alive) {
			s.maxsize++;
			//if (G[j].type == MONOMER)
			//	s.monosize++;
			//else
			//	s.dimersize++;
		}
	}
	return s;
}

void init_attack(double p)
{
	srand(rand());

	for (int j = 0; j < N; j++) {
		double r = (double)rand() / RAND_MAX;
		if (r < p) {
			A[j].alive = false;
			if (A[j].inter != -1)	   //if has interdependent node
				B[A[j].inter].alive = false;
		}
		//img_print(A, true);
		//img_print(B, false);
	}
}

//每个线程一个独立的sps stack 用来计算dfs，线程内部通用，线程间独立
int dfs(node *G, int pt, int lable, stack *sps)
{
	int size = 1;
	stack_clear(sps);
	stack_push(sps, pt);

	while (stack_len(sps) != 0) {
		node *top = G + *(sps->tail);
		if (!top->alive)
			continue;
		top->cluster = lable;
		bool changed = false;
		for (int i=0; i<stack_len(top->base); i++) {
			node *neighbor = G + top->base->head[i];
			if(neighbor->alive && neighbor->cluster == 0) {//alive and not visited
				neighbor->cluster = lable;
				size++;
				stack_push(sps, top->base->head[i]);
				changed = true;
			}
		}
		if (!changed)
			stack_pop(sps);
	}
	return size;
}
void gaint_component(node *G1, node *G2, stack *sps)
{
	int lable = 1;
	int maxsize = 0, maxcluster = -1, size = 0;//, monosize=0, dimersize=0;

	for (int j = 0; j<N && size <= N / 2 + 1; j++) {
		if (G1[j].alive && G1[j].cluster == 0) {//alive and not visited.
			size = dfs(G1, j, lable, sps);
			if (size > maxsize) {
				maxcluster = lable;
				maxsize = size;
			}
			lable++;
		}
	}


	for (int j = 0; j < N; j++) {
		if (G1[j].alive) {
			if (G1[j].cluster != maxcluster){
				G1[j].alive = false;
				if (G1[j].inter != -1)
					G2[G1[j].inter].alive = false;
				//img_print(A, true);
				//img_print(B, false);
			}
			G1[j].cluster = 0;
		}
	}
}

void avg_degree(node *G)
{
	int sum = 0;
	for (int i = 0; i < N; i++) {
		sum += stack_len(G[i].base);
	}
	printf("%f\n", (double)sum / N);
}


void* run(void *param)
{
	thread_args * arg = (thread_args*)param;
	int minl = arg->minl;
	int maxl = arg->maxl;


	char filename[100];	
//#ifdef _DEBUG
//	sprintf(filename, "data/result_%d_%d_debug.dat", minl, maxl);
//	FILE *fp = fopen(filename, "w");
//	
//#else					  	
//	sprintf(filename, "data/result_%d_%d.dat", minl, maxl);
//	FILE *fp = fopen(filename, "w");
//#endif
	sprintf(filename, "data/result_%d_%d.dat", minl, maxl);
	printf("data/result_%d_%d.dat", minl, maxl);
	FILE *fp = fopen(filename, "w");
	if (!fp)
		printf("Bad file or directory!");

	A = init();//allocate memorys for nodes in graph
	B = init();

	stack *sps = stack_init();

	//ER network with probability p
	//for speed consideration, I moved it out of loop.
	//random network should regenerate every relazation   

	//avg_degree(A);
	lattice(A);
	
	//for (double c = 0.5; c <= 10.0; c+= 0.5) {
	for (double c = minl; c <= maxl; c++) {
		//ER_length1(A, 5.0 / N, c);
		ER_length1(B, (double)AVERAGE_DEGREE / N, c);
		//A = B;
		//save_network(B);
		//break;
		
		for (double p = 0.01; p <= 1.0; p += 0.01) {
			for (int k = 0; k < NSAMPLE; k++) {
				gcsize s;
				init_attack(1 - p);
				int pre_cluster_size = 0, cluster_size = 0;
				int iter = 0;
				//gcsize s1, s2;
				//memset(&s1, 0, sizeof(s1));
				//memset(&s2, 0, sizeof(s2));
				while (1) {
					iter++;
					gaint_component(A, B, sps);
					s = get_size(A);
					//img_print(A, true);
					cluster_size = s.maxsize;
					if (cluster_size == pre_cluster_size)
						break;
					pre_cluster_size = cluster_size;
					gaint_component(B, A, sps);
					//img_print(B, false);
				}
				s = get_size(A);
				//printf("%d\t%f\t%f\t%d\t%d\t%d\t%d\n", k, c, p, cluster_size, s.monosize, s.dimersize, iter);
				fprintf(fp, "%d\t%f\t%f\t%d\t%d\t%d\t%d\n", k, c, p, cluster_size, s.monosize, s.dimersize, iter);

				clear(A);
				clear(B);
			}
			printf("%.2f\t%.2f\n", c, p);
		}
	}


	release(A);
	release(B);

	stack_release(sps);
	fclose(fp);
	return NULL;
}

int main()
{
	srand(time(NULL));
	//single thread
	//thread_args args;
	//args.minl = 1;
	//args.maxl = 2;
	//run(&args);

	//multithread
	thread_args *args = (thread_args*)malloc(sizeof(thread_args)*NTHREAD);
	args[0].minl = MINL;
	args[0].maxl = MINL + LOAD_PER_THREAD;
	for (int i = 1; i < NTHREAD; i++) {
		args[i].minl = args[i - 1].maxl;
		args[i].maxl = args[i].minl + LOAD_PER_THREAD;
	}

	pthread_t tid[NTHREAD];
	for (int i = 0; i < NTHREAD; i++) {
		memset(&tid[i], 0, sizeof(tid)); 
		int err = pthread_create(&tid[i], NULL, run, (void *)(args+i));
		if (err!=0)	{
			printf("create thread error: %s/n", strerror(err));
			exit(0);
		}
		//printf("create thread: %d\n", tid[i]);
	}

	for (int i = 0; i < NTHREAD; i++) {	 
		pthread_join(tid[i], NULL);
		//printf("waiting %d", tid[i]);	
	}


	return 0;
}