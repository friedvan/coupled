
#include "random.h"
//#include "visual.h"
//#include "save.h"
#include <pthread.h>

//node *A, *B;//[N][N], B[N][N];
//stack *sps;


couplednetwork::couplednetwork()
{
	A = new network(N);
	B = new network(N);
}

couplednetwork::~couplednetwork()
{
	if (A) {
		delete A;
	}
		
	if (B) {
		delete B;
	}  			
}

network::network(int network_size)
{
	G = new vertex[network_size];
	if (G == NULL) {
		printf("memory error");
		exit(0);
	}
	clear();
	sps = stack_init();
}

network::~network()
{
	if (G)
		delete[] G;
	stack_release(sps);

}

vertex::vertex()
{
	interconnect = stack_init();
}

vertex::~vertex()
{
	stack_release(interconnect);
}

void vertex::init(int id)
{
	interdependent = id;
	//G[k].inter = -1;			//not defaultly set
	alive = true;
	cluster = 0;

}
//clear nodes' status, and set it back to initial state
void network::clear()
{
	for (int k = 0; k < N; k++) {
		G[k].init(k);
	}
}

// lattice network 
void network::lattice()
{
	int m = sqrt((double)N);
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < m; j++) {
			stack_push(G[i*m + j].interconnect, ((i - 1 + m) % m) * m + j);//up
			stack_push(G[i*m + j].interconnect, ((i + 1 + m) % m) * m + j);//down
			stack_push(G[i*m + j].interconnect, i * m + (j - 1 + m) % m);//left
			stack_push(G[i*m + j].interconnect, i * m + (j + 1 + m) % m);//right
		}
	}
}

//ER graph with N nodes, connecting nodes with probablity p and with distance L
//replace probability with portion, 
//speed up network generating process with time costing of O(N^2)*T(rand()) to O(pN^2)*T(rand()), in which p is very small usually.
void network::ER_length1(double p, int L)
{
	for (int i = 0; i < N; i++) {
		stack_clear(G[i].interconnect);
	}

	int m = sqrt((double)N);		//square length
	int n = (int)(p*N*(N - 1) / 2);	//number of edges to be connected

	int r1, r2, x1, y1, minx, miny, maxx, maxy, x2, y2, dx, dy, distance;
	for (int k = 0; k < n; k++)	{
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
			stack_push(G[r1].interconnect, r2);
			stack_push(G[r2].interconnect, r1);
		}
		else k--;

		//	if (k % 1000 == 0)
		//		printf("%d\t%d\n", n, k);
	}
}


void couplednetwork::init_attack(network *G1, network *G2, double p)
{
	for (int j = 0; j < N; j++) {
		double r = (double)rand() / RAND_MAX;
		if (r < p) {
			G1->G[j].alive = false;
			//set interdependent vertex to death
			G2->G[G1->G[j].interdependent].alive = false;
		}
		//img_print(A, true);
		//img_print(B, false);
	}
}

//每个网络一个独立的sps，用户dfs栈
int network::dfs(int pt, int lable)
{
	int size = 1;
	stack_clear(sps);
	stack_push(sps, pt);

	while (stack_len(sps) != 0) {
		vertex *top = G + *(sps->tail);
		if (!top->alive)
			continue;
		top->cluster = lable;
		bool changed = false;
		int neighbor_len = stack_len(top->interconnect);
		for (int i = 0; i < neighbor_len; i++) {
			vertex *neighbor = G + top->interconnect->head[i];
			if (neighbor->alive && neighbor->cluster == 0) {//alive and not visited
				neighbor->cluster = lable;
				size++;
				stack_push(sps, top->interconnect->head[i]);
				changed = true;
			}
		}
		if (!changed)
			stack_pop(sps);
	}
	return size;
}
void couplednetwork::gaint_component(network *G1, network *G2)
{
	int label = 1;
	int maxsize = 0, maxcluster = -1, size = 0;//, monosize=0, dimersize=0;
	int secondsize = 0, secondcluster = -1;
	
	//for (int j = 0; j<N && size <= N / 2 + 1; j++) {	//if size of cc larger than half of the nodes, it must be the largest
	for (int j = 0; j<N && size <= N; j++) {//find second largest 
		if (G1->G[j].alive && G1->G[j].cluster == 0) {//alive and not visited.
			size = G1->dfs(j, label);
			if (size > maxsize) {				
				secondsize = maxsize;//first copy max to second, then update max
				secondcluster = maxcluster;
				maxcluster = label;
				maxsize = size;
			}
			else if(size > secondsize) {
				secondsize = size;
				secondcluster = label;
			}
			label++;
		}
	}

	G1->maxcc_size = maxsize;
	G1->secondcc_size = secondsize;
	


	//set non gaint component vertex to death
	for (int j = 0; j < N; j++) {
		if (G1->G[j].alive) {
			if (G1->G[j].cluster != maxcluster && G1->G[j].cluster != secondcluster){ //if cluster label is not max and second
				G1->G[j].alive = false;
				if (G1->G[j].interdependent != -1)	 //-1 means has no interdependent nodes
					G2->G[G1->G[j].interdependent].alive = false;
				//img_print(A, true);
				//img_print(B, false);
			}
			G1->G[j].cluster = 0;
		}
	}
}

double average_degree(network *net)
{
	int degree = 0;
	for (int i = 0; i < N; i++)
	{
		degree += stack_len((net->G + i)->interconnect);
	}
	return degree / (float)N;
}

void couplednetwork::random_couple()
{
	//get a random permulation array of length N
	int SHUFFLE_TIMES = 5;
	int arr[N];
	for (int i = 0; i < N; i++) {
		arr[i] = i;
	}
	for (int k = 0; k < SHUFFLE_TIMES; k++) {
		for (int i = 0; i < N; i++) {
			//get a random number smaller than N
			int r = rand() % N;
			//swap element in position i and r
			int temp = arr[i];
			arr[i] = arr[r];
			arr[r] = temp;
		}
	}
	//random couple(interdependent) nodes in network A and B
	for (int i = 0; i < N; i++) {
		A->G[i].interdependent = arr[i];
		B->G[arr[i]].interdependent = i;
	}
}
void* run(void *param)
{
	//prase args
	thread_args * arg = (thread_args*)param;
	int minl = arg->minl;
	int maxl = arg->maxl;

	//set output file
	char filename[100];
	sprintf(filename, "data/result_%d_%d_%d_%d.dat", N, AVERAGE_DEGREE, minl, maxl);
	FILE *fp = fopen(filename, "w");
	if (!fp)
		printf("Bad file or directory!");

	couplednetwork cn;

	cn.A->lattice(); 
	for (double c = minl; c < maxl; c++) {
		//cn.A->ER_length1((double)AVERAGE_DEGREE / N, c);
		//network B
		cn.B->ER_length1((double)AVERAGE_DEGREE / N, c);
		
		//double avg_d = average_degree(cn.B);
		for (double p = 0.01; p <= 1.0; p += 0.01) {
			for (int k = 0; k < NSAMPLE; k++) {
				//get random targets
				//cn.random_couple();
				cn.init_attack(cn.A, cn.B, 1-p);

				int pre_cluster_size = 0, cluster_size = 0;
				int iter = 0;

				while (1) {
					iter++;
					//largest cc in A live
					cn.gaint_component(cn.A, cn.B);
					cluster_size = cn.A->maxcc_size;
					if (cluster_size == pre_cluster_size)
						break;
					pre_cluster_size = cluster_size;
					//largest cc in B live
					cn.gaint_component(cn.B, cn.A);			 				
				}
				
				//printf("%d\t%f\t%f\t%d\t%d\t%d\t%d\n", k, c, p, cluster_size, s.monosize, s.dimersize, iter);
				fprintf(fp, "%d\t%f\t%f\t%d\t%d\t%d\n", k, c, p, cn.A->maxcc_size, cn.A->secondcc_size, iter);
				//fprintf(fp, "%d\t%f\t%f\t%d\t%f\t%d\n", k, c, p, cn.A->maxcc_size, avg_d, iter);

				cn.A->clear();
				cn.B->clear();
			}
			printf("%.3f\t%.3f\n", c, p);
		}
	}

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

	//couplednetwork cn;
	pthread_t tid[NTHREAD];
	for (int i = 0; i < NTHREAD; i++) {
		//memset(&tid[i], 0, sizeof(tid)); 
		int err = pthread_create(&tid[i], NULL, run, (void *)(args + i));
		if (err != 0)	{
			printf("create thread error: %s\n", strerror(err));
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