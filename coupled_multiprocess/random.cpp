
#include "random.h"
//#include "visual.h"

node *A, *B;//[N][N], B[N][N];
stack *sps;

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

//ER graph with N nodes and connection probablity p
void ER_length(node *G, double p, int l)
{
	for (int i = 0; i < N; i++) {
		stack_clear(G[i].base);
	}
	double r = 0.0;
	int m = sqrt((double)N);
	for (int i = 0; i < N; i++) {
		//clock_t dt1 = 0, dt2 = 0, dt3 = 0, dt4 = 0;
		for (int j = i+1; j < N; j++) {	
			//clock_t t1, t2, t3, t4;
			//t1 = clock();
			int dx=0, dy=0, dd = 0, distance = 0;
			dd = abs(i - j);
			dx = dd / m;
			dy = dd % m;
			distance = MIN(dx, m-dx) + MIN(dy, m-dy);
			if (distance > l)
				continue;
			//t2 = clock();
			r = (double)rand() / RAND_MAX;
			//t3 = clock();
			if (r < p) {
				stack_push(G[i].base, j);
				stack_push(G[j].base, i);
			} 
			//t4 = clock();
			//dt1 += t2 - t1;
			//dt2 += t3 - t2;
			//dt3 += t4 - t3;
			//dt4 += t4 - t1; 
		}
		//printf("%ld\t%ld\t%ld\t%ld\n", dt1, dt2, dt3, dt4);
		if (i % 1000 == 0)
			printf("%d\n", i);
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


int dfs(node *G, int pt, int lable)
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
void gaint_component(node *G1, node *G2)
{
	int lable = 1;
	int maxsize = 0, maxcluster = -1, size = 0;//, monosize=0, dimersize=0;

	for (int j = 0; j<N && size <= N / 2 + 1; j++) {
		if (G1[j].alive && G1[j].cluster == 0) {//alive and not visited.
			size = dfs(G1, j, lable);
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


int main()
{
	srand(time(NULL));
	FILE *fp = fopen("data/result.dat", "w");
	if (!fp)
		printf("Bad file or directory!");

	A = init();//allocate memorys for nodes in graph
	B = init();

	sps = stack_init();

	//ER network with probability p
	//!!!!!!!!for speed consideration, I moved it out of loop.
	//However, random network should regenerate every relazation   

	//avg_degree(A);
	lattice(A);
	for (int c = 130; c < 1000; c+=10) {
		//ER_length(A, 5.0 / N, c);
		ER_length(B, 5.0 / N, c);
		for (double p=0.6; p<=1.0; p+=0.01) {
			for (int k=0; k<NSAMPLE; k++) {
				gcsize s;	
				init_attack(1-p);
				int pre_cluster_size = 0, cluster_size = 0;				
				int iter = 0;
				//gcsize s1, s2;
				//memset(&s1, 0, sizeof(s1));
				//memset(&s2, 0, sizeof(s2));
				while (1) {
					iter++;
					gaint_component(A, B);
					s=get_size(A);
					//img_print(A, true);
					cluster_size = s.maxsize;
					if (cluster_size == pre_cluster_size)
						break;
					pre_cluster_size = cluster_size;
					gaint_component(B, A);
					//img_print(B, false);
				}
				s=get_size(A);
				printf("%d\t%d\t%f\t%d\t%d\t%d\t%d\n", k, c, p, cluster_size, s.monosize, s.dimersize, iter);
				fprintf(fp, "%d\t%d\t%f\t%d\t%d\t%d\t%d\n", k, c, p, cluster_size, s.monosize, s.dimersize, iter);

				clear(A);
				clear(B);
			}
			printf("%d\t%.2f\n", c, p);
		} 
	}

	release(A);
	release(B);

	stack_release(sps);
	fclose(fp);
}