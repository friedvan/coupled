
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
	//scps->stack_init();
}

couplednetwork::~couplednetwork()
{
	if (A) {
		delete A;
		A = NULL;
	}
		
	if (B) {
		delete B;
		B = NULL;
	}
	//if (scps) {
	//	scps->stack_release();
	//	delete scps;
	//	scps = NULL;
	//}
}

network::network(int network_size)
{
	G = new vertex[network_size];
	if (G == NULL) {
		printf("memory error");
		exit(0);
	}
	init();
	sps = new stack();
	sps->stack_init(100, 100);
}

network::~network()
{
	if (G)
		delete[] G;
	if (sps) {
		sps->stack_release();
		delete sps;	
	}
}

vertex::vertex()
{
	//interconnect = stack_init();
	interconnect = new stack();	
	interconnect->stack_init();
}

vertex::~vertex()
{
	//stack_release(interconnect);
	if (interconnect) {
		interconnect->stack_release();
		delete interconnect;
		interconnect = NULL;
	}
		
}

void vertex::init(int id)
{
	interdependent = id;
	//G[k].inter = -1;			//not defaultly set
	alive = true;
	cluster = 0;
	//for (auto it : hm)	{
	//	hm[it.first] = true;
	//}

}
//clear nodes' status, and set it back to initial state
void network::init()
{
	for (int k = 0; k < N; k++) {
		G[k].init(k);
	}
}

// lattice network 
void network::lattice()
{
	//for (int i = 0; i < N; i++) {
	//	//stack_clear(G[i].interconnect);
	//	G[i].interconnect->stack_clear();
	//}
	int m = sqrt((double)N);
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < m; j++) {
			//stack_push(G[i*m + j].interconnect, ((i - 1 + m) % m) * m + j);//up
			//stack_push(G[i*m + j].interconnect, ((i + 1 + m) % m) * m + j);//down
			//stack_push(G[i*m + j].interconnect, i * m + (j - 1 + m) % m);//left
			//stack_push(G[i*m + j].interconnect, i * m + (j + 1 + m) % m);//right

			G[i*m + j].interconnect->stack_push(((i - 1 + m) % m) * m + j);//up
			G[i*m + j].interconnect->stack_push(((i + 1 + m) % m) * m + j);//down
			G[i*m + j].interconnect->stack_push(i * m + (j - 1 + m) % m);//left
			G[i*m + j].interconnect->stack_push(i * m + (j + 1 + m) % m);//right 
			//G[i*m + j].hs.insert(((i - 1 + m) % m) * m + j);
			//G[i*m + j].hs.insert(((i + 1 + m) % m) * m + j);
			//G[i*m + j].hs.insert(i * m + (j - 1 + m) % m);
			//G[i*m + j].hs.insert(i * m + (j + 1 + m) % m);
			//G[i*m + j].hm[((i - 1 + m) % m) * m + j] = true;
			//G[i*m + j].hm[((i + 1 + m) % m) * m + j] = true;
			//G[i*m + j].hm[i * m + (j - 1 + m) % m] = true;
			//G[i*m + j].hm[i * m + (j + 1 + m) % m] = true;
		}
	}
}

//void network::random_geometric_graph(double p, int L) {
//	for (int i = 0; i < N; i++) {
//		//stack_clear(G[i].interconnect);
//		G[i].interconnect->stack_clear();
//	}
//	bool alive[N] = { false };
//	for (int i = 0; i<N; i++) {
//		if (double(rand()) / RAND_MAX < (double)L/LOAD_PER_THREAD)  //randomly choose p*N nodes
//			alive[i] = true;
//	}
//	
//	int m = LATTICE_SIZE;
//	int r1, r2, x1, y1, minx, miny, maxx, maxy, dx, dy, distance;
//
//	for (int k = 0; k < N; k++)	{
//		if (!alive[k]) continue;
//		r1 = k;
//		x1 = r1 / m;
//		y1 = r1 % m;
//		minx = x1 - L;
//		maxx = x1 + L;
//		miny = y1 - L;
//		maxy = y1 + L;
//
//		int area = 4 * L*L;
//		for (int x2 = minx; x2 < maxx; x2++){
//			int x = (x2 + m) % m;
//			for (int y2 = miny; y2 < maxy; y2++) {
//				int y = (y2 + m) % m;
//				r2 = x * m + y;
//				if (!alive[r2] || r1 >= r2)
//					continue;	//not alive or already connected
//				dx = abs(r1 / m - r2 / m);
//				dy = abs(r1 % m - r2 % m);
//				distance = MIN(dx, m - dx) + MIN(dy, m - dy);	   //periodic boundry, 
//				if (distance <= 2) {			//connection distance shorter than L
//					//if ((double)rand() / RAND_MAX>10.0 / area) continue;
//					//connect r2 as r1's neighbr.
//					/*
//					NOTICE: if r2 is already r1's neighbor, this will add r2 again,
//					but it does not matter because the probability of geting two identical edge
//					is about n^2/N^4=p^2, usually p is very small.
//					*/
//					G[r1].interconnect->stack_push(r2);
//					G[r2].interconnect->stack_push(r1);
//				}
//			}
//		}
//	}
//}
//void network::random_geometric_graph(double p, int L) {
//	for (int i = 0; i < N; i++) {
//		stack_clear(G[i].interconnect);
//	}
//	bool alive[N] = {false};
//	for (int i = 0; i<N; i++) {
//		if (double(rand()) / RAND_MAX < 1.01)  //randomly choose p*N nodes
//			alive[i]=true;
//	}
//	
//	int m = LATTICE_SIZE;
//	int r1, r2, x1, y1, minx, miny, maxx, maxy, dx, dy, distance;
//
//	for (int k = 0; k < N; k++)	{
//		if (!alive[k]) continue;
//		r1 = k;
//		x1 = r1 / m;
//		y1 = r1 % m;
//		minx = x1 - L;
//		maxx = x1 + L;
//		miny = y1 - L;
//		maxy = y1 + L;
//
//		int area = 4*L*L; 		
//		for (int x2 = minx; x2 < maxx; x2++){
//			int x = (x2 + m) % m;
//			for (int y2 = miny; y2 < maxy; y2++) {
//				int y = (y2 + m) % m;
//				r2 = x * m + y;
//				if (!alive[r2] || r1 >= r2)
//					continue;	//not alive or already connected
//				dx = abs(r1 / m - r2 / m);
//				dy = abs(r1 % m - r2 % m);
//				distance = MIN(dx, m - dx) + MIN(dy, m - dy);	   //periodic boundry, 
//				if (distance <= L) {			//connection distance shorter than L
//					if ((double)rand() / RAND_MAX>10.0/area) continue;
//					//connect r2 as r1's neighbr.
//					/*
//					NOTICE: if r2 is already r1's neighbor, this will add r2 again,
//					but it does not matter because the probability of geting two identical edge
//					is about n^2/N^4=p^2, usually p is very small.
//					*/
//					stack_push(G[r1].interconnect, r2);
//					stack_push(G[r2].interconnect, r1);
//				}
//			}
//		}
//	}
//}

//ER graph with N nodes, connecting nodes with probablity p and with distance L
//replace probability with portion, 
//speed up network generating process with time costing of O(N^2)*T(rand()) to O(pN^2)*T(rand()), in which p is very small usually.
void network::ER_length1(double p, int L)
{
	for (int i = 0; i < N; i++) {
		G[i].interconnect->stack_clear();
		//G[i].hm.clear();
	}

	int m = LATTICE_SIZE;		//square length
	int n = (int)(p*N*(N - 1) / 2);	//number of edges to be connected

	int r1, r2, x1, y1, minx, miny, maxx, maxy, x2, y2, dx, dy, distance;
	for (int k = 0; k < n; k++)	{
		r1 = rand() % N;
		x1 = r1 / m;
		y1 = r1 % m;
		minx = (x1 - L);// > 0 ? (x1 - L) : 0;
		maxx = (x1 + L);// > m ? m : (x1 + L);
		miny = (y1 - L);// > 0 ? (y1 - L) : 0;
		maxy = (y1 + L);// > m ? m : (y1 + L);
		x2 = minx + rand() % (maxx - minx);
		y2 = miny + rand() % (maxy - miny);
		r2 = ((x2 + m) % m) * m + (y2 + m) % m;

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
			G[r1].interconnect->stack_push(r2);
			//G[r1].hm[r2] = true;
			//G[r2].hm[r1] = true;
						  			
			G[r2].interconnect->stack_push(r1);
		}
		else k--;

		//	if (k % 1000 == 0)
		//		printf("%d\t%d\n", n, k);
	}
}

//void network::store_nodes_within_distance(stack *pool, bool exist[], int node, int L) {
//	int minx, miny, maxx, maxy, rx, ry, dx, dy, distance, pt;
//	rx = node / LATTICE_SIZE;
//	ry = node % LATTICE_SIZE;
//	minx = (rx - L) > 0 ? (rx - L) : 0;
//	maxx = (rx + L) > LATTICE_SIZE ? LATTICE_SIZE : (rx + L);
//	miny = (ry - L) > 0 ? (ry - L) : 0;
//	maxy = (ry + L) > LATTICE_SIZE ? LATTICE_SIZE : (ry + L);
//	for (int x = minx; x <= maxx; x++) {
//		for (int y = miny; y <= maxy; y++) {
//			dx = abs(x - rx);
//			dy = abs(y - ry);
//			distance = MIN(dx, LATTICE_SIZE - dx) + MIN(dy, LATTICE_SIZE - dy);	   //periodic boundry
//			if (distance <= L) {
//				pt = x*LATTICE_SIZE + y;
//				if (!exist[pt])	{
//					stack_push(pool, pt);
//					exist[pt] = true;
//				}
//			}
//		}
//	}
//}

///*the basic idea is:
//	1. randomly pick SF_M0 inital nodes to form a initial Graph, put their neighbors within distance L to Pool
//	2. pick one node from Pool everytime, add it to Graph by:
//		(a) connect it to SF_M nodes in graph with probability of SF_ALPHA
//		(b) or connect it to SF_M nodes in Pool with probability of 1-SF_ALPHA 
//		then delete it from Pool and put its neighbor within distance L to Pool
//	3. loop step2 until there is no more nodes in Pool
//*/
//void network::scale_free_length(int L)
//{
//	//price model
//	int nodes[SF_M*N];	//array to store nodes' degree infomation, add nodes to array of every new edge.
//	int index = 0;			//index for nodes array, point to first 
//
//	stack *pool = stack_init();		//store all possible nodes to be attched to
//	bool pool_exist[N] = { false };		//help pool to check if node already in pool
//	bool nodes_exsit[N] = { false };	//help nodes to check if node already in graph
//	//randomly add m_0 initial nodes to graph
//	int r = rand() % N;
//	for (int i = 0; i < SF_M0; i++) {
//		store_nodes_within_distance(pool, pool_exist, r, L);
//		nodes_exsit[r] = true;
//		for (int j = 0; j<SF_M0; j++) {
//			if (i == j) continue;
//			stack_push(G[r+i].interconnect, r+j);
//		}
//		
//		nodes[index++] = r++;
//	}
//
//	while (stack_len(pool) > 0) {
//		int node1 = stack_pop(pool);
//		if (nodes_exsit[node1])
//			continue;
//		else
//			nodes_exsit[node1] = true;
//		for (int j = 0; j < SF_M; j++) {
//			//double rr1 = (double)rand() / RAND_MAX;
//			//double rr2 = (double)SF_M / (SF_M + SF_ALPHA);
//			//if (rr1 < rr2) { //prefrential attachment
//				int node2 = nodes[rand() % index];		//choose a exsit node with probability proportional to its degree
//				nodes[index++] = node2;	 //from node
//				nodes[index++] = node1;	  //to node
//				stack_push(G[node1].interconnect, node2);
//				stack_push(G[node2].interconnect, node1);
//				store_nodes_within_distance(pool, pool_exist, node1, L);
//			//}
//			//else {				//random attchment
//			//	int node2 = pool->head[rand() % stack_len(pool)];	//randomly choose a node from pool
//			//	if (nodes_exsit[node2]) {
//			//		j--;
//			//		continue;
//			//	}
//			//	nodes[index++] = node2;	 //from node
//			//	nodes[index++] = node1;	  //to node
//			//	stack_push(G[node1].interconnect, node2);
//			//	stack_push(G[node2].interconnect, node1);
//			//	store_nodes_within_distance(pool, pool_exist, node1, L);
//			//	store_nodes_within_distance(pool, pool_exist, node2, L);
//			//}
//		}
//	}
//}
//
//vector<int> network::node_sequence_scale_free(double gamma) {
//	//random nodes
//	int nodes[N];
//	for (int i = 0; i < N; i++) {
//		nodes[i] = i;
//	}
//	random_shuffle(nodes, nodes + N);
//
//	int index = 0;
//	//degree sequence
//	vector<int> node_seq;
//	for (int i = 1; i < N; i++) {
//		/*normalized by a Riemann zeta function to become probability
//		gamma	normalization factor
//		3.0	1.202
//		2.9	1.223
//		2.8 1.247
//		2.7 1.274
//		2.6 1.305
//		2.5 1.341
//		*/
//		int n = (int)(N*(pow((double)i, -gamma) / SF_GAMMA_NORM_FACTOR)); //every degree i map with n nodes
//		if (n == 0) break;
//		for (int j = 0; j < n; j++) {			
//			for (int k = 0; k < i; k++){		//every nodes repeat degree i times
//				node_seq.push_back(nodes[index]);
//			}
//			if (++index >= N)
//				break;		
//		} 		
//	}
//	//random shuffle node sequence
//	random_shuffle(node_seq.begin(), node_seq.end());
//	return node_seq;
//}
//
//
//vector<int> network::node_degree_scale_free(double gamma) {
//	//random nodes
//	//int nodes[N];
//	//for (int i = 0; i < N; i++) {
//	//	nodes[i] = i;
//	//}
//	//random_shuffle(nodes, nodes + N);
//
//	int index = 0;
//	//degree sequence
//	vector<int> node_seq(N);
//	for (int i = 1; i < N; i++) {
//		/*normalized by a Riemann zeta function to become probability
//		gamma	normalization factor
//		3.0	1.202
//		2.9	1.223
//		2.8 1.247
//		2.7 1.274
//		2.6 1.305
//		2.5 1.341
//		*/
//		int n = (int)(N*(pow((double)i, -gamma) / SF_GAMMA_NORM_FACTOR)); //every degree i map with n nodes
//		if (n == 0) break;
//		for (int j = 0; j < n; j++) {
//			node_seq[index] = i;
//			if (++index >= N)
//				break;
//		}
//	}
//	//random shuffle node sequence
//	random_shuffle(node_seq.begin(), node_seq.end());
//	return node_seq;
//}
//
//void network::configuration_model(double gamma, int L)
//{
//	//vector<int> node_seq = node_sequence_scale_free(gamma);
//	vector<int> node_seq = node_degree_scale_free(gamma);
//	int size = 0;
//	for (int i = 0; i < N; i++) {
//		size += node_seq[i];
//	}
//	int x1, y1, minx, miny, maxx, maxy, x2, y2, dx, dy, distance;
//	int m = LATTICE_SIZE;
//	for (int i = 0; i < N; i++) {
//		int node1 = i;
//		if (node_seq[node1] <= 0) continue;	
//		x1 = node1 / m;
//		y1 = node1 % m;
//		minx = (x1 - L) > 0 ? (x1 - L) : 0;
//		maxx = (x1 + L) > m ? m : (x1 + L);
//		miny = (y1 - L) > 0 ? (y1 - L) : 0;
//		maxy = (y1 + L) > m ? m : (y1 + L);
//
//		int failure = 0;
//		while (node_seq[node1] > 0) {
//			x2 = minx + rand() % (maxx - minx);
//			y2 = miny + rand() % (maxy - miny);
//			int node2 = x2 * m + y2;
//			if (node_seq[node2] <= 0 || node1 == node2) continue;
//			int dx = abs(node1 / LATTICE_SIZE - node2 / LATTICE_SIZE);
//			int dy = abs(node1%LATTICE_SIZE - node2%LATTICE_SIZE);
//			int distance = MIN(dx, LATTICE_SIZE - dx) + MIN(dy, LATTICE_SIZE - dy);	   //periodic boundry
//			if (distance <= L) {
//				stack_push(G[node1].interconnect, node2);
//				stack_push(G[node2].interconnect, node1);
//				node_seq[node1]--;
//				node_seq[node2]--;
//			}
//			else {
//				if (failure++ > (maxx-minx)*(maxy-miny)) {
//					failure;
//				}
//			}
//		}	
//	}
//}
//
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

	sps->stack_clear();
	sps->stack_push(pt);


	LARGE_INTEGER  large_interger;
	double dff;
	__int64  c1, c2;
	QueryPerformanceFrequency(&large_interger);
	dff = large_interger.QuadPart;
	QueryPerformanceCounter(&large_interger);
	c1 = large_interger.QuadPart;
	while (sps->stack_len() != 0) {
		vertex *top = G + *(sps->tail);
		if (!top->alive)
			continue;
		top->cluster = lable;
		bool changed = false;
		int neighbor_len = top->interconnect->stack_len();
		for (int i = 0; i < neighbor_len; i++) {
			//for (auto it : top->hm) {
			//	if (!it.second)
			//		continue;
			//	vertex *neighbor = G + it.first;
			vertex *neighbor = G + top->interconnect->head[i];
			if (neighbor->alive && neighbor->cluster == 0) {//alive and not visited
				neighbor->cluster = lable;
				size++;
				sps->stack_push(top->interconnect->head[i]);
				changed = true;
			}
		}
		if (!changed)
			sps->stack_pop();
	}
	QueryPerformanceCounter(&large_interger);
	c2 = large_interger.QuadPart;
	//printf("本机高精度计时器频率%lf\n", dff);
	//printf("第一次计时器值%I64d 第二次计时器值%I64d 计时器差%I64d\n", c1, c2, c2 - c1);
	//printf("%f\n", (c2 - c1) * 1000 / dff);

	return size;
}
//bool couplednetwork::attack(network *G1, network *G2)
//{
//	bool changed = false;
//	for (int i = 0; i < N; i++) {
//		//scps->stack_clear();
//		vertex *top = G1->G + i;
//		if (!top->alive)
//			continue;
//		vertex *inter = G2->G + top->interdependent;
//		//hash_set<int> temp(top->hs);
//		for (auto it : top->hm){
//			if (it.second && inter->hm.end() == inter->hm.find(G1->G[it.first].interdependent)) {
//				changed = true;
//				G1->G[i].hm[it.first] = false;
//				G1->G[it.first].hm[i] = false;
//			}
//			else {
//				printf("");
//			}
//		}
//	}
//	return changed;
//}

void couplednetwork::gaint_component(network *G1, network *G2)
{								  
	int label = 1;
	int maxsize = 0, maxcluster = -1, size = 0;//monosize=0, dimersize=0;
	int secondsize = 0, secondcluster = -1;
	DWORD t1 = GetTickCount();
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
	DWORD t2 = GetTickCount();

	G1->maxcc_size = maxsize;
	G1->secondcc_size = secondsize;
	//set non gaint component vertex to death
	for (int j = 0; j < N; j++) {
		if (G1->G[j].alive) {
			if (G1->G[j].cluster != maxcluster && G1->G[j].cluster != secondcluster){ //if cluster label is not max cluster
				G1->G[j].alive = false;
				if (G1->G[j].interdependent != -1)	 //-1 means has no interdependent nodes
					G2->G[G1->G[j].interdependent].alive = false;
				//img_print(A, true);
				//img_print(B, false);
			}
			G1->G[j].cluster = 0;
		}
	}
	DWORD t3 = GetTickCount();
	//printf("dfs:%ld\t%ld\n", t2 - t1, t3 - t2);
}

double average_degree(network *net)
{
	int degree = 0;
	for (int i = 0; i < N; i++)
	{
		/*for (auto it : net->G[i].hm) {
			if (it.second)
				degree++;
		}*/
		degree += net->G[i].interconnect->stack_len();
		
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
	sprintf(filename, "data/result_%d_%d_%d_%d.dat", N, ER_AVERAGE_DEGREE, minl, maxl);
	FILE *fp = fopen(filename, "w");
	if (!fp)
		printf("Bad file or directory!");

	couplednetwork cn;

	cn.A->lattice(); 
	for (double c = minl; c < maxl && c < LATTICE_SIZE; c++) {
		//cn.A->ER_length1((double)ER_AVERAGE_DEGREE / N, c);
		//network B
		cn.B->ER_length1((double)ER_AVERAGE_DEGREE / N, c);
		//cn.B->random_geometric_graph(0.5, c);
		//cn.A->random_geometric_graph(0.5, c);
		//cn.A->G = cn.B->G;
		//cn.B->scale_free_length(c);
		//cn.B->configuration_model(3.0, c)		
		double avg_d = average_degree(cn.B);
		for (double p = 0.71; p <= 1.0; p += 0.01) {
			for (int k = 0; k < NSAMPLE; k++) {
				//get random targets
				//cn.random_couple();
				cn.A->init();
				cn.B->init();
				cn.init_attack(cn.A, cn.B, 1-p);
				//printf("%f\n", average_degree(cn.B));
				int pre_cluster_size = 0, cluster_size = 0;
				int iter = 0;

				while (1) {
					iter++;
					//largest cc in A live
					//long t1 = clock();
					//
					//LARGE_INTEGER  large_interger;
					//double dff;
					//__int64  c1, c2;
					//QueryPerformanceFrequency(&large_interger);
					//dff = large_interger.QuadPart;
					//QueryPerformanceCounter(&large_interger);
					//c1 = large_interger.QuadPart;

					cn.gaint_component(cn.A, cn.B);
					//if (!cn.attack(cn.A, cn.B))
					//	break;
					//QueryPerformanceCounter(&large_interger);
					//c2 = large_interger.QuadPart;
					//printf("本机高精度计时器频率%lf\n", dff);
					//printf("第一次计时器值%I64d 第二次计时器值%I64d 计时器差%I64d\n", c1, c2, c2 - c1);
					//printf("%d\t%lf\t%d\n", k, (c2 - c1) * 1000 / dff, cn.A->maxcc_size);

					//long t2 = clock();
					cluster_size = cn.A->maxcc_size;
					if (cluster_size == pre_cluster_size)
						break;
					pre_cluster_size = cluster_size;
					//largest cc in B live
					cn.gaint_component(cn.B, cn.A);	
					//if (!cn.attack(cn.B, cn.A))
					//	break;
					//long t3 = clock();
					//printf("%ld\t%ld\n", t2 - t1, t3 - t2);
				}
				
				//cn.gaint_component(cn.A, cn.B);
				//printf("%d\t%f\t%f\t%d\t%d\t%d\t%d\n", k, c, p, cluster_size, s.monosize, s.dimersize, iter);
				//printf("%d\t%f\t%f\t%d\t%d\t%d\n", k, c, p, cn.A->maxcc_size, cn.A->secondcc_size, iter);
				fprintf(fp, "%d\t%f\t%f\t%d\t%d\t%d\n", k, c, p, cn.A->maxcc_size, cn.A->secondcc_size, iter);
				//printf("%d\t%f\t%f\t%d\t%d\t%d\n", k, c, p, cn.A->maxcc_size, cn.A->secondcc_size, iter);
				//fprintf(fp, "%d\t%f\t%f\t%d\t%f\t%d\n", k, c, p, cn.A->maxcc_size, avg_d, iter);
				
				
			}
			printf("%.3f\t%.3f\t%.3f\n", c, p, avg_d);
		
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
	vector<thread*> tid;
	for (int i = 0; i < NTHREAD; i++) {
		//memset(&tid[i], 0, sizeof(tid)); 
		//int err = pthread_create(&tid[i], NULL, run, (void *)(args + i));
		//if (err != 0)	{
		//	printf("create thread error: %s\n", strerror(err));
		//	exit(0);
		//}
		//printf("create thread: %d\n", tid[i]);
		thread *t = new thread(run, &args[i]);
		//t.join();
		tid.push_back(t);
	}

	//for (int i = 0; i < NTHREAD; i++) {
	//	pthread_join(tid[i], NULL);
	//	//printf("waiting %d", tid[i]);	
	//}
	for (auto t : tid) {
		if (t->joinable())
			t->join();
	}
	//if (args) {
	//	free(args);
	//	args = NULL;
	//}

	//getchar();
	return 0;
}