
#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

static const int LATTICE_SIZE = 800;
static const int N = LATTICE_SIZE*LATTICE_SIZE;
#define NSAMPLE 10	 
#define MONOMER 0
#define DIMER 1

#define MIN(a,b) ((a)<(b)?(a):(b))


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
/******************************************************************/


typedef struct gcsize
{
	int maxsize;
	int monosize;
	int dimersize;
}gcsize;

typedef struct node
{
	//point id;
	int inter;
	stack *base;
	bool alive;
	int cluster;
	//int type;
	//int gaint;
}node;
