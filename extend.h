#ifndef __EXTEND__
#define __EXTEND__

#include "example.h"
#include <stdio.h>

#define TRACK_END 0
#define TRACK_DEL 1
#define TRACK_INS 2
#define TRACK_SUB 3
#define MAX(x,y) ((x)>(y)?(x):(y))
#define NUM_PAIRSW(a)   ((a)*((a)+1)>>1)
#define FLAT_INDEX(a,b) (((a)>(b))?(NUM_PAIRSW(a)+(b)):(NUM_PAIRSW(b)+(a)))

typedef struct tag_entry {
  char tck;
  int gap, len;
  double sub;
} type_entry;

typedef struct tag_node {
  unsigned int b, p;
} type_node;

typedef struct tag_stack {
  int top;
  type_node *list;
} type_stack;

typedef struct tag_table {
  int row;
  int col;
  double **V; 
  double **H;
  double **D;
  char   **v;
  /* candidate lists */
  type_stack rs;
  type_stack *cs;
  /* recovery path */
  char **t;
  char **hV;
  char **hH;
  int  **iV;
  int  **iH;
} type_table;

void create_stack(type_stack *S, int n);
void delete_stack(type_stack *S); 
void reset_stack(type_stack *S); 
int  empty_stack(type_stack *S); 
int  full_stack(type_stack *S);
void push_stack(type_stack *S, const type_node *e);
void pop_stack(type_stack *S, type_node *e);
void create_table(type_table *T, int row, int col);
void delete_table(type_table *T);
void get_positions_table
(type_table *T, unsigned int *LA, unsigned int *LB, unsigned int m, unsigned int n, int *i);
void write_string(FILE *stream, char *A, unsigned int *S, unsigned int len);
int work_the_oracle_1
(type_table *T, double **O, const alignment_t *S,
 const alignment_t *predictions, const alignment_t *confidences,
 const dist_t ***distributions, const match_t *M,
 const double *E, unsigned int id, unsigned int as, unsigned int bs, int ae, int be, double delta,
 const double *P, double **F, char type, char bg, char eg, double *rate);
int work_the_oracle_2
(type_table *T, double **O, const alignment_t *S,
 const alignment_t *predictions, const alignment_t *confidences,
 const dist_t ***distributions, const match_t *M,
 const double *E, unsigned int id, unsigned int as, unsigned int bs, int ae, int be, double delta,
 const double *P, double **F, char type, char bg, char eg, double *rate);
int work_the_oracle_3
(type_table *T, double **O, const alignment_t *S,
 const alignment_t *predictions, const alignment_t *confidences,
 const dist_t ***distributions, const match_t *M,
 const double *E, unsigned int id, unsigned int as, unsigned int bs, int ae, int be, double delta,
 const double *P, double **F, char type, char bg, char eg, double *rate);
void complete_examples
(type_table *T, double **H, const alignment_t *S,
 const alignment_t *confidences, const alignment_t *predictions,
 const dist_t ***distributions, const match_t *M,
 const block_t *B, unsigned int num_examples, const double *P, double **F);

#endif
