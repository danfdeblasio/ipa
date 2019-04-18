#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "error.h"
#include "extend.h"

extern char opt_debg;
extern char opt_term_gaps;

/*****************************************************************************
 * Simple stack implementation.  
 * This stack will be used as a candidate list. 
 *****************************************************************************/
void create_stack(type_stack *S, int n) {
  S->list = (type_node *)malloc(n*sizeof(type_node)); 
  assert(S->list); 
  S->top = 0; 
}
void delete_stack(type_stack *S) { 
  free(S->list); 
  S->top = 0;
}
void reset_stack(type_stack *S) { 
  S->top = 0; 
}
int empty_stack(type_stack *S) { 
  return S->top == 0; 
}
int full_stack(type_stack *S) { 
  return S->top == INT_MAX; 
}
void push_stack(type_stack *S, const type_node *e) { 
  memcpy(&(S->list[S->top++]), e, sizeof(type_node)); 
}
void pop_stack(type_stack *S, type_node *e) { 
  memcpy(e, &(S->list[--S->top]), sizeof(type_node)); 
}

void create_table(type_table *T, int row, int col) {
  /* Create stacks for the candidate forwarding. */
  create_stack(&(T->rs), col+1);
  T->cs = (type_stack *)calloc(col+1, sizeof(type_stack)); 
  assert(T->cs);
  for (int i = 0; i <= col; i++) {
    create_stack(&(T->cs[i]), row+1); 
    assert(&(T->cs[i]));
  }
  /* Allocate cells in each entry of the table. */
  double **V = (double **)malloc((row+1)*sizeof(double *)); assert(V);
  double **H = (double **)malloc((row+1)*sizeof(double *)); assert(H);
  double **D = (double **)malloc((row+1)*sizeof(double *)); assert(D);
  char **v   = (char **)malloc((row+1)*sizeof(char *)); assert(v);
  char **t   = (char **)malloc((row+1)*sizeof(char *)); assert(t);
  char **hV  = (char **)malloc((row+1)*sizeof(char *)); assert(hV);
  char **hH  = (char **)malloc((row+1)*sizeof(char *)); assert(hH);
  int **iV   = (int **)malloc((row+1)*sizeof(int *)); assert(iV);
  int **iH   = (int **)malloc((row+1)*sizeof(int *)); assert(iH);
  for (int i = 0; i <= row; i++) {
    V[i]  = (double *)malloc((col+1)*sizeof(double)); assert(V[i]);
    H[i]  = (double *)malloc((col+1)*sizeof(double)); assert(H[i]);
    D[i]  = (double *)malloc((col+1)*sizeof(double)); assert(D[i]);
    v[i]  = (char *)malloc((col+1)*sizeof(char)); assert(v[i]);
    t[i]  = (char *)malloc((col+1)*sizeof(char)); assert(t[i]);
    hV[i] = (char *)malloc((col+1)*sizeof(char)); assert(hV[i]);
    hH[i] = (char *)malloc((col+1)*sizeof(char)); assert(hH[i]);
    iV[i] = (int *)malloc((col+1)*sizeof(int)); assert(iV[i]);
    iH[i] = (int *)malloc((col+1)*sizeof(int)); assert(iH[i]);
  }  
  T->row = row; 
  T->col = col;
  T->V = V; 
  T->H = H; 
  T->D = D;
  T->v = v; 
  T->t = t;
  T->hV = hV; 
  T->hH = hH; 
  T->iV = iV; 
  T->iH = iH;
}

void delete_table(type_table *T) {
  /* Free stacks for the candidate forwarding. */
  delete_stack(&(T->rs));
  for (int i = 0; i <= T->col; i++) 
    delete_stack(&(T->cs[i]));
  free(T->cs);
  /* Free cells in each entry. */
  for (int i = 0; i <= T->row; i++) {
    free(T->V[i]);  
    free(T->H[i]);  
    free(T->D[i]);
    free(T->v[i]);  
    free(T->t[i]); 
    free(T->hV[i]); 
    free(T->hH[i]);
    free(T->iV[i]); 
    free(T->iH[i]);
  }
  T->row = T->col = 0;
  free(T->V);  
  free(T->H);  
  free(T->D);
  free(T->v);  
  free(T->t);
  free(T->hV); 
  free(T->hH);
  free(T->iV); 
  free(T->iH);
  T->V = T->H = T->D = 0;
  T->v = T->t = T->hV = T->hH = 0;
  T->iV = T->iH = 0;
}

void get_positions_table
(type_table *T, unsigned int *LA, unsigned int *LB, unsigned int m, unsigned int n, int *i) {
  int l = *i;
  while (TRACK_END != T->t[m][n]) {
    switch (T->t[m][n]) {
    case TRACK_SUB:
      LA[*i] = m-1; LB[(*i)++] = n-1; 
      m--; n--;
      break;
    case TRACK_DEL:
      for (int k = m; k >= T->iV[m][n]+1; k--) 
	{ LA[*i] = k-1; LB[(*i)++] = -1; }
      m = T->iV[m][n];
      break;
    case TRACK_INS:
      for (int k = n; k >= T->iH[m][n]+1; k--) 
	{ LA[*i] = -1; LB[(*i)++] = k-1; }
      n = T->iH[m][n];
      break;
    default: break;
    }
  }
  /* alright, the table has an optimal alignment of two empty strings. */
  if (*i == l) return;
  /* otherwise, reverse the order. */
  for (int k = l; k < ((*i-1+l)>>1)+((*i-l)%2 == 0); k++) {
    LA[k] ^= LA[*i-1-k+l]; LA[*i-1-k+l] ^= LA[k]; LA[k] ^= LA[*i-1-k+l]; 
    LB[k] ^= LB[*i-1-k+l]; LB[*i-1-k+l] ^= LB[k]; LB[k] ^= LB[*i-1-k+l];
  }
}

void write_string(FILE *stream, char *A, unsigned int *S, unsigned int len) {
  char *letter = "ARNDCQEGHILKMFPSTWYVBZX";
  unsigned int k;
  for (k = 0; k < len; k++) 
    fprintf(stream, "%c", S[k] == (unsigned int)-1 ? '-' : letter[(int)A[S[k]]]);
  printf("\n");
}


void 
complete_examples
(
 type_table *T, 
 double **H, 
 const alignment_t *S,
 const alignment_t *confidences, 
 const alignment_t *predictions,
 const dist_t ***distributions, 
 const match_t *M,
 const block_t *B, 
 unsigned int num_examples, 
 const double *P, 
 double **F
 ) {
  extern int (*work_the_oracle)
    (type_table *, double **, const alignment_t *,
     const alignment_t *, const alignment_t *, const dist_t ***, 
     const match_t *, const double *, unsigned int, unsigned int, unsigned int, int, int, 
     double, const double *, double **, char, char, char, double *);
  extern unsigned int feature_size;

  double *V = (double *)malloc(feature_size*sizeof(double)); 
  assert(V != NULL);

  for (int i = 0; i < (int)num_examples; i++) {
    /* Save the current cost of the example. */
    double old_cost = 0;
    for (int j = 0; j < (int)feature_size; j++)
      old_cost += F[i][j]*P[j];
    memcpy(F[i], B[i].fea, feature_size*sizeof(double));

    /* Complete the example optimally. */
    for (int j = 0; j < (int)B[i].count; j++) {
      char bg = 0;
      char eg = 0;
      if (B[i].pos[j].x == 0 ||
	  B[i].pos[j].y == 0) 
	/* The condition means that one of substrings to be aligned is a prefix. */
	bg = opt_term_gaps;
      if (B[i].len[j].x+B[i].pos[j].x == S[M[i].a].length[M[i].s] ||
	  B[i].len[j].y+B[i].pos[j].y == S[M[i].a].length[M[i].t]) 
	/* The condition means that one of substrings to be aligned is a suffix. */
	eg = opt_term_gaps;
      
      (*work_the_oracle)
	(T, H, S, confidences, predictions, distributions,
	 M, 0, i, B[i].pos[j].x, B[i].pos[j].y, 
	 B[i].len[j].x+B[i].pos[j].x-1,
	 B[i].len[j].y+B[i].pos[j].y-1, 0, P, &V, 0, bg, eg, 0);
   
      for (int k = 0; k < (int)feature_size; k++)
	F[i][k] += V[k];
    }

    double new_cost = 0;
    for (int j = 0; j < (int)feature_size; j++)
      new_cost+= F[i][j]*P[j];
    if (opt_debg) {
      fprintf(stderr, "example %d: old cost: %.4e >= new cost: %.4e", 
	      i, old_cost, new_cost);
      if (old_cost+1e-12 >= new_cost) 
	fprintf(stderr, "\n");
      else error("\n");
    }
    if (new_cost > old_cost+1e-12)
      error("the old score(%f) is better than the new score(%f).\n", 
	    old_cost, new_cost);
  }
  free(V);
}
