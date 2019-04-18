#include <float.h>
#include <stdlib.h>
#include <string.h>
#include "error.h"
#include "extend.h"

extern char   opt_debg;
extern char   opt_disr_lvls;
extern char   opt_norm_cost;
extern char   opt_prot_pred;
extern double opt_conv_wght;
extern double threshold_v;
extern unsigned int feature_size, gap_offset, len_offset, sub_offset, add_offset;
//extern double mean_length;
extern unsigned int add_size;

static void get_features_table
(type_table *T, const char *A, const char *B, const char *HA, const char *HB,
 const char *PA, const char *PB, const dist_t *DA, const dist_t *DB,
 unsigned int m, unsigned int n, double *F, int *len) {
  if (len) *len = 0;
  memset(F, 0, feature_size*sizeof(double));
  double *G = &F[gap_offset]; 
  double *L = &F[len_offset];
  double *S = &F[sub_offset];
  double *R = &F[add_offset];
  const char term = opt_disr_lvls;
  while (TRACK_END != T->t[m][n]) {
    switch (T->t[m][n]) {
    case TRACK_SUB:
      if (T->v[m][n] && len) (*len)++;
      unsigned int index = FLAT_INDEX(A[m], B[n]);
      S[index]++;
      switch (opt_prot_pred) {
      case 1: 
	index = FLAT_INDEX(PA[m], PB[n]);
	if (index < add_size) R[index]++;
	break;
      case 2:
	R[0] += DA[m].a*DB[n].a;
	R[1] += DA[m].a*DB[n].b+DA[m].b*DB[n].a;
	R[2] += DA[m].b*DB[n].b;
	R[3] += DA[m].a*DB[n].c+DA[m].c*DB[n].a;
	R[4] += DA[m].b*DB[n].c+DA[m].c*DB[n].b;
	R[5] += DA[m].c*DB[n].c;
	break;
      }
      m--; n--;
      break;
    case TRACK_DEL:
      G[(int)T->hV[m][n]]++;
      if (T->hV[m][n] == term)
	L[(int)term] += m-T->iV[m][n];
      else
	for (unsigned int i = T->iV[m][n]+1; i <= m; i++)
	  //L[(int)HA[i]]++;
	  L[0]++;
      m = T->iV[m][n];
      break;
    case TRACK_INS:
      G[(int)T->hH[m][n]]++;
      if (T->hH[m][n] == term)
	L[(int)term] += n-T->iH[m][n];
      else 
	for (unsigned int i = T->iH[m][n]+1; i <= n; i++)
	  //L[(int)HB[i]]++;
	  L[0]++;
      n = T->iH[m][n];
      break;
    default: ;
    }
  }
}

static inline double cost
(double s, const double *L, unsigned int i, unsigned int j) {
  //return s+query_sum(L, i-1, j-1);
  return s+L[0]*(j-i+1);
}

static inline void forward
(char t, double **D, const double *L, type_stack *S, unsigned int i, unsigned int j, unsigned int n) {
  /* forward candidate j to other entries if t = 0 (i.e. a row is fixed).
     otherwise, forward candidate i to other entries if t != 0. */
  type_node e;
  memcpy(&e, &S->list[S->top-1], sizeof(type_node));
  unsigned int k = t ? i : j;
  double subst_cost = t ? D[e.b][j] : D[i][e.b];
  unsigned int last_known = 0;
  if (cost(D[i][j], L, k+1, k+1) < 
      cost(subst_cost, L, e.b+1, k+1)) {
    last_known = k+1;
    /* whenever the current candidate wins against the candidate on the
       stack at the top, pop the stack. */
    while (!empty_stack(S) &&
	   cost(D[i][j], L, k+1, e.p) <
	   cost(subst_cost, L, e.b+1, e.p)) { 
      last_known = e.p;
      pop_stack(S, &e);
      if (!empty_stack(S)) {
	memcpy(&e, &S->list[S->top-1], sizeof(type_node));
	subst_cost = t ? D[e.b][j] : D[i][e.b];
      }
    }
    if (empty_stack(S)) e.p = n;
    else {
      /* lo is the known position at which the current candidate wins and
	 hi is the known position at which the current candidate loses. */
      unsigned int lo = last_known;
      unsigned int hi = S->list[S->top-1].p;
      while (lo+1 < hi) {
	unsigned int md = (lo+hi)>>1;
	if (cost(D[i][j], L, k+1, md) <
	    cost(subst_cost, L, e.b+1, md)) lo = md;
	else hi = md;
      }
      e.p = lo;
    }
    e.b = k;
    push_stack(S, &e);
  }
  /* check the validity of stack. */
  if (S->list[0].p != n) error("stack is corrupted.");
  for (k = 0; k < S->top-1; k++) {
    if (S->list[k].b >= S->list[k+1].b) error("stack is corrupted.");
    if (S->list[k].p <= S->list[k+1].p) error("stack is corrupted.");
  }
}

int work_the_oracle_3
(type_table *T, double **O, const alignment_t *S,
 const alignment_t *confidences, const alignment_t *predictions,
 const dist_t ***distributions, const match_t *M,
 const double *E, unsigned int id, unsigned int as, unsigned int bs, int ae, int be, double delta,
 const double *P, double **F, char flag_sep, char bg, char eg, double *rate) {
  /* flag_sep = 0 : best alignment
     flag_sep = 1 : separation */
  /************************************************************** short cuts */
  /* alignment parameters */  
  const double *G = &P[gap_offset];  /* gap-open      costs */
  const double *L = &P[len_offset];  /* gap-extension costs */
  const double *Q = &P[sub_offset];  /* substitution  costs */
  const double *R = &P[add_offset];  /* modifier      costs */
  const char term = opt_disr_lvls;
  /* substrings to be aligned */
  const char  *A = S[M[id].a].string[M[id].s]-1;
  const char  *B = S[M[id].a].string[M[id].t]-1;

  const dist_t *DA = 0;
  const dist_t *DB = 0;
  const char *PA = 0;
  const char *PB = 0;
  switch (opt_prot_pred) {
  case 1:
    PA = predictions[M[id].a].string[M[id].s]-1;
    PB = predictions[M[id].a].string[M[id].t]-1;
    break;
  case 2:
    DA = distributions[M[id].a][M[id].s]-1;
    DB = distributions[M[id].a][M[id].t]-1;
    break;
  }

  double norm_const = 1;
  if (opt_norm_cost == 1)
    norm_const = (S[M[id].a].length[M[id].s]+S[M[id].a].length[M[id].t])/2.0;
  else if (opt_norm_cost == 2)
    norm_const = MAX(S[M[id].a].length[M[id].s], S[M[id].a].length[M[id].t]);
  assert(norm_const > 0);

  double recv_const = 0;
  if (flag_sep == 1 && M[id].length > 0) 
    recv_const = norm_const*(1/opt_conv_wght-1)/M[id].length;
  assert(recv_const >= 0);
  
  /*************************************** (1) data structure initialization */
  /* short cuts for the dynamic programming table */
  double **V = T->V;    /* cost ending with deletion at the cell */
  double **H = T->H;    /* cost ending with insertion at the cell */
  double **D = T->D;    /* minimum cost ending at the cell */
  char   **v = T->v;    /* 1 if the substitution is in a core block */
  char   **t = T->t;    /* track to cell giving the cheapest cost */
  char  **hV = T->hV;  
  char  **hH = T->hH;
  int  **iV = T->iV;
  int  **iH = T->iH;
  /* initialize the origin. */
  int i, j;
  V[as][bs] = DBL_MAX;
  H[as][bs] = DBL_MAX;
  D[as][bs] = 0;
  t[as][bs] = TRACK_END;
  /* initialize the first column boundary entries. */
  for (i = as+1; i <= ae+1; i++) {
    iV[i][bs] = as; 
    /* If terminal gap, else nonterminal gap */
    if (bg || (eg && (int)bs > be)) {
      hV[i][bs] = term;
      V[i][bs]  = G[(int)term]+L[(int)term]*(i-(as+1)+1);
    }
    else { 
      hV[i][bs] = 0;
      V[i][bs]  = G[(int)hV[i][bs]]+L[0]*(i-(as+1)+1);
    }
    iH[i][bs] = as; hH[i][bs] = 0; H[i][bs] = DBL_MAX;
    D[i][bs] = V[i][bs]; t[i][bs] = TRACK_DEL;
  }
  /* initialize the first row boundary entries. */
  for (j = bs+1; j <= be+1; j++) {
    iV[as][j] = bs; hV[as][j] = 0; V[as][j] = DBL_MAX;
    iH[as][j] = bs; 
    /* If terminal gap, else nonterminal gap */
    if (bg || (eg && (int)as > ae)) {
      hH[as][j] = term;
      H[as][j]  = G[(int)term]+L[(int)term]*(j-(bs+1)+1);
    }
    else { 
      hH[as][j] = 0;
      H[as][j]  = G[(int)hH[as][j]]+L[0]*(j-(bs+1)+1);
    }
    D[as][j] = H[as][j]; t[as][j] = TRACK_INS;
  }
  /* initialize stacks for forwarding. */
  reset_stack(&(T->rs));
  for (i = bs; i <= be+1; i++) {
    type_node e;
    e.b = as; e.p = ae+1;
    reset_stack(&(T->cs[i]));
    push_stack(&(T->cs[i]), &e);
  }
  /* initialize the hit-miss table for recovery rate. */
  for (i = as+1; i <= ae+1; i++) 
    memset(&v[i][bs+1], 0, (be-bs+1)*sizeof(char));
  for (i = 0; i < (int)M[id].length; i++)
    v[M[id].list[i].x+1][M[id].list[i].y+1] = 1;

  /**************************************************** (3) string alignment */
  /* align two strings. */
  for (i = as+1; i <= ae+1; i++) {
    type_node e;
    e.b = bs; e.p = be+1;
    type_stack *RS = &(T->rs);
    reset_stack(RS);
    push_stack(RS, &e);
    for (j = bs+1; j <= be+1; j++) {
      type_stack *CS = &(T->cs[j]);
      /* Assertion: The column stack is not corrupted. */
      assert(CS->list[CS->top-1].p >= (unsigned int)i); 
      /* Assertion: The row stack is not corrupted. */
      assert(RS->list[RS->top-1].p >= (unsigned int)j);
      double x;
      /*--------------------------------------------- (3.1) cost computation */
      /* the cost of an alignment ending with a deletion */
      iV[i][j] = CS->list[CS->top-1].b;
      V[i][j]  = D[iV[i][j]][j];
      /* If terminal gap, else nonterminal gap */
      if (j == be+1 && eg) {
	hV[i][j] = term;
	V[i][j] += G[(int)term]+L[(int)term]*(i-(iV[i][j]+1)+1);
      }
      else {
	hV[i][j] = 0;
	V[i][j] += G[(int)hV[i][j]]+L[0]*(i-(iV[i][j]+1)+1);
      }
      /* the cost of an alignment ending with an insertion */
      iH[i][j] = RS->list[RS->top-1].b;
      H[i][j]  = D[i][iH[i][j]];
      /* If terminal gap, else nonterminal gap */
      if (i == ae+1 && eg) {
	hH[i][j] = term;
	H[i][j] += G[(int)term]+L[(int)term]*(j-(iH[i][j]+1)+1); 
      }
      else {
	hH[i][j] = 0;
	H[i][j] += G[(int)hH[i][j]]+L[0]*(j-(iH[i][j]+1)+1); 
      }
      /* the cost of an alignment ending with a substitution */
      unsigned int index = FLAT_INDEX(A[i], B[j]);
      x = D[i-1][j-1]+Q[index];
      switch (opt_prot_pred) {
      case 1:
	index = FLAT_INDEX(PA[i], PB[j]); 
	assert(index < add_size);
	x += R[index];
	break;
      case 2:
	x += R[0]* DA[i].a*DB[j].a+
	  R[1]*(DA[i].a*DB[j].b+DA[i].b*DB[j].a)+
	  R[2]* DA[i].b*DB[j].b+
	  R[3]*(DA[i].a*DB[j].c+DA[i].c*DB[j].a)+
	  R[4]*(DA[i].b*DB[j].c+DA[i].c*DB[j].b)+
	  R[5]* DA[i].c*DB[j].c;
	break;
      }
      if (v[i][j]) 
	x += recv_const;
      /* contest among three possible alignments ending at (i, j). */
      if (H[i][j] < V[i][j]) {
	if (H[i][j] < x) { D[i][j] = H[i][j]; t[i][j] = TRACK_INS; }
	else { D[i][j] = x; t[i][j] = TRACK_SUB; }
      } 
      else {
	if (V[i][j] < x) { D[i][j] = V[i][j]; t[i][j] = TRACK_DEL; }
	else { D[i][j] = x; t[i][j] = TRACK_SUB; }
      }
      /* ------------------------- (3.2) eliminate past candidates in stacks */
      while (!empty_stack(CS) && CS->list[CS->top-1].p <= (unsigned int)i) 
	pop_stack(CS, &e);
      /* Assertion: The column stack is not empty. */
      assert(i == ae+1 || !empty_stack(CS));
      while (!empty_stack(RS) && RS->list[RS->top-1].p <= (unsigned int)j) 
	pop_stack(RS, &e);
      /* Assertion: The row stack is not empty. */
      assert(j == be+1 || !empty_stack(RS));
      /*----------------------------------------- (3.3) candidate forwarding */
      /* forward candidate onto the current column (i.e. deletion). */
      if (i != ae+1) {
	if (j == be+1 && eg) {
	  memcpy(&e, &CS->list[CS->top-1], sizeof(type_node));
	  if (D[i][j] < D[e.b][j]+L[(int)term]*(i-(e.b+1)+1))
	    { CS->list[CS->top-1].b = i; }
	}
	else forward(1, D, L, CS, i, j, ae+1);
      }
      /* forward candidate onto the current row (i.e. insertion). */
      if (j != be+1) {
	if (i == ae+1 && eg) {
	  memcpy(&e, &RS->list[RS->top-1], sizeof(type_node));
	  if (D[i][j] <= D[i][e.b]+L[(int)term]*(j-(e.b+1)+1))
	    { RS->list[RS->top-1].b = j; }
	} 
	else forward(0, D, L, RS, i, j, be+1);
      }
    }
  }
  /************************************ 
   * 4. Recover an optimal alignment. *
   ************************************/
  double *feature = (double *)malloc(feature_size*sizeof(double));
  assert(feature != NULL);
  int length;
  //get_features_table(T, A, B, HA-1, HB-1, PA, PB, DA, DB, 
  get_features_table(T, A, B, 0, 0, PA, PB, DA, DB, 
		     ae+1, be+1, feature, &length);
  double recv_rate = M[id].length ? (double)length/M[id].length : 1;
  if (F) 
    memcpy(F[0], feature, feature_size*sizeof(double));
  if (rate) 
    *rate = recv_rate;
  /***************************
   * 5. Check the violation. *
   ***************************/
  int return_value = 0;
  if (flag_sep == 1) {
    double cost = 0;
    for (int i = 0; i < feature_size; i++)
      cost += (E[i]-feature[i])*P[i];
    double lhs = cost+norm_const*(1/opt_conv_wght-1)*(1-recv_rate);
    double rhs = norm_const/opt_conv_wght*delta;
    if (lhs > rhs+threshold_v)
      return_value = 1;
  }
  free(feature);
  return return_value;
}
