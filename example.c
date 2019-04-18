#include <ctype.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "inverse.h"
#include "fasta.h"
//#include "nonstd.h"
#include "example.h"
 
#define NUM_PAIRSO(a)   ((a)*((a)-1)>>1)
#define NUM_PAIRSW(a)   ((a)*((a)+1)>>1)
#define FLAT_INDEX(a,b) (((a)>(b))?(NUM_PAIRSW(a)+(b)):(NUM_PAIRSW(b)+(a)))

void 
get_alignments
(
 const char **path, 
 unsigned int num_alignments, 
 alignment_t **A
 ) 
{
  *A = (alignment_t *)malloc(num_alignments*sizeof(alignment_t));
  assert(*A);
  unsigned int i, j;
  for (i = 0; i < num_alignments; i++) {
    read_fasta(path[i], 1, &(*A)[i].name, &(*A)[i].string, &(*A)[i].count);
    (*A)[i].length = (unsigned int *)malloc((*A)[i].count*sizeof(unsigned int));
    assert((*A)[i].length);
    for (j = 0; j < (*A)[i].count; j++)
      (*A)[i].length[j] = strlen((*A)[i].string[j]);
  }
}

void free_alignments
(alignment_t **A, unsigned int num_alignments)
{
  unsigned int i, j;
  for (i = 0; i < num_alignments; i++) {
    for (j = 0; j < (*A)[i].count; j++) {
      free((*A)[i].name[j]);
      free((*A)[i].string[j]);
    }
    free((*A)[i].name);
    free((*A)[i].string);
    free((*A)[i].length);
  }
  free(*A);
  *A = NULL;
}

void make_alphabet(const char *letter, const char *space, char **alphabet) {
  unsigned int length = strlen(letter)+strlen(space);
  *alphabet = (char *)malloc((length+1)*sizeof(char));
  assert(*alphabet);
  strcpy(*alphabet, letter);
  strcat(*alphabet, space);
  (*alphabet)[length] = 0;
}

void number_alphabet(const char *letter, const char *space, int **number) {
  *number = (int *)malloc(CHAR_MAX*sizeof(int));
  assert(*number);
  memset(*number, -1, CHAR_MAX*sizeof(int));
  unsigned int i, size = strlen(letter);
  for (i = 0; i < size; i++) 
    (*number)[toupper(letter[i])] = i;
  for (i = 0; i < strlen(space); i++) 
    (*number)[toupper(space[i])] = size;
}

void convert_alignments
(const int *map, const alignment_t *A, unsigned int num_alignments, char space, alignment_t **S)
{
  *S = (alignment_t *)malloc(num_alignments*sizeof(alignment_t));
  assert(*S);
  unsigned int i, j;
  for (i = 0; i < num_alignments; i++) {
    (*S)[i].count  = A[i].count;
    (*S)[i].name   = (char **)malloc(A[i].count*sizeof(char *));
    assert((*S)[i].name);
    (*S)[i].string = (char **)malloc(A[i].count*sizeof(char *));
    assert((*S)[i].string);
    (*S)[i].length = (unsigned int *)malloc(A[i].count*sizeof(unsigned int));
    assert((*S)[i].length);
    for (j = 0; j < A[i].count; j++) {
      (*S)[i].name[j] = (char *)malloc(strlen(A[i].name[j])+1);
      assert((*S)[i].name[j]);
      strcpy((*S)[i].name[j], A[i].name[j]);
      //(*S)[i].name[j]   = strdup(A[i].name[j]);
      (*S)[i].length[j] = A[i].length[j];
      (*S)[i].string[j] = (char *)malloc(A[i].length[j]*sizeof(char));
      assert((*S)[i].string[j]);
      unsigned int k, l;
      for (l = k = 0; k < A[i].length[j]; k++) 
	if (A[i].string[j][k] == '-' && !space) continue;
	else (*S)[i].string[j][l++] = map[toupper(A[i].string[j][k])];
      (*S)[i].length[j] = l;
    }
  }
}

void get_matchlists
(const alignment_t *A, unsigned int num_alignments, match_t **M, unsigned int *num_examples)
{
  unsigned int i;
  for (*num_examples = i = 0; i < num_alignments; i++)
    *num_examples += NUM_PAIRSO(A[i].count);
  *M = (match_t *)malloc(*num_examples*sizeof(match_t));
  assert(*M);

  unsigned int j, k, l;
  for (l = i = 0; i < num_alignments; i++) 
    for (j = 0; j < A[i].count; j++)
      for (k = 0; k < j; k++) {
	(*M)[l].a = i;
	(*M)[l].s = j;
	(*M)[l].t = k;
	(*M)[l].length = strlen(A[i].string[0]);
	(*M)[l].list = (point_t *)malloc((*M)[l].length*sizeof(point_t));
	assert((*M)[l].list);
	unsigned int m, n, p1, p2;
	p1 = p2 = m = 0;
	for (n = 0; n < (*M)[l].length; n++) {
	  if (isupper(A[i].string[j][n]) && isupper(A[i].string[k][n])) {
	    (*M)[l].list[m].x = p1;
	    (*M)[l].list[m].y = p2;
	    m++;
	  }
	  if (A[i].string[j][n] != '-') p1++;
	  if (A[i].string[k][n] != '-') p2++;
	}
	(*M)[l].length = m;
	(*M)[l].list = (point_t *)realloc
	  ((*M)[l].list, (m ? m : 1)*sizeof(point_t));
	assert((*M)[l].list);
	l++;
      }
}

void free_matchlists
(match_t **M, unsigned int num_blocks)
{
  unsigned int i;
  for (i = 0; i < num_blocks; i++)
    free((*M)[i].list);
  free(*M);
  *M = NULL;
}

//extern switch_t sw;
extern unsigned int feature_size;
extern unsigned int gap_offset, len_offset, sub_offset, add_offset;

void
get_coreblocks
(const alignment_t *S, /* sequence */
 const alignment_t *C, /* prediction confidence on S */
 const alignment_t *P, /* single prediction on S */
 const dist_t ***D,    /* multiple prediction on S */
 const match_t *M, 
 unsigned int num_examples, 
 block_t **B)
{
  *B = (block_t *)malloc(num_examples*sizeof(block_t));
  assert(*B);  

  for (unsigned int i = 0; i < num_examples; i++) {
    /* These two pointers are notational shortcuts. */
    block_t *block = &(*B)[i];
    point_t *list = M[i].list;
    /* Find the number of non-core blocks or unreliable regions. */
    unsigned int count = 0;
    if (M[i].length == 0) 
      count = 1;
    else {
      if (list[0].x != 0 || 
	  list[0].y != 0) 
	count++;
      unsigned int j;
      for (j = 1; j < M[i].length; j++)
	if (list[j].x-list[j-1].x != 1 || 
	    list[j].y-list[j-1].y != 1) 
	  count++;
      if (S[M[i].a].length[M[i].s]-list[j-1].x != 1 ||
	  S[M[i].a].length[M[i].t]-list[j-1].y != 1) 
	count++;
    }
    /* gather information about core and non-core blocks. */
    block->count = count;
    block->pos = (point_t *)malloc(block->count*sizeof(point_t));
    assert(block->pos);
    block->len = (point_t *)malloc(block->count*sizeof(point_t));
    assert(block->len);
    block->fea = (double *)calloc(feature_size, sizeof(double));
    assert(block->fea);
    double *Q = &block->fea[sub_offset];
    double *R = &block->fea[add_offset];
    char *X = S[M[i].a].string[M[i].s];
    char *Y = S[M[i].a].string[M[i].t];
    /*
     * Process reliable regions.
     */
    /* Count substitution features inside core blocks */
    for (unsigned int j = 0; j < M[i].length; j++) 
      Q[FLAT_INDEX(X[list[j].x], Y[list[j].y])]++;
    extern char opt_prot_pred;
    if (opt_prot_pred == 1) {
      extern unsigned int add_size;
      //const char *CX = C[M[i].a].string[M[i].s];
      //const char *CY = C[M[i].a].string[M[i].t];
      const char *PX = P[M[i].a].string[M[i].s];
      const char *PY = P[M[i].a].string[M[i].t];
      for (unsigned int j = 0; j < M[i].length; j++) {
	unsigned int index = FLAT_INDEX(PX[list[j].x], PY[list[j].y]);
	if (index < add_size) 
	  R[index]++;
      }
    }
    else if (opt_prot_pred == 2) {
      const dist_t *DX = D[M[i].a][M[i].s];
      const dist_t *DY = D[M[i].a][M[i].t];
      for (unsigned int j = 0; j < M[i].length; j++) {
	R[0] += DX[list[j].x].a*DY[list[j].y].a;
	R[1] += DX[list[j].x].a*DY[list[j].y].b
	  +DX[list[j].x].b*DY[list[j].y].a;
	R[2] += DX[list[j].x].b*DY[list[j].y].b;
	R[3] += DX[list[j].x].a*DY[list[j].y].c
	  +DX[list[j].x].c*DY[list[j].y].a;
	R[4] += DX[list[j].x].b*DY[list[j].y].c
	  +DX[list[j].x].c*DY[list[j].y].b;
	R[5] += DX[list[j].x].c*DY[list[j].y].c;
      }
    }
    /*
     * Process unreliable regions. 
     */
    if (M[i].length == 0) {
      block->pos[0].x = 0;
      block->pos[0].y = 0;
      block->len[0].x = S[M[i].a].length[M[i].s]-1;
      block->len[0].y = S[M[i].a].length[M[i].t]-1;
    }
    else {
      count = 0;
      if (list[0].x != 0 || list[0].y != 0) {
	block->pos[count].x = 0;
	block->pos[count].y = 0;      
	block->len[count].x = list[0].x;
	block->len[count].y = list[0].y; 
	count++;
      }
      unsigned int j;
      for (j = 1; j < M[i].length; j++) {
	unsigned int m = list[j].x-list[j-1].x;
	unsigned int n = list[j].y-list[j-1].y;
	if (m == 1 && n == 1) 
	  continue;
	block->pos[count].x = list[j-1].x+1;
	block->pos[count].y = list[j-1].y+1;
	block->len[count].x = m-1;
	block->len[count].y = n-1; 
	count++;
      }
      unsigned int m = S[M[i].a].length[M[i].s]-list[j-1].x;
      unsigned int n = S[M[i].a].length[M[i].t]-list[j-1].y;
      if (m != 1 || n != 1) {
	block->pos[count].x = list[j-1].x+1;
	block->pos[count].y = list[j-1].y+1;
	block->len[count].x = m-1;
	block->len[count].y = n-1;
      }
    }
  }
}

void free_coreblocks
(block_t **B, unsigned int num_blocks)
{
  unsigned int i;
  for (i = 0; i < num_blocks; i++) {
    free((*B)[i].fea);
    free((*B)[i].pos);
    free((*B)[i].len);
  }
  free(*B);
  *B = NULL;
}

void get_wholeblocks
(const alignment_t *A, unsigned int num_alignments, const int *alphabet, double ***F, unsigned int *num_examples)
{
  unsigned int i;
  for (*num_examples = i = 0; i < num_alignments; i++)
    *num_examples += NUM_PAIRSO(A[i].count);
  *F = (double **)malloc(*num_examples*sizeof(double *));
  assert(*F);
  for (i = 0; i < *num_examples; i++) {
    (*F)[i] = (double *)calloc(feature_size, sizeof(double));
    assert((*F)[i]);
  }
  unsigned int j, k, l, m;
  for (m = i = 0; i < num_alignments; i++)
    for (j = 0; j < A[i].count; j++)
      for (k = 0; k < j; k++) {
	char edge = 0;
	double *G = &(*F)[m][gap_offset];
	double *S = &(*F)[m][sub_offset];
	for (l = 0; l < A[i].length[0]; l++) 
	  if (A[i].string[j][l] != '-') {
	    if (A[i].string[k][l] != '-') {
	      edge = 0; 
	      unsigned int a = alphabet[toupper(A[i].string[j][l])];
	      unsigned int b = alphabet[toupper(A[i].string[k][l])];
	      S[FLAT_INDEX(a, b)]++;
	    }
	    else 
	      { G[1]++; if (edge != 1) G[0]++; edge = 1; }
	  }
	  else {
	    if (A[i].string[k][l] != '-') 
	      { G[1]++; if (edge != 2) G[0]++; edge = 2; }
	  }
	m++;
      }
}

void get_wholeblocks_2
(const alignment_t *A, const alignment_t *T, unsigned int num_alignments, double ***F, unsigned int *num_examples)
{
  extern char opt_prot_pred;

  *num_examples = 0;
  for (int i = 0; i < num_alignments; i++)
    *num_examples += NUM_PAIRSO(A[i].count);
  *F = (double **)malloc(*num_examples*sizeof(double *));
  assert(*F);
  for (int i = 0; i < *num_examples; i++) {
    (*F)[i] = (double *)calloc(feature_size, sizeof(double));
    assert((*F)[i]);
  }
  /* remove this. */
  return;

  /*
  unsigned int m = 0;
  for (int i = 0; i < num_alignments; i++)
    for (int j = 0; j < A[i].count; j++)
      for (int k = 0; k < j; k++) {
	double *G = &(*F)[m][gap_offset];
	double *L = &(*F)[m][len_offset];
	double *S = &(*F)[m][sub_offset];
	double *R = &(*F)[m][add_offset];
	extern unsigned int opt_disc_lvls;
	char terminal = opt_disc_lvls;
	int s = 0;
	int t = 0;
	int n = 0;
	for (int l = 0; l < A[i].length[0]; l++) {
	  if (A[i].string[j][l] != '-')
	    if (A[i].string[k][l] != '-') { 
	    // Count substitutions. 
	      char a = T[i].string[j][s];
	      char b = T[i].string[k][t];
	      S[FLAT_INDEX(a, b)]++;
	      // Count .
	      switch (opt_prot_pred) {
	      case 1:
		extern unsigned int add_size;

		//const char *CX = C[M[i].a].string[M[i].s];
		//const char *CY = C[M[i].a].string[M[i].t];
		const char *PX = P[M[i].a].string[M[i].s];
		const char *PY = P[M[i].a].string[M[i].t];
		


		for (unsigned int j = 0; j < M[i].length; j++) {
		  unsigned int index = FLAT_INDEX(PX[list[j].x], PY[list[j].y]);
		  if (index < add_size) 
		    R[index]++;
		}

		break;
	      case 2:
		const dist_t *DX = D[M[i].a][M[i].s];
		const dist_t *DY = D[M[i].a][M[i].t];
		R[0] += DX[list[j].x].a*DY[list[j].y].a;
		R[1] += DX[list[j].x].a*DY[list[j].y].b
		  +DX[list[j].x].b*DY[list[j].y].a;
		R[2] += DX[list[j].x].b*DY[list[j].y].b;
		R[3] += DX[list[j].x].a*DY[list[j].y].c
		  +DX[list[j].x].c*DY[list[j].y].a;
		R[4] += DX[list[j].x].b*DY[list[j].y].c
		  +DX[list[j].x].c*DY[list[j].y].b;
		break;
	      }
	      s++;
	      t++;
	    }
	    else { 
	    // deletion 
	      for (n = l; n < A[i].length[0]; n++) 
		if (A[i].string[j][n] == '-' || A[i].string[k][n] != '-') 
		  break;
	      n--;
	      s += n-l+1;
	      if (l == 0 || n == A[i].length[0]-1) {
		G[(int)terminal]++;
		L[(int)terminal] += n-l+1;
	      }
	      else {
		G[(int)query_max(&T[i].heap[j], s-1-(n-l), s-1)]++;
		for (; l <= n; l++) L[(int)T[i].hydstr[j][s-1-(n-l)]]++;
	      }
	      l = n;
	    }
	  else
	    if (A[i].string[k][l] != '-') { 
	      // insertion 
	      for (n = l; n < A[i].length[0]; n++) 
		if (A[i].string[j][n] != '-' || A[i].string[k][n] == '-') 
		  break;
	      n--;
	      t += n-l+1;
	      if (l == 0 || n == A[i].length[0]-1) {
		G[(int)terminal]++;
		L[(int)terminal] += n-l+1;
	      }
	      else {
		G[(int)query_max(&T[i].heap[k], t-1-(n-l), t-1)]++;
		for (; l <= n; l++) L[(int)T[i].hydstr[k][t-1-(n-l)]]++;
	      }
	      l = n;
	    }
	    else ;
	}
	m++;
      }
  */
}

void free_features
(double ***F, int num_features)
{
  if (!*F) return;
  for (int i = 0; i < num_features; i++) 
    free((*F)[i]);
  free(*F);
  *F = 0;
}
