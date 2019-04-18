#include <math.h>
#include <limits.h>
#include <ctype.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <glpk.h>
#include "error.h"
#include "matrix.h"
#include "extend.h"
#include "inverse.h"

/***************************************************************************** 
 * Global variable 
 *****************************************************************************/
/* visible controllable options for IPA */
char   opt_debg = 0;
char   opt_stat = 0;
char   opt_verb = 0;
char   opt_part_exam = 1;
char   opt_alph_type = 0;
char   opt_cont_gaps = 3;
char   opt_disr_meas = 0;
char   opt_disr_lvls = 1;
char   opt_fixd_subs = 0;
char   opt_init_triv = 0;
char  *opt_just_algn = 0;
char   opt_maxm_itrs = 10;
char   opt_norm_cost = 0;
char   opt_prot_pred = 0;
char   opt_term_gaps = 0;
char   opt_wind_size = 1;
double opt_conv_wght = 1;
/* Invisible controllable options for IPA */
int   num_bests  = 2;
unsigned int   opt_cull_fold = 2;

static double mean_score = 0;
static double mean_recv  = 0;

int (*work_the_oracle)
(type_table *, double **, const alignment_t *,
 const alignment_t *, const alignment_t *, const dist_t ***, 
 const match_t *, const double *, unsigned int, unsigned int, unsigned int, int, int, 
 double, const double *, double **, char, char, char, double *) = 0;

double threshold_v = 1e-4;
const unsigned int threshold_e = 10;
double threshold_dist;

const double ratio = 1-0.0005;

char flag_term;
char term_count;

double *prev_params;

/* globals not visible to users */
unsigned int alphabet_size;
unsigned int feature_size;

unsigned int gap_size, gap_offset;
unsigned int len_size, len_offset;
unsigned int sub_size, sub_offset;
unsigned int add_size, add_offset;
static unsigned int err_size, err_offset;
static unsigned int lp_gap_size, lp_gap_offset;
static unsigned int lp_len_size, lp_len_offset;
static unsigned int lp_sub_size, lp_sub_offset;
static unsigned int lp_add_size, lp_add_offset;
static unsigned int lp_err_size, lp_err_offset;
static unsigned int lp_param_size;

static FILE *report;
extern float protein_prior[][PRIOR_SIZE];
extern float protein_joint[][MATRIX_SIZE]; 
extern short protein_matrix[][MATRIX_SIZE]; 

#define NUM_PAIRSO(a)   ((a)*((a)-1)>>1)
#define NUM_PAIRSW(a)   ((a)*((a)+1)>>1)
#define FLAT_INDEX(a,b) (((a)>(b))?(NUM_PAIRSW(a)+(b)):(NUM_PAIRSW(b)+(a)))
#define timerclear(tvp) ((tvp)->tv_sec = (tvp)->tv_usec = 0)

typedef struct tag_equality {
  int    *ind, len;
  double *val, bnd;
} type_equality;

typedef struct tag_distance { 
  unsigned int ind; 
  double dist; 
} type_distance;

/*
 * Compute the elapsed time between 'starting' and 'ending'.
 */
#include <sys/time.h> 
typedef struct timeval tv_t;
void timing(tv_t *elapsed, const tv_t *a, const tv_t *b) {
  tv_t temp;
  temp.tv_sec  = b->tv_sec -a->tv_sec;
  temp.tv_usec = b->tv_usec-a->tv_usec;
  if (temp.tv_usec < 0) { 
    temp.tv_sec--; 
    temp.tv_usec += 1000000; 
  }
  elapsed->tv_sec  += temp.tv_sec;
  elapsed->tv_usec += temp.tv_usec;
  if (elapsed->tv_usec >= 1000000) { 
    elapsed->tv_sec++; 
    elapsed->tv_usec -= 1000000; 
  }
}

/*
 * Translate the specified conventional substitution score matrix 
 * to the corresponding substitution cost matrix whose elements are in [0, 1]. 
 */
double *translate_to_cost(unsigned int matrix_index, unsigned int matrix_size)
{
  extern short protein_matrix[][MATRIX_SIZE]; 
  const short *score_matrix = protein_matrix[matrix_index];
  double *cost_matrix = (double *)malloc(matrix_size*sizeof(double));
  assert(cost_matrix);
  double min = DBL_MAX;
  double max = DBL_MIN;
  for (unsigned int i = 0; i < matrix_size; i++) {
    if (min > score_matrix[i]) min = score_matrix[i];
    if (max < score_matrix[i]) max = score_matrix[i];
  }
  for (unsigned int i = 0; i < matrix_size; i++)
    cost_matrix[i] = (max-score_matrix[i])/(max-min);
  return cost_matrix;
}

/*****************************************************************************
 *  
 *       
 *                                                    cutting plane algorithm 
 *****************************************************************************/

/* setup_linear_program
 * This function sets the number of variables, the limits of each variable,
 * and the coefficients of the objective function for the linear program.
 * LP. */
static void setup_linear_program(LPX *LP) {
  lpx_set_int_parm(LP, LPX_K_MSGLEV, 0);

  /* Set the number of variables in the linear program. */
  assert(lp_param_size > 0);
  lpx_add_cols(LP, (int)lp_param_size);

  /* Set the limits on each variables. */
  int bound_type = opt_fixd_subs ? LPX_LO : LPX_DB;
  for (unsigned int i = 1; i <= lp_param_size; i++)
    lpx_set_col_bnds(LP, i, bound_type, 0, 1);
  for (unsigned int i = lp_err_offset; i < lp_err_offset+lp_err_size; i++)
    lpx_set_col_bnds(LP, i, LPX_LO, 0, 0);  

  /* Set the coefficients of the objective function. */
  lpx_set_obj_dir(LP, LPX_MIN);
  for (unsigned int i = 1; i <= lp_param_size; i++)
    lpx_set_obj_coef(LP, i, 0);
  for (unsigned int i = lp_err_offset; i < lp_err_offset+lp_err_size; i++)
    lpx_set_obj_coef(LP, i, 1);
}

void constraint_substitution_parameters(LPX *LP, type_equality *I) {
  if (opt_fixd_subs) return; 

  /* Increase the number of inequalities in the linear program by
     the number of constrains. */
  assert(alphabet_size > 1);
  unsigned int number_constraints = ((alphabet_size-1)*alphabet_size)<<1;
  int  row = lpx_add_rows(LP, number_constraints);
  
  I->val[1] = -1;                                       /* match    coeff */
  I->val[2] = +1;                                       /* mismatch coeff */
  for (unsigned int k = 0, i = 0; i < alphabet_size; i++, k++)
    for (unsigned int j = 0; j < i; j++, k++) {
      I->ind[2] = lp_sub_offset+k;                      /* mismatch index */
      I->ind[1] = lp_sub_offset+FLAT_INDEX(i, i);       /* match    index */
      lpx_set_mat_row(LP, row, 2, I->ind, I->val);
      lpx_set_row_bnds(LP, row++, LPX_LO, 0, 0);
      I->ind[1] = lp_sub_offset+FLAT_INDEX(j, j);       /* match    index */
      lpx_set_mat_row(LP, row, 2, I->ind, I->val);	
      lpx_set_row_bnds(LP, row++, LPX_LO, 0, 0);
    }
}

void constraint_gap_parameters(LPX *LP, type_equality *I) {
  if (opt_disr_lvls <= 1) return;
  assert(opt_disr_lvls >= 2);

  int number_constraints = (opt_disr_lvls-1)<<1;
  int row = lpx_add_rows(LP, number_constraints);

  I->val[1] = -1;       /* more hydrophilic gap-open or gap-ext coefficient */
  I->val[2] = +1;       /* more hydrophobic gap-open or gap-ext coefficient */
  for (int i = 0; i < opt_disr_lvls-1; i++) {
    I->ind[1] = lp_gap_offset+i;                  /* -philic gap-open index */
    I->ind[2] = lp_gap_offset+i+1;                /* -phobic gap-open index */
    lpx_set_mat_row(LP, row, 2, I->ind, I->val);
    lpx_set_row_bnds(LP, row++, LPX_LO, 0, 0);
    I->ind[1] = lp_len_offset+i;                   /* -philic gap-ext index */
    I->ind[2] = lp_len_offset+i+1;                 /* -phobic gap-ext index */
    lpx_set_mat_row(LP, row, 2, I->ind, I->val);
    lpx_set_row_bnds(LP, row++, LPX_LO, 0, 0);
  }
}

void constraint_modification_parameters(LPX *LP, type_equality *I) {
  if (!opt_prot_pred) return;

  assert(add_size);
  int up_bound = opt_fixd_subs ? LPX_UP : LPX_DB;
  int lo_bound = opt_fixd_subs ? LPX_LO : LPX_DB; 
  lpx_set_col_bnds(LP, (int)lp_add_offset+0, up_bound, -1, 0);  /* Ah-Ah */
  lpx_set_col_bnds(LP, (int)lp_add_offset+1, lo_bound,  0, 1);  /* Ah-Bs */
  lpx_set_col_bnds(LP, (int)lp_add_offset+2, up_bound, -1, 0);  /* Bs-Bs */
  lpx_set_col_bnds(LP, (int)lp_add_offset+3, lo_bound,  0, 1);  /* Ah-Lp */
  lpx_set_col_bnds(LP, (int)lp_add_offset+4, lo_bound,  0, 1);  /* Bs-Lp */
  lpx_set_col_bnds(LP, (int)lp_add_offset+5, up_bound, -1, 0);  /* Lp-Lp */

	if (opt_fixd_subs) return;

  int row = lpx_add_rows(LP, sub_size*add_size);
  I->val[1] = 1;  /* substitution parameter */
  I->val[2] = 1;  /* modification parameter */
  for (int i = 0; i < sub_size; i++) {
    I->ind[1] = lp_sub_offset+i;
    for (int j = 0; j < add_size; j++) {
      I->ind[2] = lp_add_offset+j;
      lpx_set_mat_row(LP, row, 2, I->ind, I->val);
      lpx_set_row_bnds(LP, row++, LPX_LO, 0, 0);
    }
  }
}

void constraint_on_degeneracy
(LPX *LP, type_equality *equal, unsigned int matrix_index, const char *used_params) {
  if (opt_fixd_subs) return;
  /* Compute sums of priors and the weighted costs for each of identity and 
     substitution types where the costs are translated from a standard 
     substitution score matrix. */
  extern float protein_joint[][MATRIX_SIZE]; 
  const float *joint_dist = protein_joint[0];
  double *std_matrix = translate_to_cost(matrix_index, sub_size);
  const char *used_subst = &used_params[sub_offset];
  double subst_dist, subst_cost;
  double ident_dist, ident_cost;
  subst_dist = subst_cost = 0;
  ident_dist = ident_cost = 0;
  for (unsigned int i = 0, k = 0; i < alphabet_size; i++, k++) {
    for (unsigned int j = 0; j < i; j++, k++)
      if (used_subst[k]) {
	subst_dist += joint_dist[k];
	subst_cost += joint_dist[k]*std_matrix[k];
      }
    if (used_subst[k]) {
      ident_dist += joint_dist[k];
      ident_cost += joint_dist[k]*std_matrix[k];
    }
  }
  free(std_matrix);

  /* Compute the expected costs for each of used (used means the feature is 
     present in some input example) identity and used substitution types. */
  double degen_threshold = subst_cost/subst_dist-ident_cost/ident_dist;
  if (degen_threshold <= 1e-6)
    error("Threshold in nondegeneracy inequality is virtually nonpositive \
(%.6f).", degen_threshold);

  /* Add the nondegeneracy inequality to the linear program. */
  equal->len = 1;  /* GLPK indexing scheme, which starts at 1 not 0. */
  for (unsigned int i = 0, k = 0; i < alphabet_size; i++, k++) {
    for (unsigned int j = 0; j < i; j++, k++) 
      if (used_subst[k]*joint_dist[k] && subst_dist) {
	equal->val[equal->len]   = used_subst[k]*joint_dist[k]/subst_dist;
	equal->ind[equal->len++] = lp_sub_offset+k;
      }
    if (used_subst[k]*joint_dist[k] && ident_dist) {
      equal->val[equal->len]   = -used_subst[k]*joint_dist[k]/ident_dist;
      equal->ind[equal->len++] = lp_sub_offset+k;
    }
  }
  equal->len--;
  unsigned int row = lpx_add_rows(LP, 1);
  lpx_set_mat_row (LP, row, equal->len, equal->ind, equal->val);
  lpx_set_row_bnds(LP, row, LPX_LO, degen_threshold, 0);
  lpx_set_row_name(LP, row, "nondegeneracy");
}


void make_equality
(const double norm_const, const double *B, const double *A, double recv_rate, unsigned int id, type_equality *I, const double *P, char *trained_params) {
  double subst_const = 0;
  int i = 1;
  for (int j = 0; j < feature_size; j++) 
    if (A[j] != B[j]) {
      trained_params[j] = 1;
      if (opt_fixd_subs) {
	if (sub_offset <= j && j < sub_offset+sub_size) 
	  subst_const += (A[j]-B[j])*P[j];
	else {
	  I->ind[i] = j+1-(j >= sub_offset ? sub_size : 0);
	  I->val[i] = B[j]-A[j];
	  i++;
	}
      }
      else {
	I->ind[i] = j+1;
	I->val[i] = B[j]-A[j];
	i++;
      }
    }
  I->ind[i] = lp_err_offset+id;
  I->val[i] = norm_const/opt_conv_wght;
  I->len = i;
  I->bnd = norm_const*(1/opt_conv_wght-1)*(1-recv_rate)+subst_const;
}

/* With this function, qsort will sort things in increasing order. */
int compare_distance(const void *x, const void *y) {
  return ((type_distance *)x)->dist > ((type_distance *)y)->dist ? 1 
    : ((type_distance *)x)->dist < ((type_distance *)y)->dist ? -1 : 0; 
}

double projected_distance
(const double *V, const type_equality *I, const double *P) {
  /* Measure the magnitude of V. */
  double m = 0;
  for (unsigned int i = 0; i < lp_param_size; i++)
    m += V[i]*V[i];
  /* Get the dot product x of V and N, where N is normal to plane I.
     Also get the dot product y of P and N, where P is a position vector. */
  double x = 0;
  double y = -I->bnd;
  for (unsigned int i = 1; i <= (unsigned int)I->len; i++) {
    unsigned int ind = I->ind[i];
    x += V[ind-1]*I->val[i];
    if (opt_fixd_subs && ind >= lp_sub_offset) 
      ind += sub_size;
    y += P[ind-1]*I->val[i];
  }
  /* If dot product of V and N is zero, then V and I is parallel with the
     infinite projected distance. */
  if (x == 0) return DBL_MAX;
  return y/x*sqrt(m);
}

double shortest_distance(const type_equality *I, const double *P) {
  /* This is a speicial version of projected_distance function, special
     in a sense that the projection vector is the normal of the plane. */
  if (!I->len) return DBL_MAX;
  double x = 0;
  double y = -I->bnd;
  for (unsigned int i = 1; i <= (unsigned int)I->len; i++) {
    unsigned int ind = I->ind[i];
    if (opt_fixd_subs && ind >= lp_sub_offset) 
      ind += sub_size;
    x += I->val[i]*I->val[i];
    y += P[ind-1]*I->val[i];
  }
  return y/sqrt(x);
}

static void cull_inequalities
(LPX *LP, unsigned int base, unsigned int *count, const double *P, type_equality *I) {
  /* Choosing constraints:  This is a heuristic attempt to reduce the number of
   * constraints in the current linear program for speeding up the computation
   * time by GLPK.  For this, we choose a few constraints to keep and remove 
   * all the others in the current linear program.  Which ones to be kept?
   * That should be the art of this.  Here we just simply keep those close to 
   * the current point P. */

  /* Retrive the objective vector vec. */
  double *obj_vector = (double *)malloc(lp_param_size*sizeof(double));
  assert(obj_vector);
  for (unsigned int i = 0; i < lp_param_size; i++) 
    obj_vector[i] = lpx_get_obj_coef(LP, i+1);
  
  /* Measure the distance between the point P and the halfspace plane I. */
  unsigned int row = lpx_get_num_rows(LP);
  unsigned int number_constraints = row-base+1;
  type_distance *list = (type_distance *)malloc
    (number_constraints*sizeof(type_distance));
  assert(list);
  for (unsigned int i = base; i <= row; i++) {
    I->len = lpx_get_mat_row(LP, i, I->ind, I->val);
    assert(I->len);
    I->bnd = lpx_get_row_lb (LP, i);
    list[i-base].ind = i;
    //list[i-base].dist = shortest_distance(I, P);
    list[i-base].dist = projected_distance(obj_vector, I, P);
  }  
  free(obj_vector);

  /* Sort all constraints in an increasing order of distance from P. */
  qsort(list, number_constraints, sizeof(type_distance), compare_distance);
  /* Find the first postive distance. */
  unsigned int end = 0;
  for (; end < number_constraints; end++)
    if (list[end].dist > 10e-4) break;
  assert(*count >= lp_param_size);
  if (*count < end) *count = end;
  else if (*count > lp_param_size) *count = lp_param_size;

  if (*count == number_constraints) return;
  assert(*count >= end);

  /* Delete the constraints far from P, which should be at the end of list. */
  unsigned int number_deleted = number_constraints-(*count);
  int *deletion_list = (int *)malloc((number_deleted+1)*sizeof(int));
  assert(deletion_list);
  for (unsigned int i = 1; i <= number_deleted; i++)
    deletion_list[i] = list[*count+i-1].ind;
  lpx_del_rows(LP, number_deleted, deletion_list);
  free(deletion_list);
  free(list);
  lpx_std_basis(LP);
  if (opt_debg >= 2) {
    fprintf(stderr, "(%d)", opt_debg);
    fprintf(stderr, "Constraints: %d==>%d\n", number_constraints, *count);
  }
}

static inline void update_parameters(LPX *LP, double *P) {
  const unsigned int example_count = err_size;
  const unsigned int feature_count = feature_size;
  memcpy(prev_params, P, (feature_count+example_count)*sizeof(double));
  
  if (opt_fixd_subs) {
    assert(!lp_sub_size);
    unsigned int i = 1;
    for (unsigned int j = 0; i < lp_add_offset; i++, j++) {
      P[gap_offset+j] = lpx_get_col_prim(LP, (int)i);
      if (P[gap_offset+j] < 0) P[gap_offset+j] = 0;
    }
    for (unsigned int j = 0; i < lp_err_offset; i++, j++)
      P[add_offset+j] = lpx_get_col_prim(LP, (int)i);
    for (unsigned int j = 0; i <= lp_param_size; i++, j++) {
      P[err_offset+j] = lpx_get_col_prim(LP, (int)i);
      if (P[err_offset+j] < 0) P[err_offset+j] = 0;
    }
  }
  else {  
    for (unsigned int i = 1; i <= lp_param_size; i++) {
      P[i-1] = lpx_get_col_prim(LP, (int)i);
      if (P[i-1] < 0) P[i-1] = 0;
      if (i < lp_err_offset && P[i-1] > 1) P[i-1] = 1;
    }
    for (unsigned int i = lp_add_offset; i < lp_add_offset+lp_add_size; i++) {
      P[i-1] = lpx_get_col_prim(LP, (int)i);
      if (P[i-1] > 1)  P[i-1] = 1;
      if (P[i-1] < -1) P[i-1] = -1;
    }
  }

  double distance = 0;
  for (int i = 0; i < feature_count; i++)
    distance += (P[i]-prev_params[i])*(P[i]-prev_params[i]);
  distance = sqrt(distance/feature_count);
  if (distance > 1e-4)
    term_count = 0;
  else 
    term_count++;
  if (term_count == 5)
    flag_term = 1;

  if (opt_debg == 0) return;
  const double *G = &P[gap_offset]; double g = 0;
  const double *L = &P[len_offset]; double l = 0;
  const double *S = &P[sub_offset]; double s = 0;
  const double *M = &P[add_offset]; double m = 0;
  const double *E = &P[err_offset]; double e = 0;
  for (unsigned int i = 0; i < gap_size; i++) g += G[i];
  for (unsigned int i = 0; i < len_size; i++) l += L[i];
  for (unsigned int i = 0; i < sub_size; i++) s += S[i];
  for (unsigned int i = 0; i < add_size; i++) m += M[i];
  for (unsigned int i = 0; i < err_size; i++) e += E[i];
  fprintf(stderr, "params: %.2e(o) ", lpx_get_obj_val(LP)/example_count);
  fprintf(stderr, "%.2e(g) %.2e(l) ", g/gap_size, l/len_size);

  fprintf(stdout, "\n");
  for (int i = 0; i < gap_size; i++)
    fprintf(stdout, "%f ", G[i]);
  fprintf(stdout, "\n");



  if (!opt_fixd_subs) fprintf(stderr, "%.2e(s) ", s/sub_size);
  if (opt_prot_pred)  fprintf(stderr, "%+.2e(m) ", m/add_size);
  fprintf(stderr, "%.2e(e) ", e/err_size);
  fprintf(stderr, "%.2e(i)\n", distance);
}

static char *obtain_used_params
(double **feature, unsigned int example_count, unsigned int feature_count) 
{
  char *used_params = (char *)calloc(feature_count, sizeof(char));
  assert(used_params);
  /* Find out which features are present in the feature vectors and 
     record the information in a boolean vector and return it. */
  for (unsigned int i = 0; i < feature_count; i++)
    for (unsigned int j = 0; j < example_count; j++)
      if (feature[j][i]) { used_params[i] = 1; break; }
  return used_params;
}

static void handle_untrained_params
(double *params, const char *used_params, const char *trained_params, 
 unsigned int matrix_index)
{
  /* Find out if there exists an untrained substitution parameter or not. */
  const char *trained_subst = &trained_params[sub_offset];
  char untrained_exist = 0;
  for (unsigned int i = 0; i < sub_size; i++)
    if (!trained_subst[i]) { untrained_exist = 1; break; }
  if (!untrained_exist) return;

  /* Find out if there exists an used but not trained substitution parameter 
     or not. */
  const char *used_subst = &used_params[sub_offset];
  char intersection = 0;
  for (unsigned int i = 0; i < sub_size; i++)
    if (!trained_subst[i] && used_subst[i]) { intersection = 1; break; }

  /* */
  if (intersection) {
    /* */
    extern float protein_joint[][MATRIX_SIZE]; 
    const float *joint_dist = protein_joint[0];
    double *subst = &params[sub_offset];
    double *std_matrix = translate_to_cost(matrix_index, sub_size);
    double subst_dist, subst_comp, subst_cost;
    double ident_dist, ident_comp, ident_cost;
    subst_dist = subst_comp = subst_cost = 0;
    ident_dist = ident_comp = ident_cost = 0;
    double z = DBL_MIN;
    unsigned int k = 0;
    for (unsigned int i = 0; i < alphabet_size; i++, k++) {
      for (unsigned int j = 0; j < i; j++, k++)
	if (used_subst[k]) {
	  subst_dist += joint_dist[k];
	  subst_cost += joint_dist[k]*std_matrix[k];
	  if (trained_subst[k]) subst_cost -= joint_dist[k]*subst[k];
	  else                  subst_comp += joint_dist[k];
	}
      if (used_subst[k]) {
	ident_dist += joint_dist[k];
	ident_cost += joint_dist[k]*std_matrix[k];
	if (trained_subst[k]) ident_cost -= joint_dist[k]*subst[k];
	else                  ident_comp += joint_dist[k];
	if (trained_subst[k] && z < subst[k]) z = subst[k];
      }
    }
    free(std_matrix);
    /* */
    double a = subst_comp/subst_dist;
    double b = ident_comp/ident_dist;
    double c = subst_cost/subst_dist-ident_cost/ident_dist;
    double new_subst, new_ident;
    new_subst = new_ident = 0;
    if (a == 0 && b != 0 && c >= 0) new_ident = -c/b < 1 ? -c/b : 1;
    else if (a != 0 && b == 0 && c <= a) new_subst = 1;
    else if (a != 0 && b != 0 && c <= a) 
      { new_subst = 1; new_ident = (a-c)/b < z ? (a-c)/b : z; }
    else error("Something is wrong. a=%f b=%f c=%f", a, b, c);
    /* */
    k = 0;
    for (unsigned int i = 0; i < alphabet_size; i++, k++) {
      for (unsigned int j = 0; j < i; j++, k++)
	if (used_subst[k] && !trained_subst[k]) subst[k] = new_subst;
      if (used_subst[k] && !trained_subst[k]) subst[k] = new_ident;
    }
    fprintf(report, "subst=%f ident=%f\n", new_subst, new_ident);
  }
}

static int 
find_error_parameters
(
 const alignment_t *S, 
 const alignment_t *confidences, 
 const alignment_t *predictions, 
 const dist_t ***distributions, 
 const match_t *M, 
 const double **F, 
 type_table *T, 
 double **H, 
 double *P, 
 double *rate, 
 const char *used_params
 ) {
  const unsigned int example_count = err_size;
  const unsigned int feature_count = feature_size;

  type_equality I;
  I.ind =    (int *)malloc((1+feature_count+1)*sizeof(int));
  assert(I.ind);
  I.val = (double *)malloc((1+feature_count+1)*sizeof(double));
  assert(I.val);
  
  double **V = (double **)malloc(num_bests*sizeof(double *));
  assert(V);
  for (int i = 0; i < num_bests; i++) {
    V[i] = (double *)malloc(feature_count*sizeof(double));
    assert(V[i]);
  }
  
  char *trained_params = (char *)calloc(feature_count, sizeof(char));
  assert(trained_params);

  prev_params = (double *)malloc((feature_count+example_count)*sizeof(double));
  assert(prev_params != NULL);

  tv_t st, et, time_lp, time_so;
  timerclear(&time_lp);     /* the time spent by GLPK solver */
  timerclear(&time_so);     /* the time spent by the separation oracle */

  double *delta = &P[err_offset];

  unsigned int num_added   = 0; /* the number of constraints added */
  unsigned int num_reduced = 0; /* the number of the reduction was done */
  unsigned int num_solved  = 0; /* the number of the GLPK solver was called */
  unsigned int num_called  = 0; /* the number of the separation oracle was called */

  /***************************************************************************/
  /*                                                                         */
  /*                                the cutting plane algorithm starts here. */ 
  /***************************************************************************/
  LPX *LP = lpx_create_prob();
  setup_linear_program(LP);
  constraint_substitution_parameters(LP, &I);
  constraint_modification_parameters(LP, &I);
  constraint_gap_parameters(LP, &I);
  constraint_on_degeneracy(LP, &I, MATRIX_BLOSUM62_3, used_params);
  int base = lpx_get_num_rows(LP)+1;

  /* Seeding the linear program:  This is a heuristic attempt for speeding up
   * finding the optimal parameter solution.  We experinced that when the 
   * empty linear program starts with constraints each of which is obtained by
   * the example and an optimal alignment under the optimal (or default) 
   * parameter choice from the previous iteration (or for the first iteration).
   * We do not care if these constraints are violated or not under this
   * parameter choice.  We add them to the linear program and solve it. */

  /*
  for (unsigned int i = 0; i < feature_size; i++)
    fprintf(stderr, "%f ", P[i]);
  fprintf(stderr, "\n");
  */

  //if(0) {
  for (unsigned int i = 0; i < example_count; i++) {
    /* Ask the separation oracle for violated inequalities. */
    gettimeofday(&st, 0);
    unsigned int num_violated = (*work_the_oracle)
      (T, H, S, confidences, predictions, distributions, M, F[i], i, 0, 0,
       S[M[i].a].length[M[i].s]-1, S[M[i].a].length[M[i].t]-1,
       0, P, V, 1, opt_term_gaps, opt_term_gaps, &rate[i]);
    gettimeofday(&et, 0);
    timing(&time_so, &st, &et); 
    num_called++;
    
    /* The following condition meets if and only if the example score is the
     * same as the optimal score and the recovery rate measure on this example 
     * is perfect.  Otherwise, the condition must fail. */
    if (!num_violated) continue;
    
    /* This example is not satisfied.  Add the violated inequalities. */
    num_added += num_violated;
    unsigned int row = lpx_add_rows(LP, num_violated);
    for (unsigned int j = 0; j < num_violated; j++, row++) {
      double norm_const = 1;
      if (opt_norm_cost == 1)
	norm_const = (S[M[i].a].length[M[i].s]+S[M[i].a].length[M[i].t])/2.0;
      else if (opt_norm_cost == 2)
	norm_const = MAX(S[M[i].a].length[M[i].s], S[M[i].a].length[M[i].t]);
      assert(norm_const > 0);
      make_equality(norm_const, V[j], F[i], rate[i], i, &I, P, trained_params);
      lpx_set_mat_row(LP, row, I.len, I.ind, I.val);
      lpx_set_row_bnds(LP, row, LPX_LO, I.bnd, 0);
    }
    /*
  int status = lpx_get_status(LP);
  if (status == LPX_OPT) update_parameters(LP, P);
  else return status;
  /*
  else {fprintf(stderr, "%d Error!\n", i); exit(1); }
  fprintf(stderr, "Solved.\n");
  */
  }
  //}
  //const unsigned int example_count = err_size;
  //const unsigned int feature_count = feature_size;  
  if (!opt_fixd_subs)
    memset(P, 0, (feature_count+example_count)*sizeof(double));

  /* start the timer */ gettimeofday(&st, 0);
  lpx_simplex(LP);
  /* stop  the timer */ gettimeofday(&et, 0);
  timing(&time_lp, &st, &et);
  num_solved++;


  flag_term = 0;
  term_count = 0;

  int status = lpx_get_status(LP);
  if (status == LPX_OPT) 
    update_parameters(LP, P);
    
  /* Example management:  This is a heuristic algorithm to effectively manage 
   * examples.  There are two lists of examples.  One is called the "current" 
   * list which contains examples under consideration and the other called the 
   * "other" list.  Once an example in the current list is not interesting any 
   * more, it is moved to the other list.  And the current list becomes empty,
   * we switch the lists.  If all examples in the current list are satisfied
   * with the current choice of parameters, then we are done. */

  /* Initialize variables related to the example lists.  Initialize the other
   * list which becomes the current list right away. */
  unsigned int *list[2], num_members[2] = { example_count, 0 };
  list[0] = (unsigned int *)malloc(example_count*sizeof(unsigned int));
  assert(list[0]);
  list[1] = (unsigned int *)malloc(example_count*sizeof(unsigned int));
  assert(list[1]);
  for (int i = 0; i < (int)example_count; i++) 
    list[0][i] = i;
  
  unsigned int curr = 1;
  char done = 0;
  unsigned int cull_size = lp_param_size;
  while (!done && status == LPX_OPT && flag_term == 0) {
    done = 1;
    /* Switch to the other list and let's call it current.  Then the current  
     * list has all examples to be examined while the other list doesn't have 
     * any. */
    curr = (curr+1)%2;

    /* While there are examples to be examined, do the following. */
    unsigned int num_unsatisfied = 0;
    unsigned int first_id = (unsigned int)-1;
    while (num_members[curr] && status == LPX_OPT && flag_term == 0) { 
      /* For each example in the current list, ask the separation oracle to
       * give (hopely good) violated inequalities if there's any.  If 
       * the oracle does not give back any, it means that the example is 
       * satisfied under the current parameter choice.  Otherwise, the example 
       * is not satisfied and a new parameter choice should be sought. */
      for (unsigned int i = 0; i < example_count && status == LPX_OPT && 
	     flag_term == 0; i++) {
	/* If this is not a valid example, skip and move on. */
	unsigned int id = list[curr][i];
	if (id == (unsigned int)-1) continue;
	
	/* Solve the linear program either if enough additional constraints
	 * are collected or if there's not enough unsatisfied examples left
	 * contributing any constraints to meet this "enoughness". */
	if ((num_unsatisfied >= threshold_e) || (first_id == id)) {
	  gettimeofday(&st, 0);
	  //fprintf(stdout, "Solving... "); fflush(stdout);
	  lpx_simplex(LP);
	  //fprintf(stdout, "Solved...\n"); fflush(stdout);
	  gettimeofday(&et, 0);
	  timing(&time_lp, &st, &et);
	  num_solved++;
	  num_unsatisfied = 0;
	  first_id = (unsigned int)-1;
	  status = lpx_get_status(LP);
	  /* If the problem is not feasible, break this loop to quit. */
	  if (LPX_OPT != status) continue;
	  /* Update the alignment and error parameters. */
	  update_parameters(LP, P);

	  /* If too many constraints are present in the linear program, 
	   * let's reduce the number of constraints in some way. */
	  if ((unsigned int)(lpx_get_num_rows(LP)-base+1) > cull_size*opt_cull_fold) {
	    cull_inequalities(LP, base, &cull_size, P, &I);
	    num_reduced++;
	  }
	}
	
	/* Ask the separation oracle for violated inequalities. */
	/* start the timer */ gettimeofday(&st, 0);
	unsigned int num_violated = (*work_the_oracle)
	  (T, H, S, confidences, predictions, distributions, M, F[id], id, 0, 
	   0, S[M[id].a].length[M[id].s]-1, S[M[id].a].length[M[id].t]-1, 
	   delta[id], P, V, 1, opt_term_gaps, opt_term_gaps, &rate[id]);
	/* stop  the timer */ gettimeofday(&et, 0);
	timing(&time_so, &st, &et); 
	num_called++;
	
	/* This example is satisfied.  Let's move on to the next example. */
	if (!num_violated) { 
	  /* Move the example to the other list. */
	  list[curr][i] = (unsigned int)-1; 
	  num_members[curr]--;
	  unsigned int other = (curr+1)%2;
	  list[other][num_members[other]] = id;
	  num_members[other]++;
	  continue;
	}
	
	/* This example is not satisfied.  Add the violated inequalities. */
	done = 0;
	num_added += num_violated;
	num_unsatisfied += num_violated;
	if (first_id == (unsigned int)-1) first_id = id;
	/*
	fprintf(stdout, "%2d(%2d)", id, num_unsatisfied);
	fflush(stdout);
	*/
	/* Add the constraints to the linear program. */
	unsigned int row = lpx_add_rows(LP, num_violated);
	for (unsigned int k = 0; k < num_violated; k++, row++) {
	  double norm_const = 1;
	  if (opt_norm_cost == 1)
	    norm_const = (S[M[i].a].length[M[i].s]+S[M[i].a].length[M[i].t])/2.0;
	  else if (opt_norm_cost == 2)
	    norm_const = MAX(S[M[i].a].length[M[i].s], S[M[i].a].length[M[i].t]);
	  assert(norm_const > 0);
	  make_equality(norm_const, V[k], F[id], rate[id], id, &I, P, trained_params);
	  lpx_set_mat_row(LP, row, I.len, I.ind, I.val);
	  lpx_set_row_bnds(LP, row, LPX_LO, I.bnd, 0);
	}
      }
    }
    /* At this point, all examples are satisfied at least once at a certain
     * moment.  The current list should be empty and the other list full. */
  }
  
  /***************************************************************************/
  /*                                                                         */
  /*                                  The cutting plane algorithm ends here. */
  /***************************************************************************/

  mean_score = mean_recv = 0;
  if (LPX_OPT == status) {
    /* Try and make all parameters reasonable. */
    /*
    handle_untrained_params
      (P, used_params, trained_params, MATRIX_BLOSUM62_3);
    */
    /* Measure the recovery rate of an optimal scoring alignment under the
     * optimal parameter choice with respect to its corresponding example. */
    for (unsigned int i = 0; i < example_count; i++) {
      (*work_the_oracle)
	(T, H, S, confidences, predictions, distributions,
	 M, 0, i, 0, 0, S[M[i].a].length[M[i].s]-1, 
	 S[M[i].a].length[M[i].t]-1, 0, P, V, 0, opt_term_gaps, opt_term_gaps, &rate[i]);
      double x = 0;
      double this_mean = (S[M[i].a].length[M[i].s]+S[M[i].a].length[M[i].t])/2.0;
      for (unsigned int j = 0; j < feature_count; j++)
	x += (F[i]-V[j])*P[j];
      mean_score += x/this_mean;
      mean_recv += rate[i];
    }
    mean_score /= example_count;
    mean_recv /= example_count;
  }

  /* Report some statistics on the GLPK solver and the separation oracle. */
  fprintf(report, "Constraints added : %d\n", num_added);
  fprintf(report, "Reduction executed: %d\n", num_reduced);
  fprintf(report, "GLPK solver called: %d\n", num_solved); 
  fprintf(report, "Separation called : %d\n", num_called);
  fprintf(report, "GLPK solver spent : %d.%02d\n", 
	  (int)time_lp.tv_sec, (int)(time_lp.tv_usec/10000.0+0.5));
  fprintf(report, "Separation spent  : %d.%02d\n", 
	  (int)time_so.tv_sec, (int)(time_so.tv_usec/10000.0+0.5));

  /* Clean up the mess, I mean, the allocated memories. */
  free(prev_params);
  free(trained_params);
  lpx_delete_prob(LP);
  free_features(&V, num_bests);
  free(I.val);
  free(I.ind);
  free(list[0]);
  free(list[1]);
  return status;
}

/*****************************************************************************
 *
 *
 *                                                    iterations, iterations.
 *****************************************************************************/
static inline void write_parameters
(const char *path, double error, double rate, double *P) {
  FILE *stream = fopen(path, "w");
  assert(stream);
  fprintf(stream, "%d\n", opt_disr_lvls);
  fprintf(stream, "%d\n", opt_term_gaps);
  fprintf(stream, "%d\n", (int)alphabet_size);
  unsigned int i;
  for (i = 0; i < feature_size; i++) 
    fprintf(stream, "%.40e\n", P[i]);
  fprintf(stream, "%.40e\n", rate);
  fprintf(stream, "%.40e\n", error);
  fclose(stream);
}

static inline void read_parameters
(const char *path, double *error, double *rate, double *P) {
  FILE *stream = fopen(path, "r");
  assert(stream);
  int x;
  fscanf(stream, "%d\n", &x); opt_disr_lvls = x;
  fscanf(stream, "%d\n", &x); opt_term_gaps = x;
  fscanf(stream, "%d\n", &x); alphabet_size = x;
  for (unsigned int i = 0; i < feature_size; i++)
    fscanf(stream, "%lf\n", &P[i]);
  if (rate)  fscanf(stream, "%lf\n", rate);
  if (error) fscanf(stream, "%lf\n", error);
  fclose(stream);
}

static int find_parameters
(const alignment_t *S, const alignment_t *confidences, const alignment_t *predictions, const dist_t ***distributions, unsigned int num_alignments, const match_t *M, const block_t *B, const double **W, const char *output) { 
  /* make up the output file names. */
  char param_file[2048+1];
  sprintf(param_file, "param.%c", opt_fixd_subs ? 'f' : 'v');
  if (output) { strcat(param_file, "."); strcat(param_file, output); }
  char report_file[2048+1];
  sprintf(report_file, "score.%c", opt_fixd_subs ? 'f' : 'v');
  if (output) { strcat(report_file, "."); strcat(report_file, output); }  
  /* allocate heap for parameters. */
  double *P = (double *)malloc((feature_size+err_size)*sizeof(double));
  assert(P);
  double *rate = (double *)malloc(err_size*sizeof(double));
  assert(rate);
  unsigned int i, j;
  double **F = (double **)malloc(err_size*sizeof(double *));
  assert(F);
  for (i = 0; i < err_size; i++) {
    F[i] = (double *)malloc(feature_size*sizeof(double));
    assert(F[i]);
    for (j = 0; j < feature_size; j++)
      F[i][j] = DBL_MAX;
  }
  /* prepare separation procedure. */ 

  double temp = 0;
  int temp_count = 0;
  unsigned int k = 0;
  for (i = 0; i < num_alignments; i++) {
    temp_count += S[i].count; //
    for (j = 0; j < S[i].count; j++) {
      temp += S[i].length[j]; //
      if (S[i].length[j] > k) 
	k = S[i].length[j];
    }
  }
  temp /= temp_count; //
  //printf("%f\n", temp); //
  type_table D;
  create_table(&D, k, k);
  double *H[2];
  for (i = 0; i < 2; i++) {
    H[i] = (double *)malloc(k*sizeof(double));
    assert(H[i]);
  }

  const unsigned int example_count = err_size;
  const unsigned int feature_count = feature_size;  

  /* Seeding parameter values:  Initialize the alignment parameters with 
     the default values (currently, Opal default parameters are in use) 
     and the error parameters with zeros. */
  memset(P, 0, (feature_count+example_count)*sizeof(double));

  /* Default values for parameters. */
  double *PP = &P[sub_offset];
  switch (opt_alph_type) {
  case 0:
    for (int i = 0; i < sub_size; i++)
      PP[i] = (4.0*(16-protein_matrix[MATRIX_BLOSUM62_3][i]))/88;
    break;
  case 1: case 2:
    PP[0] = 40; 
    PP[1] = 68; PP[2] = 12;
    PP[3] = 64; PP[4] = 80; PP[5] = 32;
    PP[6] = 64; PP[7] = 68; PP[8] = 72; PP[9] = 36;
    for (int i = 0; i < sub_size; i++)
      PP[i] /= 80.0;
    break;
  }
  P[gap_offset+gap_size-1] = 15.0/88;
  P[len_offset+len_size-1] = 36.0/88;
  for (int i = 0; i < gap_size-opt_term_gaps; i++) 
    P[gap_offset+i] = 60.0/88;
  for (int i = 0; i < len_size-opt_term_gaps; i++) 
    P[len_offset+i] = 38.0/88;

/*
  if (sw.p) {
    P[add_offset+0] = -16.0/88; //-0.14;
    P[add_offset+2] = -34.0/88; //-0.28;
    P[add_offset+3] =   5.0/88; //+0.04;
  }
*/

  /* search the smallest average absolute error. */
  int status;
  bool found = false;
  double prev_error, curr_error = DBL_MAX;
  int num_iterations = 0;
  int best_iteration = 1;
  const double *delta = &P[err_offset];
  double best_cover = 0; 
  report = fopen(report_file, "w");
  assert(report);
  fprintf(report, "(Start)\n");
  fprintf(report, "Trivial completion: %d\n", opt_init_triv);
  fprintf(report, "Convextity weight: %.4e\n", opt_conv_wght);
  fprintf(report, "Terminal: %d\n", opt_term_gaps);
  fprintf(report, "Number of hydrophobicity level: %d\n", opt_disr_lvls);
  fprintf(report, "Window size: %d\n", opt_wind_size);
  fprintf(report, "Violation threshold: %f\n", threshold_v);
  if (opt_verb) {
    fprintf(stderr, "Trivial completion: %d\n", opt_init_triv);
    fprintf(stderr, "Convextity weight: %.4e\n", opt_conv_wght);
    fprintf(stderr, "Terminal gap cost: %d\n", opt_term_gaps);
    fprintf(stderr, "Number of disruption level: %d\n", opt_disr_lvls);
    fprintf(stderr, "Window size: %d\n", opt_wind_size);
    fprintf(stderr, "Violation threshold: %f\n", threshold_v);
  }

  char *used_params = 0;
  
  do {
    prev_error = curr_error;
    num_iterations++;
    fprintf(report, "* Iteration (%d)\n", (int)num_iterations);
    if (opt_verb) fprintf(stderr, "* Iteration (%d)\n", (int)num_iterations);
    if (num_iterations == 1 && opt_init_triv) 
      for (i = 0; i < err_size; i++) 
	memcpy(F[i], W[i], feature_size*sizeof(double));
    else 
      complete_examples(&D, H, S, confidences, predictions, distributions,
			M, B, err_size, P, F);
    if (num_iterations == 1)
      used_params = obtain_used_params(F, err_size, feature_count);
   
    tv_t st, et, rt;
    timerclear(&rt);
    gettimeofday(&st, 0);
    status = find_error_parameters
      (S, confidences, predictions, distributions, 
       M, (const double **)F, &D, H, P, rate, (const char *)used_params);
    /* modify here. 
    opt_recv_weight = 1*(mean_score/mean_recv);
    status = find_error_parameters
      (S, confidences, predictions, distributions, 
       M, (const double **)F, &D, H, P, rate, (const char *)used_params);
    */
    gettimeofday(&et, 0);
    timing(&rt, &st, &et);
#define LPX_STATUS(s)							\
  (LPX_OPT==(s)    ? "optimal"     : LPX_FEAS==(s)   ? "feasible"   :	\
   LPX_NOFEAS==(s) ? "no feasible" : LPX_INFEAS==(s) ? "infeasible" :	\
   LPX_UNBND==(s)  ? "unbounded"   : "undefined")  
    fprintf(report, "Linear program    : %s\n", LPX_STATUS(status));
    if (opt_verb) fprintf(stderr, "Linear program    : %s\n", LPX_STATUS(status));
#undef LPX_STATUS
    if (status != LPX_OPT) break;
    found = true;
    double cover = 0;
    curr_error = 0;
    for (i = 0; i < err_size; i++) 
      { cover += rate[i]; curr_error += delta[i]; }
    curr_error /= err_size;
    cover /= err_size/100.0;
    /* save the best recovery rate so far. */
    if (best_cover < cover) { 
      best_cover = cover;
      best_iteration = num_iterations;
    }
    if (opt_verb) fprintf(stderr, "Average error     : %.4e\n", curr_error);
    if (opt_verb) fprintf(stderr, "Average recovery  : %.2f\n", cover);
    char buffer[2048+1];
    sprintf(buffer, "%s.%02d", param_file, (int)num_iterations);
    write_parameters(buffer, curr_error, cover, P);
    fprintf(report, "Time spent        : %d.%02d\n", 
	    (int)rt.tv_sec, (int)(rt.tv_usec/10000.0+0.5));
    fprintf(report, "Average error     : %.4e\n", curr_error);
    fprintf(report, "Average recovery  : %.2f\n", cover);
    fflush(report);
    if (cover == 100) break;
  } while (curr_error/prev_error <= ratio && 
	   num_iterations < opt_maxm_itrs);

  fprintf(report, "(End)\n");
  fclose(report);
  free(used_params);
  for (i = 0; i < 2; i++) 
    free(H[i]);
  delete_table(&D);
  free_features(&F, err_size);
  free(rate);
  free(P);
  if (found == true) {
    sprintf(report_file, "\\cp -f %s.%02d %s", param_file, 
	    (int)best_iteration, param_file);
    system(report_file);
  }
  return (found == true) ? LPX_OPT : status;
}

 
static int apply_parameters
(const alignment_t *S, const alignment_t *confidences, const alignment_t *predictions, const dist_t ***distributions, unsigned int num_alignments, const match_t *M, const block_t *B, const double **W, const char *output) { 
  if (!opt_just_algn) return 0;

  const unsigned int num_examples = err_size;
  const unsigned int num_features = feature_size;
  
  /* I need to write code parsing the parameter file. */

  /* Allocate heap for parameters. */
  double *P = (double *)malloc((num_features+num_examples)*sizeof(double));
  assert(P);
  double *F = (double *)malloc(num_features*sizeof(double));
  assert(F);
  double *rate = (double *)malloc(num_examples*sizeof(double));
  assert(rate);
  double *cost = (double *)malloc(num_examples*sizeof(double));
  assert(cost);
  
  /* Prepare separation procedure. */ 
  unsigned int k = 0;
  for (unsigned int i = 0; i < num_alignments; i++)
    for (unsigned int j = 0; j < S[i].count; j++) 
      if (S[i].length[j] > k) k = S[i].length[j];
  type_table D;
  create_table(&D, k, k);
  double *H[2];
  for (int i = 0; i < 2; i++) {
    H[i] = (double *)malloc(k*sizeof(double));
    assert(H[i]);
  }

  read_parameters(opt_just_algn, 0, 0, P);
  /*
  for (unsigned int i = 0; i < feature_size; i++) {
    printf("%f\n", P[i]);
  }
  */


  for (unsigned int i = 0; i < num_features; i++)
    P[i] = (double)((int)(0.5+P[i]*1000));


  /*
  memset(P, 0, (num_features+num_examples)*sizeof(double));
  P[gap_offset+gap_size-1] = 15.0/88;
  P[len_offset+len_size-1] = 36.0/88;
  for (unsigned int i = 0; i < gap_size-opt_term_gaps; i++) P[gap_offset+i] = 60.0/88;
  for (unsigned int i = 0; i < len_size-opt_term_gaps; i++) P[len_offset+i] = 38.0/88;
  for (unsigned int i = 0; i < sub_size; i++)
    P[sub_offset+i] = (4.0*(16-protein_matrix[MATRIX_BLOSUM62_3][i]))/88;
  */


  unsigned int *SA = (unsigned int *)malloc((D.row+D.col)*sizeof(unsigned int));
  assert(SA);
  unsigned int *SB = (unsigned int *)malloc((D.row+D.col)*sizeof(unsigned int));
  assert(SB);
  double total = 0;
  for (unsigned int i = 0; i < num_examples; i++) {
    (*work_the_oracle)
      (&D, H, S, confidences, predictions, distributions,
       M, 0, i, 0, 0, S[M[i].a].length[M[i].s]-1,
       S[M[i].a].length[M[i].t]-1, 0, P, &F, 0, opt_term_gaps, opt_term_gaps, &rate[i]);
    total += rate[i];
    for (unsigned int j = 0; j < num_features; j++)
      cost[i] += F[j]*P[j];
    int l = 0;
    get_positions_table(&D, SA, SB, S[M[i].a].length[M[i].s],
			S[M[i].a].length[M[i].t], &l);
    /*
    fprintf(stdout, ";rate=%f, cost=%f\n", rate[i], cost[i]);
    fprintf(stdout, ">%s\n", S[M[i].a].name[M[i].t]);
    write_string(stdout, S[M[i].a].string[M[i].t], SB, l);
    fprintf(stdout, ">%s\n", S[M[i].a].name[M[i].s]);
    write_string(stdout, S[M[i].a].string[M[i].s], SA, l);
    */
  }
  //printf("%f\n", total/num_examples);

  free(SA);
  free(SB);  
  free(H[0]); 
  free(H[1]);
  delete_table(&D);
  free(cost);
  free(rate);
  free(F);
  free(P);
  return 1;
}

int learn_parameters
(const alignment_t *A, 
 const alignment_t *confidences,
 const alignment_t *predictions,
 const dist_t ***distributions,
 unsigned int num_alignments, const char *output) {

  /* Set up the alphabet. */
  int *number = NULL;
  switch (opt_alph_type) {
  case 0: /* Alignments of protein sequences */
    number_alphabet(ALPHABET_PROTEIN_EXT, ALPHABET_SPACE, &number);
    alphabet_size = ALPHABET_SIZE_PROTEIN_EXT;
    break;
  case 1: /* Alignments of DNA sequences */
    number_alphabet(ALPHABET_DNA, ALPHABET_SPACE, &number);
    alphabet_size = ALPHABET_SIZE_DNA;
    break;
  case 2: /* Alignments of RNA sequences */
    number_alphabet(ALPHABET_RNA, ALPHABET_SPACE, &number);
    alphabet_size = ALPHABET_SIZE_RNA;
    break;
  }
  /* Strip off spaces from alignments, and convert to number strings. */
  alignment_t *S = NULL;
  convert_alignments((const int *)number, A, num_alignments, 0, &S);
  free(number);


  /* retrieve information on match points. */
  match_t *M = NULL;
  unsigned int num_examples;
  get_matchlists(A, num_alignments, &M, &num_examples);
  /* set the parameter vector component offsets and sizes.  the diagrams of 
     the feature vector and the parameter vector follow. 
     * feature vector:
     |-----------|-----------|--------------------|
     gap-open    gap-ext     substitutions
     * parameter vector:
     |-----------|-----------|--------------------|---------------|
     gap-open    gap-ext     substitutions        errors
     if terminal gap costs are used, the last ones of gap-open and gap-ext
     are terminal gap-open and gap-ext costs. */
  /* bases (offsets) and sizes for parameter vectors */
  gap_size = opt_disr_lvls+opt_term_gaps;
  len_size = opt_disr_lvls+opt_term_gaps;
  sub_size = NUM_PAIRSW(alphabet_size);
  add_size = opt_prot_pred ? 5+1 : 0;
  err_size = num_examples;
  gap_offset = 0;
  len_offset = gap_offset+gap_size;
  sub_offset = len_offset+len_size;
  add_offset = sub_offset+sub_size;
  err_offset = add_offset+add_size;
  feature_size = gap_size+len_size+sub_size+add_size;

  /*
  printf("%d %d %d %d %d\n", gap_size, len_size, sub_size, add_size, err_size);
  exit(1);
  */
  /* bases (offsets) and sizes for LP */
  lp_gap_size = gap_size;
  lp_len_size = len_size;
  lp_sub_size = (1-opt_fixd_subs)*sub_size;
  lp_add_size = add_size;
  lp_err_size = err_size;
  lp_gap_offset = 1;
  lp_len_offset = lp_gap_offset+lp_gap_size;
  lp_sub_offset = lp_len_offset+lp_len_size;
  lp_add_offset = lp_sub_offset+lp_sub_size;
  lp_err_offset = lp_add_offset+lp_add_size;
  lp_param_size = lp_gap_size+lp_len_size+lp_sub_size+lp_add_size+lp_err_size;
  
  switch (opt_cont_gaps) {
  //case 1: work_the_oracle = &work_the_oracle_1; break;
  //case 2: work_the_oracle = &work_the_oracle_2; break;
  case 3: work_the_oracle = &work_the_oracle_3; break;
  }

  /*
  fprintf(stdout, "============================\n");
  fprintf(stdout, "%d\n", gap_size);
  fprintf(stdout, "%d\n", len_size);
  fprintf(stdout, "%d\n", sub_size);
  fprintf(stdout, "%d\n", add_size);
  fprintf(stdout, "%d\n", err_size);
  fprintf(stdout, "%d\n", gap_offset);
  fprintf(stdout, "%d\n", len_offset);
  fprintf(stdout, "%d\n", sub_offset);
  fprintf(stdout, "%d\n", err_offset);
  fprintf(stdout, "%d\n", feature_size);
  fprintf(stdout, "============================\n");
  fprintf(stdout, "%d\n", lp_gap_size);
  fprintf(stdout, "%d\n", lp_len_size);
  fprintf(stdout, "%d\n", lp_sub_size);
  fprintf(stdout, "%d\n", lp_add_size);
  fprintf(stdout, "%d\n", lp_err_size);
  fprintf(stdout, "%d\n", lp_gap_offset);
  fprintf(stdout, "%d\n", lp_len_offset);
  fprintf(stdout, "%d\n", lp_sub_offset);
  fprintf(stdout, "%d\n", lp_err_offset);
  fprintf(stdout, "%d\n", lp_param_size);
  fprintf(stdout, "============================\n");
  */

  /* retrieve information on core blocks. */
  block_t *B = NULL; 
  if (opt_part_exam) 
    get_coreblocks(S, confidences, predictions, distributions,
		   M, num_examples, &B);
/*
  if (opt_prot_pred)
    attach_hydroinfo_3
      (S, confidences, predictions, distributions, num_alignments);
*/

  /* retrieve information on the given examples. */
  double **W = 0;
  if (opt_init_triv || !opt_part_exam)
    get_wholeblocks_2(A, S, num_alignments, &W, &num_examples);

  /*
    extern double *hval;
    unsigned int i;
  typedef struct _blat {
    double x; char a;
  } blat;

  blat xx[20];
  for (i = 0; i < 20; i++) {
    xx[i].a = (char)i;
    xx[i].x = hval[i];
  }

  int compare(const void *x, const void *y) {
    return 
      (((blat *)x)->x>((blat *)y)->x) ? 1 :
      (((blat *)x)->x<((blat *)y)->x) ? -1 : 0;
  }
  qsort(xx, 20, sizeof(blat), compare);
  */
  /*
  char letter[300];
  memcpy(letter, ALPHABET_PROTEIN_EXT, 30);
  
  
  unsigned int wq;
  extern double *hval;
  for (wq = 0; wq < 20; wq++) {
    printf("%c %d\n", letter[wq], hmap(hval[wq]));
  }


  exit(1);
  

  unsigned int i, j;
  for (i = 0; i < num_examples; i++) {
    for (j = 0; j < feature_size; j++)
      printf("%d ", W[i][j]-B[i].fea[j]);
    printf("\n");
  }
  */
  
  /* this is it. */
  int result;
  if (opt_just_algn) 
    result = apply_parameters
      (S, confidences, predictions, distributions, num_alignments,
       M, B, (const double **)W, output);
  else
    result = find_parameters
      (S, confidences, predictions, distributions, num_alignments, 
       M, B, (const double **)W, output);
  
  /* clean up. */
  free_features(&W, num_examples);
/*
  if (opt_prot_pred)
    detach_hydroinfo(S, num_alignments);
*/
  if (opt_part_exam)
    free_coreblocks(&B, num_examples);
  free_matchlists(&M, num_examples);
  free_alignments(&S, num_alignments);
  return result;
}
