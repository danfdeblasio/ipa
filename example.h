#ifndef __EXAMPLE__
#define __EXAMPLE__

//#include "nonstd.h"

typedef struct alignment_type {
  unsigned int count;      /* the number of sequences */
  unsigned int *length;    /* list of sequence lengths */
  char         **name;     /* list of sequence names */
  char         **string;   /* list of sequences */
  char         **disrupt;  /* list of sequences of local disruption levels */
  char         **useful;   /* list of sequences of some disruption levels */
} alignment_t;

typedef struct point_type {
  unsigned int x;
  unsigned int y;
} point_t;

typedef struct dist_type {
  float c, a, b;
} dist_t;

typedef struct match_type {
  unsigned int a;
  unsigned int s;
  unsigned int t;
  unsigned int length;
  point_t     *list;
} match_t;

typedef struct block_type {
  unsigned int count; /* number of non-core blocks */
  point_t *pos;   /* non-core block beginning positions */
  point_t *len;   /* non-core block lengths */
  double  *fea;   /* feature counts of core blocks */
} block_t;
 
void get_alignments
(const char **path, unsigned int num_alignments, alignment_t **A);
void get_alignments_hydro
(const char **path, unsigned int num_alignments, alignment_t **A);
void free_alignments
(alignment_t **A, unsigned int num_alignments);

void make_alphabet
(const char *letter, const char *space, char **alphabet);
void number_alphabet
(const char *letter, const char *space, int **number);
void convert_alignments
(const int *map, const alignment_t *A, unsigned int num_alignments, char space, alignment_t **S);
void get_matchlists
(const alignment_t *A, unsigned int num_alignments, match_t **M, unsigned int *num_examples);
void free_matchlists
(match_t **M, unsigned int num_blocks);

void get_coreblocks
(const alignment_t *S, const alignment_t *confidences, const alignment_t *predictions, const dist_t ***distributions, const match_t *M, unsigned int num_examples, block_t **B);
//(const alignment_t *S, const match_t *M, uint num_examples, block_t **B);
void free_coreblocks
(block_t **B, unsigned int num_blocks);
void get_wholeblocks
(const alignment_t *A, unsigned int num_alignments, const int *alphabet, double ***F, unsigned int *num_examples);
void get_wholeblocks_2
(const alignment_t *A, const alignment_t *T, unsigned int num_alignments, double ***F, unsigned int *num_examples);
void free_features
(double ***F, int num_features);

#endif
