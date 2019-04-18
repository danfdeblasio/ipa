#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include "inverse.h"
#include "error.h"

char *input  = 0;
char *output = 0;
int parse_options(int *argc, char **argv[]);

int main(int argc, char *argv[]) {
  parse_options(&argc, &argv);
  if (!argc) error("specify multiple alignment(s).");

  alignment_t *alignments = 0;
  get_alignments((const char **)argv, (unsigned int)argc, &alignments);
  alignment_t *conf = 0;
  alignment_t *pred = 0;
  dist_t    ***dist = 0;

  extern char opt_prot_pred;
  extern char opt_disr_meas;

  /* call the IPA. */
  learn_parameters(alignments, conf, pred, (const dist_t ***)dist,
		   (unsigned int)argc, (const char *)output);
  
  /* clean up. */
  free_alignments(&alignments, (unsigned int)argc);
	
  return 0;
}
