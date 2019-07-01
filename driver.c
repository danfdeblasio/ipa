#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include "inverse.h"
#include "error.h"

char *input  = 0;
char *output = 0;
char *opt_flist = 0;
int parse_options(int *argc, char **argv[]);

int main(int argc, char *argv[]) {
  parse_options(&argc, &argv);
  if (!argc) error("specify multiple alignment(s).");

  alignment_t *alignments = 0;
  int num_aln = argc;
  if(opt_flist != 0){
    num_aln = atoi(argv[0]);
    FILE *list;
    list = fopen(opt_flist, "r");
    char* file_names[num_aln];
    for(int i=0; i<num_aln; i++){
      printf("i: %d\n", i);
      file_names[i] = malloc(sizeof(char) * 200);
      fscanf(list, "%s", file_names[i]);
      printf("Read fname %s\n",file_names[i]);
    }
    get_alignments((const char **) file_names, (unsigned int) num_aln, &alignments);
  }else{
    get_alignments((const char **)argv, (unsigned int)num_aln, &alignments);
  }
  alignment_t *conf = 0;
  alignment_t *pred = 0;
  dist_t    ***dist = 0;

  extern char opt_prot_pred;
  extern char opt_disr_meas;

  /* call the IPA. */
  learn_parameters(alignments, conf, pred, (const dist_t ***)dist,
		   (unsigned int)num_aln, (const char *)output);
  
  /* clean up. */
  //free_alignments(&alignments, (unsigned int)argc);
  free_alignments(&alignments, (unsigned int)num_aln);
	
  return 0;
}
