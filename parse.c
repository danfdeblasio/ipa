#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include "inverse.h"

int parse_options(int *argc, char **argv[]) {
  char sparam[] = "DRCg:m:s:o:w:n:l:TtN:fp:x:SI:O:dv";
  struct option lparam[] = {
//    { "DNA",                 no_argument,       0, 'D' },
//    { "RNA",                 no_argument,       0, 'R' },
//    { "complete",            no_argument,       0, 'C' },
//    { "context",             required_argument, 0, 'g' },
//    { "disruption",          required_argument, 0, 'm' },
//    { "prediction",          required_argument, 0, 's' },
//    { "convexity",           required_argument, 0, 'w' },
//    { "window",              required_argument, 0, 'n' },
//    { "levels",              required_argument, 0, 'l' },
//    { "trivial-completion",  no_argument,       0, 'T' },
    { "no-terminal-gaps",    no_argument,       0, 't' },
//    { "normalization",       required_argument, 0, 'N' },
    { "fixed-substitution",  no_argument,       0, 'f' },
    { "pairwise-alignments", required_argument, 0, 'p' },
    { "iterations",          required_argument, 0, 'x' },
    { "show-statistics",     no_argument,       0, 'S' },
    { "input",               required_argument, 0, 'I' },
    { "output",              required_argument, 0, 'O' },
    { "debug",               no_argument,       0, 'd' },
    { "verbose",             no_argument,       0, 'v' },
    { "file-list",            required_argument, 0, 'l'},
    { 0,                     0,                 0,  0  }
  };

  extern char  *input;
  extern char  *output;
  extern char   opt_verb;
  extern char   opt_debg;
  extern char   opt_stat;
  extern char  *opt_flist;
  //extern char   opt_part_exam;
  //extern char   opt_alph_type;
  //extern char   opt_disr_meas;
  //extern char   opt_cont_gaps;
  extern char   opt_fixd_subs;
  //extern char   opt_disr_lvls;
  //extern char   opt_init_triv;
  extern char  *opt_just_algn;
  extern int   opt_maxm_itrs;
  //extern char   opt_norm_cost;
  //extern double opt_conv_wght;
  //extern char   opt_prot_pred;
  extern int   opt_term_gaps;
  //extern char   opt_wind_size;
  
  int opt;
  while (EOF != (opt = getopt_long(*argc, *argv, sparam, lparam, NULL)))
    switch (opt) {
/*
    case 'D': opt_alph_type = 1; opt_prot_pred = 0;
      opt_part_exam = 0; opt_disr_meas = 0; opt_cont_gaps = 0;
      opt_fixd_subs = 0; opt_disr_lvls = 1; opt_maxm_itrs = 1;
      opt_term_gaps = 1; opt_wind_size = 1; opt_init_triv = 1;
      break;
    case 'R': opt_alph_type = 2; opt_prot_pred = 0; 
      opt_part_exam = 0; opt_disr_meas = 0; opt_cont_gaps = 0;
      opt_fixd_subs = 0; opt_disr_lvls = 1; opt_maxm_itrs = 1;
      opt_term_gaps = 1; opt_wind_size = 1; opt_init_triv = 1;
      break;
    case 'C': opt_part_exam = 0; break;
    case 'g': opt_cont_gaps = (char)atoi(optarg); break;
    case 'm': opt_disr_meas = (char)atoi(optarg); break;
    case 's': opt_prot_pred = (char)atoi(optarg); break;
    case 'w': opt_conv_wght = atof(optarg); break;
    case 'n': opt_wind_size = (char)atoi(optarg); break;
    case 'l': opt_disr_lvls = (char)atoi(optarg); break;
    case 'T': opt_init_triv = 1; break;
*/
    case 't': opt_term_gaps = 0; break;
//    case 'N': opt_norm_cost = (char)atoi(optarg); break;
    case 'f': opt_fixd_subs = 1; break;
    case 'p': opt_just_algn = optarg; break;
    case 'x': opt_maxm_itrs = atoi(optarg); break;
    case 'l': opt_flist = optarg; break;
    case 'S': opt_stat = 1; break;
    case 'I': input  = optarg; break;
    case 'O': output = optarg; break;
    case 'd': opt_debg = 1; break;
    case 'v': opt_verb = 1; break;
    }
  *argc -= optind;
  *argv += optind;
  return 1;
}
