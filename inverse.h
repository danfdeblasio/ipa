#ifndef __INVERSE__
#define __INVERSE__

#include "example.h"

int learn_parameters
(const alignment_t *A, const alignment_t *confidences, const alignment_t *predictions, const dist_t ***distributions, unsigned int num_alignments, const char *output);

#endif
