#ifndef __FASTA__
#define __FASTA__

//#include "nonstd.h"
#include <stdio.h>

void read_fasta(const char *path, char gap, char ***name, char ***string, unsigned int *count);
void write_fasta(FILE *stream, char gap, size_t col, const char **name, const char **string, unsigned int count);

#endif
