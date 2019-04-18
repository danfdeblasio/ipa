#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include "error.h"
//#include "nonstd.h"

long filesize(FILE *stream) {
  /* do some error checking. */
  long curr = ftell(stream);
  fseek(stream, 0L, SEEK_END);
  long size = ftell(stream);
  fseek(stream, curr, SEEK_SET);
  return size-curr;
}

void read_fasta(const char *path, char gap, char ***name, char ***string, 
		unsigned int *count) {
  FILE *stream = fopen(path, "r");
  assert(stream);
  long  size   = filesize(stream);
  char *temp   = (char *)malloc(size);
  assert(temp);
  char *buffer = (char *)malloc(size);
  assert(buffer);
  long  offset = 0;

  *count = 0;  
  while (!feof(stream)) {
    /* read the name of sequence */
    char *p, ch = fgetc(stream);
    if (isspace(ch)) continue;
    if (ch != '>') error("%s is not fasta formatted.", path);
    fgets(&temp[offset], size-offset, stream);
    if ((p = strchr(&temp[offset], '\n')) != 0) *p = 0;
    offset += strlen(&temp[offset]);    
    temp[offset++] = 0;

    /* read the sequence */
    long curr = offset;
    while (!feof(stream)) {
      ch = fgetc(stream);
      ungetc(ch, stream);
      if (ch == '>') break;
      fgets(buffer, size-offset, stream);
      if (ch == ';') continue;
      if ((p = strchr(buffer, '\n')) != 0) *p = 0;      
      unsigned int j = strlen(buffer);
      for (unsigned int i = 0; i < j; i++) {
	if (feof(stream)) break;
	if (isspace(buffer[i])) continue;
	if (buffer[i] == '-' && !gap) continue;
	temp[offset++] = buffer[i];
      }
    }
    if (curr == offset) error("there's no sequence in %s.", path);
    temp[offset++] = 0;
    (*count)++;
  }
  fclose(stream);
  free(buffer);

  /* allocate memory for names and sequences. */
  *name   = (char **)malloc(*count*sizeof(char *));
  assert(*name);
  *string = (char **)malloc(*count*sizeof(char *));
  assert(*string);
  unsigned int i;
  for (offset = i = 0; i < *count; i++) {
    unsigned int len = strlen(&temp[offset])+1;
    (*name)[i] = (char *)malloc(len);
    assert((*name)[i]);
    strcpy((*name)[i], &temp[offset]);
    offset += len;
    len = strlen(&temp[offset])+1;
    (*string)[i] = (char *)malloc(len);
    assert((*string)[i]);
    strcpy((*string)[i], &temp[offset]);
    offset += len;
    //(*name)[i]   = strdup(&temp[offset]);
    //offset += strlen(&temp[offset])+1;
    //(*string)[i] = strdup(&temp[offset]);
    //offset += strlen(&temp[offset])+1;
  }
  free(temp);
  if (!gap) return;

  return;

  /* check the validity of alignment. */
  if (*count < 2) error("%s is a invalid alignment.", path);
  for (i = 1; i < *count; i++)
    if (strlen((*string)[0]) != strlen((*string)[i]))
      error("not all lengths of sequences are the same in %s.", path);
}

void write_fasta(FILE *stream, char gap, unsigned int col, const char **name, const char **string, unsigned int count) {
  for (unsigned int i = 0; i < count; i++) {
    fprintf(stream, ">%s\n", name[i]);
    unsigned int k = 0;
    for (unsigned int j = 0; j < strlen(string[i]); j++) {
      if (string[i][j] == '-' && !gap) continue;
      fprintf(stream, "%c", string[i][j]);
      k++;
      if (k%col == 0) fprintf(stream, "\n");
    }
    if (k%col) fprintf(stream, "\n");
  }
}
