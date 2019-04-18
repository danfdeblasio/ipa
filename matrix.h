#ifndef __MATRIX__
#define __MATRIX__

#define ALPHABET_SPACE            "-"
#define ALPHABET_PROTEIN          "ARNDCQEGHILKMFPSTWYV"
#define ALPHABET_PROTEIN_EXT      "ARNDCQEGHILKMFPSTWYVBZX" 
#define ALPHABET_DNA              "ACGT"
#define ALPHABET_RNA              "ACGU" 
#define ALPHABET_SIZE_PROTEIN     20
#define ALPHABET_SIZE_PROTEIN_EXT 23
#define ALPHABET_SIZE_DNA         4
#define ALPHABET_SIZE_RNA         4

#define HYDRO_SIZE            ALPHABET_SIZE_PROTEIN_EXT
#define HYDRO_ZHOUZHOU        0
#define HYDRO_ZHOUZHOU_2      1
#define HYDRO_ZHOUZHOU_3      2
#define HYDRO_KYTEDOOLITTLE   3
#define HYDRO_KYTEDOOLITTLE_2 4
#define HYDRO_PASCARALLEARGOS 5
#define HYDRO_CRUSTALW_2      6

#define PRIOR_SIZE         ALPHABET_SIZE_PROTEIN_EXT
#define PRIOR_PIR          0
#define PRIOR_BLOSUM62     1

#define JOINT_BLOSUM62     0

#define MATRIX_SIZE						\
  (ALPHABET_SIZE_PROTEIN_EXT*(ALPHABET_SIZE_PROTEIN_EXT+1)>>1)
#define MATRIX_BLOSUM22_7  0
#define MATRIX_BLOSUM22    0
#define MATRIX_BLOSUM35    1
#define MATRIX_BLOSUM40    2
#define MATRIX_BLOSUM45    3
#define MATRIX_BLOSUM50    4
#define MATRIX_BLOSUM50_2  5
#define MATRIX_BLOSUM55    6
#define MATRIX_BLOSUM60    7
#define MATRIX_BLOSUM62    8
#define MATRIX_BLOSUM62_3  9
#define MATRIX_BLOSUM62_4  10
#define MATRIX_BLOSUM65    11
#define MATRIX_BLOSUM70    12
#define MATRIX_BLOSUM75    13
#define MATRIX_BLOSUM80    14
#define MATRIX_BLOSUM80_3  15
#define MATRIX_BLOSUM85    16
#define MATRIX_BLOSUM90    17
#define MATRIX_BLOSUM95    18
#define MATRIX_BLOSUM100   19
#define MATRIX_BLOSUM100_3 20
#define MATRIX_PAM120      21
#define MATRIX_PAM160      22
#define MATRIX_PAM250      23

#endif
