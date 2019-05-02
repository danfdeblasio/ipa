IPA = ipa
####################################
CC = gcc -std=c99
#CC = g++
LIBS = -lm -lglpk
CFLAGS = -Wno-sign-compare -W -g 
#CFLAGS = -O3 -DNDEBUG
OBJS = lpx.o driver.o example.o inverse.o extend.o matrix.o error.o fasta.o parse.o extend3.o

all: $(IPA)

$(IPA): $(OBJS)
	${CC} $(OBJS) $(LIBS) -o $(IPA)

%.o: %.c
	$(CC) $(CFLAGS) -c $<

clean: 
	\rm -rf $(OBJS) $(IPA) $(IPA).exe *~ a.out core
