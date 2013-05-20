# no need to use SFMT.c as a target, it's included in the code
# to use the inline functions
#SRCS=daisies.c 

TARGETS=daisies #daisies2

GCC=/usr/bin/gcc
# as suggested by SFMT how-to-compile document
CFLAGS+=-O3 -fno-strict-aliasing -fgnu89-inline -std=c99
LIBS=-lm -lncurses -lpthread #-Lfcmp

#####

INCLUDES=-I SFMT/ -I fcmp/
# MEXP requested by SFMT
DEFS=-D MEXP=19937 
# add to the sources 
SRCS+=SFMT/SFMT.c

#####

# with console output
DEFS+= -D _WITH_INFO
#CFLAGS+=-g

all: $(TARGETS)

daisies: clean
	$(GCC) $(CFLAGS) $(LIBS) $(DEFS) $(INCLUDES) -o $@ $@.c

daisies2: clean
	$(GCC) $(CFLAGS) $(LIBS) $(DEFS) $(INCLUDES) -o $@ $@.c

clean: 
	rm -f $(TARGETS) 

