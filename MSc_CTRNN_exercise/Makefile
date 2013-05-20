# no need to use SFMT.c as a target, it's included in the code
# to use the inline functions
SRCS=ctrnn.c matop.c 
GCC=/usr/bin/gcc
#CFLAGS=-g
# as suggested by SFMT how-to-compile document
CFLAGS=-O3 -fno-strict-aliasing 
LIBS=-lm -lncurses 
TARGETS=ctrnn
# MEXP requested by SFMT
DEFS=-D MEXP=19937 

# with ncurses output
#DEFS+=-D _WITH_CURSES
# with console output
DEFS+= -D _WITH_INFO
# with lots of useless debugging information
#DEFS+= -D _WITH_DEBUG

all: $(TARGETS)

ctrnn: clean
	$(GCC) $(CFLAGS) $(LIBS) $(DEFS) -o $@ $(SRCS)

clean: 
	rm -f $(TARGETS) 

test_rand: clean_test_rand
	$(GCC) $(CFLAGS) $(LIBS) $(DEFS) -o $@ test_rand.c

test_mat: clean_test_mat
	$(GCC) $(CFLAGS) $(LIBS) $(DEFS) -o $@ test_mat.c matop.c

clean_test_rand: 
	rm -f test_rand

clean_test_mat:
	rm -f test_mat
