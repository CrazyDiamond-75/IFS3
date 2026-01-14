USE_OPENMP = 1

#------------------
# Compiler flags
#------------------

CC = gcc
RM = rm -f

CFLAGS = -Wall -O3 -funroll-loops -march=native -mtune=native -std=c23
LIBS   = -lm
    
ifdef USE_OPENMP
  CFLAGS += -DUSE_OPENMP -fopenmp
endif

SOURCES_LIB :=	\
	
HEADER_LIB := \
	$(SOURCES_LIB:.c=.h)

OBJECTS_LIB := \
	$(SOURCES_LIB:.c=.o)	

SOURCES_PROGRAMS :=	\
	ifs3_cpu.c

PROGRAMS := \
	$(SOURCES_PROGRAMS:.c=)

all: $(PROGRAMS)

$(OBJECTS_LIB): %.o: %.c
	@echo Compiling \"$<\"
	$(CC) $(CFLAGS) -c $< -o $@

$(PROGRAMS): %: %.c $(OBJECTS_LIB)
	@echo Compiling and linking \"$<\"
	$(CC) $< $(CFLAGS) $(HEADER_LIB) -o $@ $(LDFLAGS) $(OBJECTS_LIB) $(LIBS)

clean:
	$(RM) $(PROGRAMS) $(OBJECTS_LIB)
