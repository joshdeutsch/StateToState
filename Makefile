COMPILER = gcc
CFLAGS = -g -O3 
LIBRARIES = -lgsl -lgslcblas -lm -lc
C_OBJECTS = $(patsubst %.c,%.o,$(wildcard *.c))
LIBDIR =
all: lmin

$(C_OBJECTS): %.o: %.c
	$(COMPILER) -c $(CFLAGS)  $< -o $@

lmin:  $(C_OBJECTS)
	$(COMPILER) $(CFLAGS) -o $@  $^  $(LIBDIR) $(LIBRARIES) 

clean:
	rm -f *.o lmin
