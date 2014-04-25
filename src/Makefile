CC=gcc
LD=gcc
AR=ar

CFLAGS=-c -Wall -Werror -O3 -fPIC
LDFLAGS=-shared
ARFLAGS=rcs

NAME=sep
SOURCES=analyse.c back.c deblend.c extract.c lutz.c util.c

OBJECTS=$(SOURCES:.c=.o)

lib : $(OBJECTS)
#	$(AR) $(ARFLAGS) lib$(NAME).a $(OBJECTS)
	$(LD) $(LDFLAGS) -o lib$(NAME).so $(OBJECTS) -lc -lm

# The following use the implicit rule for making a .o file from .c file:
#     $(CC) $(CPPFLAGS) $(CFLAGS) -c {name}.c -o {name}.o’. 

analyse.o : sep.h extract.h
back.o : sep.h
deblend.o : sep.h extract.h
extract.o : sep.h extract.h
lutz.o : sep.h extract.h
util.o : sep.h

clean :
	rm -f $(OBJECTS) lib$(NAME).so
