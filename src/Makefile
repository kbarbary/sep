CC=gcc
LD=gcc
AR=ar

# Implicit rule for making .o file from .c file:
#     $(CC) $(CPPFLAGS) $(CFLAGS) -c {name}.c -o {name}.o’. 

CFLAGS=-c -Wall -Werror -O3 -fPIC
LDFLAGS=-shared
ARFLAGS=rcs

NAME=sep
SOURCES=analyse.c back.c #clean.c extract.c lutz.c plist.c refine.c

OBJECTS=$(SOURCES:.c=.o)

lib: $(OBJECTS)
#	$(AR) $(ARFLAGS) lib$(NAME).a $(OBJECTS)
	$(LD) $(LDFLAGS) -o lib$(NAME).so $(OBJECTS) -lc -lm

# The following use the implicit rule for making a .o file from .c file:
#     $(CC) $(CPPFLAGS) $(CFLAGS) -c {name}.c -o {name}.o’. 

analyse.o: sep.h defs.h extract.h
back.o: sep.h defs.h

clean:
	rm -f $(OBJECTS) lib$(NAME).so
