CC=gcc
LD=gcc
AR=ar

CFLAGS=-c -Wall -Werror -O3 -fPIC
LDFLAGS=-shared
ARFLAGS=rcs

NAME=sep
SOURCES=analyse.c #back.c clean.c extract.c lutz.c plist.c refine.c
OBJECTS=$(SOURCES:.c=.o)

all: $(OBJECTS)
#	$(AR) $(ARFLAGS) lib$(NAME).a $(OBJECTS)
	$(LD) $(LDFLAGS) -o lib$(NAME).so $(OBJECTS) -lc -lm

%.o: %.c
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f *.o *.so
