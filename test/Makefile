CC=gcc
CFLAGS=-I../src -L../src -Wall -Wl,-rpath,../src

testfits : testfits.c
	$(CC) $(CFLAGS) -o testfits testfits.c -lcfitsio -lsep

all : testfits