CC=gcc
CFLAGS=-I../src -L../src -Wall -Wl,-rpath,../src

runtests : runtests.c
	$(CC) $(CFLAGS) -o runtests runtests.c -lm -lsep
