# Makefile for CellDivide
CC = gcc
CFLAGS = -std=c99
OBJECT = Divide.o Main.o Pileup.o Parse.o

celldivide: $(OBJECT)
	$(CC) $(CFLAGS) $(OBJECT) -o $@

Divide.o: Share.h
Main.o: Share.h
Pileup.o: Share.h
Parse.o: Share.h

.PHONY : clean
clean:
	rm -f $(OBJECT)
