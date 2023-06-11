CC = g++
CFLAGS = -Wall -O3 --std=c++11

SN8.exe : SN8.o random.o
	$(CC) random.o SN8.o -o SN8.exe
SN8.o : SN8.cpp
	$(CC) -c SN8.cpp -o SN8.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o SN8.exe seed.out
