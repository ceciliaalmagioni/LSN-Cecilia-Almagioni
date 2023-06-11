CC = g++
CFLAGS = -Wall -O3 --std=c++11

SN3.exe : SN3.o random.o
	$(CC) random.o SN3.o -o SN3.exe
SN3.o : SN3.cpp
	$(CC) -c SN3.cpp -o SN3.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o SN3.exe seed.out
