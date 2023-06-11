CC = g++
CFLAGS = -Wall -O3 --std=c++11

SN1.1.exe : SN1.1.o random.o
	$(CC) random.o SN1.1.o -o SN1.1.exe
SN1.1.o : SN1.1.cpp
	$(CC) -c SN1.1.cpp -o SN1.1.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o SN1.1.exe seed.out
