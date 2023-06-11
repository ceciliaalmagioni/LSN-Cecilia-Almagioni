CC = g++
CFLAGS = -Wall -O3 --std=c++11

SN2.2.exe : SN2.2.o random.o
	$(CC) random.o SN2.2.o -o SN2.2.exe
SN2.2.o : SN2.2.cpp
	$(CC) -c SN2.2.cpp -o SN2.2.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o SN2.2.exe seed.out
