CC = g++
CFLAGS = -Wall -O3 --std=c++11

SN2.1.exe : SN2.1.o random.o
	$(CC) random.o SN2.1.o -o SN2.1.exe
SN2.1.o : SN2.1.cpp
	$(CC) -c SN2.1.cpp -o SN2.1.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o SN2.1.exe seed.out
