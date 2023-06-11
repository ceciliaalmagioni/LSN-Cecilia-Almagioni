CC = g++
CFLAGS = -Wall -O3 --std=c++11

SN5.exe : SN5.o random.o
	$(CC) random.o SN5.o -o SN5.exe
SN5.o : SN5.cpp
	$(CC) -c SN5.cpp -o SN5.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o SN5.exe seed.out
