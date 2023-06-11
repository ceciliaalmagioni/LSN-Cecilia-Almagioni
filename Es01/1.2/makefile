CC = g++
CFLAGS = -Wall -O3 --std=c++11

SN1.2.exe : SN1.2.o random.o
	$(CC) random.o SN1.2.o -o SN1.2.exe
SN1.2.o : SN1.2.cpp
	$(CC) -c SN1.2.cpp -o SN1.2.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o SN1.2.exe seed.out
