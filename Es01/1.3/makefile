CC = g++
CFLAGS = -Wall -O3 --std=c++11

SN1.3.exe : SN1.3.o random.o
	$(CC) random.o SN1.3.o -o SN1.3.exe
SN1.3.o : SN1.3.cpp
	$(CC) -c SN1.3.cpp -o SN1.3.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o SN1.3.exe seed.out
