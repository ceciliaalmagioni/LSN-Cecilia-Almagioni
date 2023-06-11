CC = g++
CFLAGS = -Wall -O3 --std=c++11

SN9.exe : SN9.o random.o
	$(CC) random.o SN9.o -o SN9.exe
SN9.o : SN9.cpp
	$(CC) -c SN9.cpp -o SN9.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o SN9.exe seed.out
