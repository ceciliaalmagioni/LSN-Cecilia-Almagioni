CC = mpic++
CFLAGS = -Wall -O3 --std=c++11 

SN10.exe : SN10.o random.o
	$(CC) random.o SN10.o -o SN10.exe $(LFLAGS)
SN10.o : SN10.cpp
	$(CC) -c SN10.cpp -o SN10.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o SN10.exe seed.out
