CC = g++
CFLAGS = -g -Wall -std=c++11 -O3
LIBS = -lntl -lgmp -lpthread

main: main.o fntru.o ssntru.o general.o
	$(CC) $(CFLAGS) -o main obj/main.o obj/fntru.o obj/ssntru.o obj/general.o $(LIBS)

main.o: src/main.cpp include/fntru.h
	$(CC) $(CFLAGS) -c src/main.cpp -o obj/main.o

fntru.o: src/fntru.cpp include/fntru.h include/ssntru.h include/general.h include/def.h
	$(CC) $(CFLAGS) -c src/fntru.cpp -o obj/fntru.o

ssntru.o: src/ssntru.cpp include/ssntru.h include/general.h include/def.h
	$(CC) $(CFLAGS) -c src/ssntru.cpp -o obj/ssntru.o  

general.o: src/general.cpp include/general.h
	$(CC) $(CFLAGS) -c src/general.cpp -o obj/general.o

clean:
	rm obj/*.o main
