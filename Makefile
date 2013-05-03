CC=g++
CFLAGS=-Wall -O3 --std=c++11

all:
	$(CC) $(CFLAGS) example_lotka.cpp -o example_lotka.exe
	$(CC) $(CFLAGS) example_lotka_events.cpp -o example_lotka_events.exe

clean:
	rm -rf *.exe
