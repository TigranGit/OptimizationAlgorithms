all: main.exe

main.exe: main.o
	g++ main.o -o main.exe

main.o: main.cpp
	g++ -c main.cpp

clean:
	rm main.o main.exe
