DEBUG?=false

all: main.exe

main.exe: ${DEBUG}-main.o
	g++ ${DEBUG}-main.o -o main.exe

${DEBUG}-main.o: main.cpp
	g++ -c -D"DEBUG=${DEBUG}" main.cpp -o ${DEBUG}-main.o

clean:
	rm *-main.o main.exe
