#Makefile
OBJS = tools.o direct.o algorithm.o io.o
CC = g++ 
DEBUG = -g
CFLAGS =   -c $(DEBUG)
LFLAGS =   $(DEBUG)

diretto.out : $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o directACE_spin_pm1.out

tools.o: ./src/tools.h ./src/tools.cpp ./src/dataStructures.h 
	$(CC) $(CFLAGS) ./src/tools.cpp

direct.o: ./src/direct.h ./src/direct.cpp ./src/tools.h ./src/dataStructures.h
	$(CC) $(CFLAGS) ./src/direct.cpp

io.o: ./src/dataStructures.h ./src/io.cpp ./src/tools.h
	$(CC) $(CFLAGS) ./src/io.cpp

algorithm.o: ./src/algorithm.cpp ./src/algorithm.h ./src/direct.h ./src/tools.h
	$(CC) $(CFLAGS) ./src/algorithm.cpp

clean:
	rm -f *.o *~
