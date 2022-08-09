FLAGS=-O2 -fopenmp

CC=gcc

RM=rm -f

EXEC=dijkstra

all: $(EXEC)

$(EXEC):
        $(CC) $(FLAGS) $(EXEC).c -c -o $(EXEC).o
        $(CC) $(FLAGS) $(EXEC).o -o $(EXEC)

run:
        ./$(EXEC)

clean:
        $(RM) $(EXEC).o $(EXEC)

new:
        make clean
        make all
        make run