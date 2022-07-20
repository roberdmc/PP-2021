FLAGS=-O2 -fopenmp

CC=gcc

RM=rm -f

EXEC=dijkstra

all: $(EXEC)

$(EXEC):
	$(CC) $(FLAGS) $(EXEC).c -c -o exc$(EXEC).o
	$(CC) $(FLAGS) exc$(EXEC).o -o exc$(EXEC)

run:
	./exc$(EXEC)

clean:
	$(RM) exc$(EXEC).o exc$(EXEC)

new:
	make clean
	make all
	make run