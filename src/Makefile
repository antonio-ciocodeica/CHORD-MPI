CC = mpicc
TARGET = chord
CFLAGS = -Wall -Wextra -O2
SRC = chord.c

all: build

build: $(SRC)
	$(CC) $(CFLAGS) $(SRC) -o $(TARGET)

clean:
	rm -f $(TARGET) *.o

.PHONY: all clean run
