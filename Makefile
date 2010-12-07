BIN	:= histogram
SRC	:= histogram_nd.c
OBJ	:= $(SRC:.c=.o)

CC	:= gcc
OFLAGS	:= -pipe -O3
CFLAGS	:= $(OFLAGS)

.PHONY:	clean

all:	$(BIN)

clean:
	rm -f $(OBJ)

cleanall: clean
	rm -f $(BIN) *~

$(BIN): $(OBJ)
	$(CC) $(CFLAGS) -lm -o $@ $^
