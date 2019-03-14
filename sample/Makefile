CC = gcc

LIBS = -ltepla -lcrypto -lgmp

LIBDIR = -L../lib
INCDIR = -I../include

CFLAGS  = -Wall -g -O4
LDFLAGS = -Wall -g

#TARGET = sample.exe speed.exe ecdsa.exe sample2.exe
TARGET = sample2.exe

all: $(TARGET)

%.exe: %.o
	$(CC) $(LDFLAGS) -o $@ $< $(INCDIR) $(LIBDIR) $(LIBS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< $(INCDIR)

.PHONY: clean
clean:
	$(RM) *~ $(TARGET)
