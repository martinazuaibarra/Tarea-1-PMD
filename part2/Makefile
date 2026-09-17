CC = gcc
CFLAGS = -O2 -Wall -Wextra
LDFLAGS_PTHREAD = -pthread

all: mapper mapper_pthreads

mapper: mapper.c
	$(CC) $(CFLAGS) -o mapper mapper.c

mapper_pthreads: mapper_pthreads.c
	$(CC) $(CFLAGS) $(LDFLAGS_PTHREAD) -o mapper_pthreads mapper_pthreads.c

clean:
	rm -f mapper mapper_pthreads

.PHONY: all clean
