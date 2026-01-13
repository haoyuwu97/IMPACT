
TARGET	 = CAP_v6.0.1

OBJS     = main.o input.o calculate.o cluster.o volume_cry.o output.o conformation.o

HEADERS  = crystal.h

CC       = g++
CFLAGS   = -O3 -Wall
LDFLAGS  =
MKFILE   = Makefile
TEST_BIN = tests/test_input
BENCH_BIN = tests/benchmark

$(TARGET) : $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS)

.cpp.o:
	$(CC) $(CFLAGS) -c $<
clean:
	@-rm -f $(OBJS)
	@-rm -f $(TEST_BIN)
	@-rm -f $(BENCH_BIN)

test: $(TEST_BIN)
	./$(TEST_BIN)

$(TEST_BIN): tests/test_input.cpp input.cpp $(HEADERS)
	$(CC) $(CFLAGS) -I. tests/test_input.cpp input.cpp -o $(TEST_BIN)

bench: $(BENCH_BIN)
	./$(BENCH_BIN)

$(BENCH_BIN): tests/benchmark.cpp input.cpp $(HEADERS)
	$(CC) $(CFLAGS) -I. tests/benchmark.cpp input.cpp -o $(BENCH_BIN)

# dependencies
main.o: main.cpp crystal.h 
input.o: input.cpp crystal.h 
calculate.o: calculate.cpp crystal.h 
cluster.o: cluster.cpp crystal.h 
volume_cry.o: volume_cry.cpp crystal.h 
output.o: output.cpp crystal.h 
conformation.o: conformation.cpp crystal.h 
