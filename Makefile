
# the compiler to use

CC=icpc
#CC=g++

#compilation flags
CFLAGS=-opt-prefetch -fno-alias -fno-exceptions -Wall -wd981 -wd193 -fp-model fast=2 -fast-transcendentals -fast 
#CFLAGS=-Wall -pedantic -O3

TARGETS=estimation \
estimation_debug \
estimation_sim \
estimation_test \
estimation_sim_full \
estimation_sim_fullrent \
estimation_sim_fullwage \
estimation_stdev  \
estimation_sim4 \
estimation_sim4_full \
estimation_sim4_fullrent \
estimation_sim4_fullwage \
estimation_test \
estimation_test_full \
estimation_test_fullwage \
estimation_wage_selection \
estimation_fullwage_selection

all: $(TARGETS)

estimation: estimation.cpp
	@echo Building $@...;
	@$(CC) $(CFLAGS) -DINFO -o estimation estimation.cpp;

estimation_debug: estimation.cpp
	@echo Building $@...;
	@$(CC) -DTRACE_LOAD -DTRACE -DINFO -Wall -g -O0 -o estimation_debug estimation.cpp;

estimation_sim: estimation.cpp
	@echo Building $@...;
	@$(CC) $(CFLAGS) -DSIMULATION -DTRACE -DINFO -o estimation_sim estimation.cpp

estimation_sim_full: estimation.cpp
	@echo Building $@...;
	@$(CC) $(CFLAGS) -DSIMULATION -DTRACE -DINFO -DFULL_TRACE -DFULL_TRACE_INDEX -o estimation_sim_full estimation.cpp

estimation_sim_fullrent: estimation.cpp
	@echo Building $@...;
	@$(CC) $(CFLAGS) -DSIMULATION -DTRACE -DINFO -DFULL_TRACE -DFULL_TRACE_RENT -o estimation_sim_fullrent estimation.cpp

estimation_sim_fullwage: estimation.cpp
	@echo Building $@...;
	@$(CC) $(CFLAGS) -DSIMULATION -DTRACE -DINFO -DFULL_TRACE -DFULL_TRACE_WAGE -o estimation_sim_fullwage estimation.cpp

estimation_sim4: estimation_sim4.cpp
	@echo Building $@...;
	@$(CC) $(CFLAGS) -DSIMULATION -DTRACE -DINFO -o estimation_sim4 estimation_sim4.cpp

estimation_sim4_full: estimation_sim4.cpp
	@echo Building $@...;
	@$(CC) $(CFLAGS) -DSIMULATION -DTRACE -DINFO -DFULL_TRACE -DFULL_TRACE_INDEX -o estimation_sim4_full estimation_sim4.cpp

estimation_sim4_fullrent: estimation_sim4.cpp
	@echo Building $@...;
	@$(CC) $(CFLAGS) -DSIMULATION -DTRACE -DINFO -DFULL_TRACE -DFULL_TRACE_RENT -o estimation_sim4_fullrent estimation_sim4.cpp

estimation_sim4_fullwage: estimation_sim4.cpp
	@echo Building $@...;
	@$(CC) $(CFLAGS) -DSIMULATION -DTRACE -DINFO -DFULL_TRACE -DFULL_TRACE_WAGE -o estimation_sim4_fullwage estimation_sim4.cpp

estimation_stdev: estimation.cpp
	@echo Building $@...;
	@$(CC) $(CFLAGS) -DCALC_STDEV -DINFO -o estimation_stdev estimation.cpp

estimation_test: estimation.cpp
	@echo Building $@...;
	@$(CC) $(CFLAGS) -DTRACE -DINFO -o estimation_test estimation.cpp

estimation_test_full: estimation.cpp
	@echo Building $@...;
	@$(CC) $(CFLAGS) -DTRACE -DINFO -DFULL_TRACE -DFULL_TRACE_INDEX -o estimation_test_full estimation.cpp

estimation_test_fullwage: estimation.cpp
	@echo Building $@...;
	@$(CC) $(CFLAGS) -DTRACE -DINFO -DFULL_TRACE -DFULL_TRACE_WAGE -o estimation_test_fullwage estimation.cpp

estimation_wage_selection: estimation_wage_selection.cpp
	@echo Building $@...;
	@$(CC) $(CFLAGS) -DTRACE -DINFO -DWAGE_SELECTION -o estimation_wage_selection estimation_wage_selection.cpp

estimation_fullwage_selection: estimation_wage_selection.cpp
	@echo Building $@...;
	@$(CC) $(CFLAGS) -DTRACE -DINFO -DWAGE_SELECTION -DFULL_TRACE -DFULL_TRACE_WAGE -o estimation_fullwage_selection estimation_wage_selection.cpp

clean:
	@echo Removing all executables...
	@rm -rf $(TARGETS)

