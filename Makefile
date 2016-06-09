# the compiler to use
ifeq ($(CC),cc)
	CC=g++
endif

ifeq ($(CC),icpc)
	CWARNING_FLAGS=-Wall -wd981 -wd193 -wd2259 -wd1572
	CFLAGS=$(CWARNING_FLAGS) -opt-prefetch -fno-alias -fno-exceptions -fp-model fast=2 -fast-transcendentals -O3 -static
else
	CWARNING_FLAGS=-Wall
	CFLAGS=$(CWARNING_FLAGS) -O3 -lm -lstdc++
	CDEBUG_FLAGS=$(CWARNING_FLAGS) -g -O0 -lm -lstdc++
endif

TARGETS=estimation \
estimation_debug \
estimation_sim \
estimation_sim_debug \
estimation_test \
estimation_sim_full \
estimation_sim_fullrent \
estimation_sim_fullwage \
estimation_stdev  \
estimation_test_full \
estimation_test_fullwage \
estimation_test_fullrent \
estimation_test_full_random \
estimation_test_fullwage_random \
estimation_test_random \
estimation_ref \
estimation_edu

all: $(TARGETS)

estimation: estimation.cpp
	@echo Building $@...;
	$(CC) $(CFLAGS) -DINFO -o $@ $<;

estimation_debug: estimation.cpp
	@echo Building $@...;
	$(CC) $(CDEBUG_FLAGS) -DTRACE -DINFO -o $@ $<;

estimation_sim: estimation.cpp
	@echo Building $@...;
	$(CC) $(CFLAGS) -DSIMULATION -DTRACE -DINFO -o $@ $<;

estimation_sim_debug: estimation.cpp
	@echo Building $@...;
	$(CC) $(CDEBUG_FLAGS) -DSIMULATION -DTRACE -DINFO -o $@ $<;

estimation_sim_full: estimation.cpp
	@echo Building $@...;
	$(CC) $(CFLAGS) -DSIMULATION -DTRACE -DINFO -DFULL_TRACE -DFULL_TRACE_INDEX -o $@ $<;

estimation_sim_fullrent: estimation.cpp
	@echo Building $@...;
	$(CC) $(CFLAGS) -DSIMULATION -DTRACE -DINFO -DFULL_TRACE -DFULL_TRACE_RENT -o $@ $<;

estimation_sim_fullwage: estimation.cpp
	@echo Building $@...;
	$(CC) $(CFLAGS) -DSIMULATION -DTRACE -DINFO -DFULL_TRACE -DFULL_TRACE_WAGE -o $@ $<;

estimation_sim4: estimation_sim4.cpp
	@echo Building $@...;
	$(CC) $(CFLAGS) -DSIMULATION -DTRACE -DINFO -o $@ $<;

estimation_sim4_full: estimation_sim4.cpp
	@echo Building $@...;
	$(CC) $(CFLAGS) -DSIMULATION -DTRACE -DINFO -DFULL_TRACE -DFULL_TRACE_INDEX -o $@ $<;

estimation_sim4_fullrent: estimation_sim4.cpp
	@echo Building $@...;
	$(CC) $(CFLAGS) -DSIMULATION -DTRACE -DINFO -DFULL_TRACE -DFULL_TRACE_RENT -o $@ $<;

estimation_sim4_fullwage: estimation_sim4.cpp
	@echo Building $@...;
	$(CC) $(CFLAGS) -DSIMULATION -DTRACE -DINFO -DFULL_TRACE -DFULL_TRACE_WAGE -o $@ $<;

estimation_stdev: estimation.cpp
	@echo Building $@...;
	$(CC) $(CFLAGS) -DCALC_STDEV -DINFO -o $@ $<;

estimation_test: estimation.cpp
	@echo Building $@...;
	$(CC) $(CFLAGS) -DTRACE -DINFO -o $@ $<;

estimation_test_full: estimation.cpp
	@echo Building $@...;
	$(CC) $(CFLAGS) -DTRACE -DINFO -DFULL_TRACE -DFULL_TRACE_INDEX -o $@ $<;

estimation_test_fullwage: estimation.cpp
	@echo Building $@...;
	$(CC) $(CFLAGS) -DTRACE -DINFO -DFULL_TRACE -DFULL_TRACE_WAGE -o $@ $<;

estimation_test_full_random: estimation.cpp
	@echo Building $@...;
	$(CC) $(CFLAGS) -DTRACE -DINFO -DFULL_TRACE -DFULL_TRACE_INDEX -DRANDOM_SELECTION -o $@ $<;

estimation_test_fullwage_random: estimation.cpp
	@echo Building $@...;
	$(CC) $(CFLAGS) -DTRACE -DINFO -DFULL_TRACE -DFULL_TRACE_WAGE -DRANDOM_SELECTION -o $@ $<;

estimation_test_random: estimation.cpp
	@echo Building $@...;
	$(CC) $(CFLAGS) -DTRACE -DINFO -DRANDOM_SELECTION -o $@ $<;

estimation_test_fullrent: estimation.cpp
	@echo Building $@...;
	$(CC) $(CFLAGS) -DTRACE -DINFO -DFULL_TRACE -DFULL_TRACE_RENT -o $@ $<;

estimation_ref: estimation.cpp
	@echo Building $@...;
	$(CC) $(CFLAGS) -DTRACE -DINFO -DREF_PARAM -o $@ $<;

estimation_edu: estimation.cpp
	@echo Building $@...;
	$(CC) $(CFLAGS) -DPRINT_EDU_LEVEL -DTRACE -DINFO -o $@ $<;
clean:
	@echo Removing all executables...
	@rm -rf $(TARGETS)

