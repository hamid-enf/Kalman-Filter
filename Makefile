# ============================================================================
# Kalman Filter library — build system
#
# Targets:
#   make             build the static library and all host tools
#   make lib         build build/libkalman.a
#   make test        build and run the unit test suite
#   make examples    build all host examples
#   make run-examples  build and run all host examples
#   make bench       build and run the benchmark
#   make clean       remove all build artefacts
#
# Configuration is via kalman/include/kalman_config.h (or -D flags). For a
# minimal float-only linear-KF build, for example:
#   make CFLAGS="-std=c11 -O2 -DKF_ENABLE_EKF=0 -DKF_ENABLE_UKF=0"
# ============================================================================

CC      ?= gcc
AR      ?= ar

CFLAGS  ?= -std=c11 -O2 -Wall -Wextra -Werror
# Include paths live in a separate variable so that a command-line CFLAGS
# override (e.g. `make CFLAGS="-O3 ..."`) does not drop the library headers.
INC     := -Ikalman/include
LDLIBS  += -lm

BUILD   := build

CORE_SRC := \
    kalman/src/kalman_matrix.c \
    kalman/src/kalman_kf.c \
    kalman/src/kalman_ekf.c \
    kalman/src/kalman_ukf.c \
    kalman/src/kalman_common.c

CORE_OBJ := $(CORE_SRC:%.c=$(BUILD)/%.o)
LIB      := $(BUILD)/libkalman.a

TEST_SRC := \
    tests/test_framework.c \
    tests/test_matrix.c \
    tests/test_kf.c \
    tests/test_ekf.c \
    tests/test_ukf.c \
    tests/test_extensions.c \
    tests/test_reference.c \
    tests/test_boundary.c \
    tests/test_stress.c \
    tests/main.c

EXAMPLES := \
    examples/01_basic_1d/main.c \
    examples/02_position/main.c \
    examples/03_position_velocity/main.c \
    examples/04_position_velocity_acceleration/main.c \
    examples/05_multiple_measurements/main.c \
    examples/06_two_sensors/main.c \
    examples/07_predict_update_rates/main.c \
    examples/08_variable_dt/main.c \
    examples/09_ekf_basic/main.c \
    examples/10_ukf_basic/main.c \
    examples/11_outlier_gating/main.c \
    examples/12_adaptive_r/main.c \
    examples/13_rts_smoother/main.c

EXE := \
    $(BUILD)/examples/01_basic_1d \
    $(BUILD)/examples/02_position \
    $(BUILD)/examples/03_position_velocity \
    $(BUILD)/examples/04_position_velocity_acceleration \
    $(BUILD)/examples/05_multiple_measurements \
    $(BUILD)/examples/06_two_sensors \
    $(BUILD)/examples/07_predict_update_rates \
    $(BUILD)/examples/08_variable_dt \
    $(BUILD)/examples/09_ekf_basic \
    $(BUILD)/examples/10_ukf_basic \
    $(BUILD)/examples/11_outlier_gating \
    $(BUILD)/examples/12_adaptive_r \
    $(BUILD)/examples/13_rts_smoother

.PHONY: all lib test examples run-examples bench check sanitize clean

all: lib examples test bench

# ---------------------------------------------------------------------------
# Library
# ---------------------------------------------------------------------------

lib: $(LIB)

$(BUILD)/kalman/src/%.o: kalman/src/%.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

$(LIB): $(CORE_OBJ)
	@mkdir -p $(BUILD)
	$(AR) rcs $@ $^

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

test: $(BUILD)/tests/test_suite
	@$(BUILD)/tests/test_suite

$(BUILD)/tests/test_suite: $(TEST_SRC) $(LIB)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(INC) -Itests $(TEST_SRC) $(LIB) $(LDLIBS) -o $@

# ---------------------------------------------------------------------------
# Examples
# ---------------------------------------------------------------------------

examples: $(EXE)

$(BUILD)/examples/%: examples/%/main.c $(LIB)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(INC) -Iexamples/common $< $(LIB) $(LDLIBS) -o $@

# ---------------------------------------------------------------------------
# Benchmark
# ---------------------------------------------------------------------------

bench: $(BUILD)/benchmarks/benchmark
	@$(BUILD)/benchmarks/benchmark

$(BUILD)/benchmarks/benchmark: benchmarks/benchmark.c $(LIB)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(INC) -Ibenchmarks $< $(LIB) $(LDLIBS) -o $@

# ---------------------------------------------------------------------------
# Static analysis / MISRA
# ---------------------------------------------------------------------------

# Runs cppcheck (if installed) and the MISRA-C:2012 addon (see tools/).
check:
	@if command -v cppcheck >/dev/null 2>&1; then \
		cppcheck --std=c11 --enable=warning,style,performance,portability \
			--error-exitcode=1 -Ikalman/include kalman/src/ || true; \
	else \
		echo "cppcheck not installed; skipping general static analysis"; \
	fi
	@if command -v python3 >/dev/null 2>&1; then \
		CPPCHECK_ADDON="$${CPPCHECK_ADDON:-}" ./tools/misra_check.sh || true; \
	else \
		echo "python3 not installed; skipping MISRA check"; \
	fi

# ---------------------------------------------------------------------------
# Dynamic analysis
# ---------------------------------------------------------------------------

# Runs the suite under ASan + UBSan (memory safety / undefined behaviour).
sanitize:
	@./tools/sanitize.sh

# ---------------------------------------------------------------------------
# Convenience
# ---------------------------------------------------------------------------

run-examples: examples
	@for e in $(EXE); do \
		if echo $$e | grep -q examples; then \
			echo "=== $$e ==="; $$e || exit 1; echo; \
		fi; \
	done

clean:
	rm -rf $(BUILD)
