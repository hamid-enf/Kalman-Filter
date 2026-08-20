# Host build for the stm32_h743 examples (uses the HAL mock via -DH743_HOST).
#
#   make            build all three examples
#   make run        build and run all three
#   make clean
#
# On real hardware these files are added to a CubeIDE project instead (see
# README.md); this Makefile exists only to let you preview the output on a PC.

CC     ?= gcc
CFLAGS ?= -std=c11 -O2 -Wall -Wextra -Werror -DH743_HOST -I. -I../../kalman/include

CORE := ../../kalman/src/kalman_matrix.c \
        ../../kalman/src/kalman_kf.c \
        ../../kalman/src/kalman_ekf.c \
        ../../kalman/src/kalman_ukf.c \
        ../../kalman/src/kalman_common.c

EXAMPLES := 01_simple_temperature 02_position_velocity 03_sensor_fusion

.PHONY: all run clean

all: $(EXAMPLES)

%: %.c common.h host_hal.h
	$(CC) $(CFLAGS) $< $(CORE) -lm -o $@

run: all
	@for ex in $(EXAMPLES); do echo; ./$$ex; done

clean:
	rm -f $(EXAMPLES)
