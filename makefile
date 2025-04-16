TARGET_EXEC := sitnikov
MKDIR_COMM := mkdir -p
RM_COMM := rm -r
IGNORE_OUTPUT := >/dev/null 2>&1
ifeq ($(OS),Windows_NT)
    TARGET_EXEC := sitnikov.exe
    MKDIR_COMM := powershell mkdir -Force
    RM_COMM := powershell rmdir
    IGNORE_OUTPUT := | powershell Out-Null
endif
BUILD_DIR := ./build
SRC_DIR := ./src
SRCS := $(wildcard $(SRC_DIR)/*.c)
OBJS := $(patsubst %.c,%.o,$(foreach f,$(SRCS),$(subst $(SRC_DIR),$(BUILD_DIR),$(f))))
CC := gcc
CFLAGS := -O2 -I $(SRC_DIR)
LDFLAGS := -lm

$(TARGET_EXEC): $(OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	$(MKDIR_COMM) $(BUILD_DIR) $(IGNORE_OUTPUT)
	$(CC) $(CFLAGS) -c $< -o $@
.PHONY: clean
clean:
	$(RM_COMM) $(BUILD_DIR)
