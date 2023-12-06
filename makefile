# Compiler and flags
CXX = clang++
CXXFLAGS = -Wall -Wextra -std=c++11
DEBUGFLAGS = -g

# Directories
SRC_DIR = src
INC_DIR = include
BUILD_DIR = build

# Source and object files
SRC_FILES := $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SRC_FILES))

# Executable name
TARGET = main.exe
TARGET_DEBUG = main_debug.exe

# Default target
all: release

# Build release executable
release: $(OBJ_FILES)
	$(CXX) $(CXXFLAGS) $^ -O3 -o $(TARGET)

# Build debug executable
debug: CXXFLAGS += $(DEBUGFLAGS)
debug: $(OBJ_FILES)
	$(CXX) $(CXXFLAGS) $^ -o $(TARGET_DEBUG) -Wall -address -fsanitize=memory 
	
# -address -fsanitize=memory 
# this is for memory leak detection

# Compile source files to object files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -I$(INC_DIR) -c $< -o $@

# Ensure build directory exists before compiling
$(shell mkdir -p $(BUILD_DIR))

# Clean up
clean:
	rm -rf $(BUILD_DIR) $(TARGET) $(TARGET_DEBUG)

p:
	python3 plot.py

cleandata:
	rm ./data/*.txt ./plot/*.png ./plot/*.gif