BIN_DEBUG = ./bin/filterblastout-debug
BIN_RELEASE = ./bin/filterblastout
OBJ_DEBUG = ./obj_debug
OBJ_RELEASE = ./obj_release
SOURCE = src

GCC = g++

CPP_FILES := $(wildcard $(SOURCE)/*/*.cpp) $(wildcard $(SOURCE)/*.cpp)
CC_FILES := $(wildcard $(SOURCE)/*/*.cc) $(wildcard $(SOURCE)/*.cc)
H_FILES := $(wildcard $(SOURCE)/*/*.h) $(wildcard $(SOURCE)/*.h)

OBJ_FILES := $(CPP_FILES:.cpp=.o) $(CC_FILES:.cc=.o)
# OBJ_FILES_FOLDER := $(addprefix $(OBJ)/,$(OBJ_FILES))
OBJ_FILES_FOLDER_DEBUG := $(addprefix $(OBJ_DEBUG)/,$(OBJ_FILES))
OBJ_FILES_FOLDER_RELEASE := $(addprefix $(OBJ_RELEASE)/,$(OBJ_FILES))
# OBJ_FILES_FOLDER_MAC := $(addprefix $(OBJ_MAC)/,$(OBJ_FILES))

LIB_DIRS = -L"/usr/local/lib"
CC_LIBS = -static-libgcc -static-libstdc++ -D__cplusplus=201103L
INCLUDE = -I"./src/" -I"/usr/include/"

CC_FLAGS_DEBUG = -O0 -g -rdynamic -c -fmessage-length=0 -ffreestanding -fopenmp -m64 -std=c++11 -Werror=return-type -pthread
CC_FLAGS_RELEASE = -DRELEASE_VERSION -O3 -fdata-sections -ffunction-sections -c -fmessage-length=0 -ffreestanding -fopenmp -m64 -std=c++11 -Werror=return-type -pthread
LD_FLAGS = -static-libgcc -static-libstdc++ -m64 -ffreestanding
LD_LIBS = -lpthread -lgomp -lm -lz



all: release

release: $(OBJ_FILES_FOLDER_RELEASE)
	mkdir -p $(dir $(BIN_RELEASE))
	$(GCC) $(LD_FLAGS) $(LIB_DIRS) -o $(BIN_RELEASE) $(OBJ_FILES_FOLDER_RELEASE) $(LD_LIBS)
	
obj_release/%.o: %.cc $(H_FILES)
	mkdir -p $(dir $@)
	$(GCC) $(CC_LIBS) $(INCLUDE) $(CC_FLAGS_RELEASE) -o $@ $<
	
obj_release/%.o: %.cpp $(H_FILES)
	mkdir -p $(dir $@)
	$(GCC) $(CC_LIBS) $(INCLUDE) $(CC_FLAGS_RELEASE) -o $@ $<

debug: $(OBJ_FILES_FOLDER_DEBUG)
	mkdir -p $(dir $(BIN_DEBUG))
	$(GCC) $(LD_FLAGS) $(LIB_DIRS) -o $(BIN_DEBUG) $(OBJ_FILES_FOLDER_DEBUG) $(LD_LIBS)
	
obj_debug/%.o: %.cc $(H_FILES)
	mkdir -p $(dir $@)
	$(GCC) $(CC_LIBS) $(INCLUDE) $(CC_FLAGS_DEBUG) -o $@ $<
	
obj_debug/%.o: %.cpp $(H_FILES)
	mkdir -p $(dir $@)
	$(GCC) $(CC_LIBS) $(INCLUDE) $(CC_FLAGS_DEBUG) -o $@ $<
	
clean:
	-rm -rf $(OBJ_RELEASE) $(BIN_RELEASE)

cleandebug:
	-rm -rf $(OBJ_DEBUG) $(BIN_DEBUG)

rebuild: clean all
