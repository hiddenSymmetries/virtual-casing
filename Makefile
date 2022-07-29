CC = gcc
CXX = g++
CXXFLAGS = -std=c++11 -fopenmp -march=native -Wall # need C++11 and OpenMP

#Optional flags
#CXXFLAGS += -O0 -fsanitize=address # debug build
CXXFLAGS += -O3 -DNDEBUG # release build

ifeq ($(shell uname -s),Darwin)
	CXXFLAGS += -g -rdynamic -Wl,-no_pie # for stack trace (on Mac)
else
	CXXFLAGS += -g -rdynamic # for stack trace
endif

#CXXFLAGS += -DSCTL_MEMDEBUG # Enable memory checks
CXXFLAGS += -DSCTL_PROFILE=5 -DSCTL_VERBOSE # Enable profiling

#CXXFLAGS += -lblas -DSCTL_HAVE_BLAS # use BLAS
#CXXFLAGS += -llapack -DSCTL_HAVE_LAPACK # use LAPACK
#CXXFLAGS += -qmkl -DSCTL_HAVE_BLAS -DSCTL_HAVE_LAPACK -DSCTL_HAVE_FFTW3_MKL # use MKL BLAS, LAPACK and FFTW (Intel compiler)
CXXFLAGS += -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -DSCTL_HAVE_BLAS -DSCTL_HAVE_LAPACK # use MKL BLAS and LAPACK (non-Intel compiler)
#CXXFLAGS += -DSCTL_HAVE_SVML

CXXFLAGS += -lfftw3_omp -DSCTL_FFTW_THREADS
CXXFLAGS += -lfftw3 -DSCTL_HAVE_FFTW
CXXFLAGS += -lfftw3f -DSCTL_HAVE_FFTWF
CXXFLAGS += -lfftw3l -DSCTL_HAVE_FFTWL

RM = rm -f
MKDIRS = mkdir -p

BINDIR = ./bin
LIBDIR = ./lib
SRCDIR = ./src
OBJDIR = ./obj
INCDIR = ./include

TARGET_LIB = $(LIBDIR)/libvirtualcasing.a
TARGET_BIN = $(BINDIR)/virtual-casing \
						 $(BINDIR)/virtual-casing-c

BIEST_INCDIR = ./BIEST/include

all : $(TARGET_BIN) $(TARGET_LIB)

$(BINDIR)/%: ./test/%.cpp
	-@$(MKDIRS) $(dir $@)
	$(CXX) -I$(INCDIR) -I$(BIEST_INCDIR) $^ -o $@ $(CXXFLAGS)

$(BINDIR)/%: ./test/%.c $(TARGET_LIB)
	-@$(MKDIRS) $(dir $@)
	$(CC) -I$(INCDIR) $^ $(TARGET_LIB) -lm -ldl -lstdc++ -o $@ $(CXXFLAGS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	-@$(MKDIRS) $(dir $@)
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -I$(BIEST_INCDIR) -c $^ -o $@

$(TARGET_LIB) : obj/virtual-casing.o
	-@$(MKDIRS) $(dir $@)
	ar rcs $@ $^

clean:
	$(RM) -r $(BINDIR)/* $(OBJDIR)/* $(LIBDIR)/*

