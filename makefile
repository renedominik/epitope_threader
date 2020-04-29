OBJ = \
source/io.o \
source/centroid.o \
source/threader.o \
source/string_functions.o \
source/manual.o \
source/score_matrix.o \
source/molecule.o \
source/forcefield.o \
source/knowledge_based_potential_calculator.o \
source/vector3N.o \
source/vector.o \
source/vector_functions.o \
source/distribution.o \
epitope_threader.o

COMPILER = g++

NAME = exe/epitope_threader.exe

CPPFLAGS = -Wall -O1 -fmessage-length=0 -Wno-deprecated

.o:
	$(COMPILER) $(CPPFLAGS) -c -o $@ $< 

all: $(NAME)

$(NAME): $(OBJ)
	$(COMPILER) -o $(NAME) $(OBJ)

clean:
#	rm -f $(NAME)
	rm -f $(OBJ)

