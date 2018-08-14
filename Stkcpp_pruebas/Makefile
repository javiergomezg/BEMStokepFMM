cc = g++-6 
BINARY = Stokes
OBJS = Stokes.o

SOURCES=src/H2_2D_Tree.cpp src/H2_2D_Node.cpp src/kernel_Base.cpp src/kernel_Types.cpp 
OBJECTSA=$(SOURCES:.cpp=.o) $(SOURCESA:.cpp=.o)

# Con esta regla se compila el programa Stokes
Stokes: $(OBJS) $(OBJECTSA)
	$(cc) $(OBJS) $(OBJECTSA) -o $(BINARY)	

Stokes.o: Stokes.cpp cabecera.hpp
	$(cc) -c Stokes.cpp $(SOURCES)

all: Stokes

clean: 
	rm -f $(BINARY) $(OBJS) $(OBJECTSA) 
