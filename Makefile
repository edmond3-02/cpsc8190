
OFILES = base/Matrix.o \
	 base/Volume.o \
	 base/Fields.o \
	 base/LinearAlgebra.o \
	 base/Camera.o \
	 base/Color.o \
	 base/ImplicitVolumeShapes.o \
	 base/ImplicitVectorShapes.o \
	 base/ImplicitColors.o \
	 base/ImplicitMatrixShapes.o \
	 base/Stamp.o \
	 base/SparseGrid.o \
	 base/RayMarcher.o \
	 base/ObjParser.o \
	 base/VolumeGeometry.o

AFILES = $(OFILES)

ROOTDIR = .
LIB = -L$(ROOTDIR)/lib -lstarter  -lm -lOpenImageIO

CXX = g++ -shared -g -O2 -fPIC -fopenmp -fopenmp -std=c++14

SWIGCXX = g++ -shared -g -O2 -fPIC -fopenmp -fopenmp -std=c++14

PYTHONINCLUDE = -I/usr/include/python3.8

SWIGEXEC = swig

INCLUDES = -I ./include/ $(PYTHONINCLUDE)  

all: $(AFILES) 
	ar rv ./lib/libstarter.a $?

.C.o: $<
	$(CXX) -c $(INCLUDES) $< -o $@

clean:
	rm -rf *.o bin/* base/*.o core ./lib/libstarter.a  *~ swig/*~ swig/*.so swig/*.o swig/*.cxx swig/*.pyc swig/bishop.py* doc/html doc/latex python/*.pyc 

genswig:	swig/bishop.i	$(OFILES)
	$(SWIGEXEC) -c++ -python -shadow -I./include/ swig/bishop.i
	$(SWIGCXX) -c swig/bishop_wrap.cxx  $(INCLUDES) -o swig/bishop_wrap.o
	$(SWIGCXX) swig/bishop_wrap.o $(LIB) -o swig/_bishop.so

render:
	g++ -g -O2 -fPIC -fopenmp -std=c++14 base/main.C $(INCLUDES) -L. -lOpenImageIO lib/libstarter.a -o bin/main
	bin/main
	display images/grid_00001.exr
