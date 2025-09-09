
#eriksson:  eriksson.cpp
#	g++ -Wno-deprecated -O3 -o eriksson eriksson.cpp
#debug: eg.cpp
#	g++ -g -o eg eg.cpp
#classify:
#	g++ -o classify2 CoxeterSystem.cpp CoxeterElement.cpp classify.cpp

# to load the liberiksson.so file on unix, must set LD_LIBRARY_PATH to include current directory.
#   in csh this is "setenv LD_LIBRARY_PATH ."


unix:  liberiksson.so classify

cygwin:  liberiksson.dll classify.exe

# Note:  shared libraries under cygwin are named .dll (not .so).
liberiksson.dll:  CoxeterSystem.o CoxeterElement.o
	g++ -shared -o liberiksson.dll CoxeterSystem.o CoxeterElement.o 

liberiksson.so:  CoxeterSystem.o CoxeterElement.o
	g++ -shared -o liberiksson.so CoxeterSystem.o CoxeterElement.o 

#  the -c flag means do not link (just make the .o files from source).
CoxeterSystem.o:  CoxeterSystem.h CoxeterSystem.cpp
	g++ -Wno-deprecated -O3 -c -fPIC -o CoxeterSystem.o CoxeterSystem.cpp

CoxeterElement.o:  CoxeterElement.h CoxeterElement.cpp
	g++ -Wno-deprecated -O3 -c -fPIC -o CoxeterElement.o CoxeterElement.cpp

classify.exe:  classify.cpp
	g++ -Wno-deprecated -O3 -o classify.exe -L. -leriksson classify.cpp

classify:  classify.cpp
	g++ -Wno-deprecated -O3 -o classify -L. -leriksson classify.cpp

clean:
	rm *.exe *.o *.dll *.so classify

