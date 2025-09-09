
# to load the liberiksson.so file on unix, must set LD_LIBRARY_PATH to include current directory.
#   in csh this is "setenv LD_LIBRARY_PATH ."


unix:  liberiksson.so deodhar verify

cygwin:  liberiksson.dll deodhar.exe verify.exe

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

verify.exe:  verify.cpp Masks.cpp
	g++ -Wno-deprecated -O3 -o verify.exe -L. -leriksson verify.cpp Masks.cpp

verify:  verify.cpp Masks.cpp
	g++ -Wno-deprecated -O3 -o verify -L. -leriksson verify.cpp Masks.cpp

deodhar:  deodhar.cpp Masks.cpp
	g++ -Wno-deprecated -O3 -o deodhar -L. -leriksson deodhar.cpp Masks.cpp

deodhar.exe:  deodhar.cpp Masks.cpp
	g++ -Wno-deprecated -O3 -o deodhar.exe -L. -leriksson deodhar.cpp Masks.cpp

clean:
	rm *.exe *.o *.dll *.so verify deodhar

