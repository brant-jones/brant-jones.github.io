#include <iostream.h>
#include <math.h>

#include <list.h>
#include <vector.h>
#include <map.h>

#include <string>
#include <sstream>
using namespace std;

//////////////////////////////////////////////////////////////////////
//
// This is an implementation of the set of 2^{length} masks used
// in Deodhar's algorithm to obtain subexpressions.
//
// The masks are stored as (length/MASK_TYPESIZE) copies of a
// MASK_DATATYPE.
//
//////////////////////////////////////////////////////////////////////

//#define MASK_DATATYPE unsigned short
//#define MASK_TYPESIZE 16.0

#define MASK_DATATYPE unsigned int
#define MASK_TYPESIZE 32.0
#define MASK_FLIP_AT -1U


class Masks
{
        public:
                int length;
                int size;
                MASK_DATATYPE* iterator;
                int pointer;
                Masks::Masks(int length);
                Masks::~Masks();
                int Masks::exhausted();
                int Masks::next();
                int Masks::get_value(int position);
                int Masks::proper();
                void Masks::print();
		void Masks::sprint(string& s);
                void Masks::test();
                int exhausted_flag;
};

