
//////////////////////////////////////////////////////////////////////
//
// This is an implementation of the numbers game by Kimmo Eriksson, 
// as described by Bjorner/Brenti in Combinatorics of Coxeter Groups.
//
//////////////////////////////////////////////////////////////////////

#include <iostream.h>
#include <math.h>

#include <list.h>
#include <vector.h>
#include <map.h>

#include <string>
#include <sstream>
using namespace std;

const int DEBUG_VERBOSE = 0;
const int DEBUG_VERBOSE_CP = 0;

const int DATATYPE_BITLENGTH = 64;  // corresponds to long long, used for calculating subsets via masks.  This should be long enough to calculate through D_11...

class CoxeterSystem
{
	public:
		int size;
		int *coxeter_matrix;

		int automorphism_group_size;
		int *automorphism_group;
		
		int CoxeterSystem::equals(CoxeterSystem* cs);

		CoxeterSystem::CoxeterSystem(int n, int cm[], int ags, int ag[]);
		CoxeterSystem::CoxeterSystem(const CoxeterSystem& cs);
		int CoxeterSystem::get_exponent(int x, int y);
		int CoxeterSystem::get_automorphism(int i, int j);
		void CoxeterSystem::print_matrix();
		virtual int right_multiply(int one_line[], int i);
		CoxeterSystem::~CoxeterSystem();
};

class TypeDCoxeterSystem : public CoxeterSystem
{
	public:
		TypeDCoxeterSystem(int n, int cm[], int ags, int ag[]) : CoxeterSystem(n, cm, ags, ag) {}
		TypeDCoxeterSystem(const CoxeterSystem& cs) : CoxeterSystem(cs) {}
		virtual int right_multiply(int one_line[], int i);
};


//////////////////////////////////////////////////////////////////////
// 
// Coxeter matricies for various groups:
//
//////////////////////////////////////////////////////////////////////
	

	// Type A:  Linear...

	// A_2:  *--*
	//       0  1
	static int aa2[] = { 0,3,
			     3,0 };
	static int ga2[] = { 1, 0 };
	static CoxeterSystem A2 = CoxeterSystem(2, aa2, 1, ga2);
	
	// A_3:  *--*--*
	//       0  1  2
	static int aa3[] = { 0,3,2,
			     3,0,3,
			     2,3,0 };
	static int ga3[] = { 2, 1, 0 };
	static CoxeterSystem A3 = CoxeterSystem(3, aa3, 1, ga3);
	
	// A_4:  *--*--*--*
	//       0  1  2  3
	static int aa4[] = { 0,3,2,2,
			     3,0,3,2,
			     2,3,0,3,
			     2,2,3,0 };
	static int ga4[] = { 3, 2, 1, 0 };
	static CoxeterSystem A4 = CoxeterSystem(4, aa4, 1, ga4);

	static int aa5[] = { 0,3,2,2,2,
			     3,0,3,2,2,
			     2,3,0,3,2,
			     2,2,3,0,3,
			     2,2,2,3,0 };
	static int ga5[] = { 4, 3, 2, 1, 0 };
	static CoxeterSystem A5 = CoxeterSystem(5, aa5, 1, ga5);

	static int aa6[] = { 0,3,2,2,2,2,
			     3,0,3,2,2,2,
			     2,3,0,3,2,2,
			     2,2,3,0,3,2,
			     2,2,2,3,0,3,
			     2,2,2,2,3,0 };
	static int ga6[] = { 5, 4, 3, 2, 1, 0 };
	static CoxeterSystem A6 = CoxeterSystem(6, aa6, 1, ga6);

	// A_7:  *--*--*--*--*--*--*
	//       0  1  2  3  4  5  6
	static int aa7[] = { 0,3,2,2,2,2,2,
			     3,0,3,2,2,2,2,
			     2,3,0,3,2,2,2,
			     2,2,3,0,3,2,2,
			     2,2,2,3,0,3,2,
			     2,2,2,2,3,0,3,
			     2,2,2,2,2,3,0 };
	static int ga7[] = { 6, 5, 4, 3, 2, 1, 0 };
	static CoxeterSystem A7 = CoxeterSystem(7, aa7, 1, ga7);

	static int aa8[] = { 0,3,2,2,2,2,2,2,
			     3,0,3,2,2,2,2,2,
			     2,3,0,3,2,2,2,2,
			     2,2,3,0,3,2,2,2,
			     2,2,2,3,0,3,2,2,
			     2,2,2,2,3,0,3,2,
			     2,2,2,2,2,3,0,3,
			     2,2,2,2,2,2,3,0 };
	static int ga8[] = { 7, 6, 5, 4, 3, 2, 1, 0 };
	static CoxeterSystem A8 = CoxeterSystem(8, aa8, 1, ga8);

	static int aa9[] = { 0,3,2,2,2,2,2,2,2,
			     3,0,3,2,2,2,2,2,2,
			     2,3,0,3,2,2,2,2,2,
			     2,2,3,0,3,2,2,2,2,
			     2,2,2,3,0,3,2,2,2,
			     2,2,2,2,3,0,3,2,2,
			     2,2,2,2,2,3,0,3,2,
			     2,2,2,2,2,2,3,0,3,
			     2,2,2,2,2,2,2,3,0 };
	static int ga9[] = { 8, 7, 6, 5, 4, 3, 2, 1, 0 };
	static CoxeterSystem A9 = CoxeterSystem(9, aa9, 1, ga9);

	static int aa10[] = { 0,3,2,2,2,2,2,2,2,2,
			      3,0,3,2,2,2,2,2,2,2,
			      2,3,0,3,2,2,2,2,2,2,
			      2,2,3,0,3,2,2,2,2,2,
			      2,2,2,3,0,3,2,2,2,2,
			      2,2,2,2,3,0,3,2,2,2,
			      2,2,2,2,2,3,0,3,2,2,
			      2,2,2,2,2,2,3,0,3,2,
			      2,2,2,2,2,2,2,3,0,3,
			      2,2,2,2,2,2,2,2,3,0 };
	static int ga10[] = { 9, 8, 7, 6, 5, 4, 3, 2, 1, 0 };
	static CoxeterSystem A10 = CoxeterSystem(10, aa10, 1, ga10);

        static int aa11[] = { 0,3,2,2,2,2,2,2,2,2,2,
                              3,0,3,2,2,2,2,2,2,2,2,
                              2,3,0,3,2,2,2,2,2,2,2,
                              2,2,3,0,3,2,2,2,2,2,2,
                              2,2,2,3,0,3,2,2,2,2,2,
                              2,2,2,2,3,0,3,2,2,2,2,
                              2,2,2,2,2,3,0,3,2,2,2,
                              2,2,2,2,2,2,3,0,3,2,2,
                              2,2,2,2,2,2,2,3,0,3,2,
                              2,2,2,2,2,2,2,2,3,0,3,
                              2,2,2,2,2,2,2,2,2,3,0 };
        static int ga11[] = { 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0 };
        static CoxeterSystem A11 = CoxeterSystem(11, aa11, 1, ga11);

        static int aa12[] = { 0,3,2,2,2,2,2,2,2,2,2,2,
                              3,0,3,2,2,2,2,2,2,2,2,2,
                              2,3,0,3,2,2,2,2,2,2,2,2,
                              2,2,3,0,3,2,2,2,2,2,2,2,
                              2,2,2,3,0,3,2,2,2,2,2,2,
                              2,2,2,2,3,0,3,2,2,2,2,2,
                              2,2,2,2,2,3,0,3,2,2,2,2,
                              2,2,2,2,2,2,3,0,3,2,2,2,
                              2,2,2,2,2,2,2,3,0,3,2,2,
                              2,2,2,2,2,2,2,2,3,0,3,2,
                              2,2,2,2,2,2,2,2,2,3,0,3,
                              2,2,2,2,2,2,2,2,2,2,3,0 };
        static int ga12[] = { 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0 };
        static CoxeterSystem A12 = CoxeterSystem(12, aa12, 1, ga12);

	// Type B:  Linear, 4-edge on the left...
	
        // B_3:  *-4-*--*--*--*--*
        //       0   1  2  3  4  4
        static int ab3[] = { 0,4,2,
                             4,0,3,
                             2,3,0 };
        static CoxeterSystem B3 = CoxeterSystem(3, ab3, 0, NULL);

        // B_4:  *-4-*--*--*--*--*
        //       0   1  2  3  4  4
        static int ab4[] = { 0,4,2,2,
                             4,0,3,2,
                             2,3,0,3,
                             2,2,3,0 };
        static CoxeterSystem B4 = CoxeterSystem(4, ab4, 0, NULL);

        // B_5:  *-4-*--*--*--*--*
        //       0   1  2  3  4  5  
        static int ab5[] = { 0,4,2,2,2,
                             4,0,3,2,2,
                             2,3,0,3,2,
                             2,2,3,0,3,
                             2,2,2,3,0 };
        static CoxeterSystem B5 = CoxeterSystem(5, ab5, 0, NULL);

        // B_6:  *-4-*--*--*--*--*
        //       0   1  2  3  4  5
        static int ab6[] = { 0,4,2,2,2,2,
                             4,0,3,2,2,2,
                             2,3,0,3,2,2,
                             2,2,3,0,3,2,
                             2,2,2,3,0,3,
                             2,2,2,2,3,0 };
        static CoxeterSystem B6 = CoxeterSystem(6, ab6, 0, NULL);

	// B_7:  *-4-*--*--*--*--*--*
	//       0   1  2  3  4  5  6
	static int ab7[] = { 0,4,2,2,2,2,2,
			     4,0,3,2,2,2,2,
			     2,3,0,3,2,2,2,
			     2,2,3,0,3,2,2,
			     2,2,2,3,0,3,2,
			     2,2,2,2,3,0,3,
			     2,2,2,2,2,3,0 };
	static CoxeterSystem B7 = CoxeterSystem(7, ab7, 0, NULL);
	
	// F_4:  *---*-4-*--*
	//       0   1   2  3
	static int af4[] = { 0,3,2,2,
			     3,0,4,2,
			     2,4,0,3,
			     2,2,3,0 };
	static CoxeterSystem F4 = CoxeterSystem(4, af4, 0, NULL);

	// G_2:  *-6-*
	//       0   1
	static int ag2[] = { 0,6,
			     6,0 };
	static CoxeterSystem G2 = CoxeterSystem(2, ag2, 0, NULL);

	// Type D:  Branch on the left, 0, 1 are branch points connected to 2, then linear...

	// D_8:  1
	//       *---
	//       *---*--*--*--*--*--*
	//       0   2  3  4  5  6  7
	
	static int ad3[] = { 0,2,3,
			     2,0,3,
			     3,3,0 };
	static int gd3[] = { 1, 0, 2 };
	static TypeDCoxeterSystem D3 = TypeDCoxeterSystem(3, ad3, 1, gd3);

	static int ad4[] = { 0,2,3,2,
			     2,0,3,2,
			     3,3,0,3,
			     2,2,3,0 };
	static int gd4[] = { 1, 0, 2, 3 };
	static TypeDCoxeterSystem D4 = TypeDCoxeterSystem(4, ad4, 1, gd4);

	static int ad5[] = { 0,2,3,2,2,
			     2,0,3,2,2,
			     3,3,0,3,2,
			     2,2,3,0,3,
			     2,2,2,3,0 };
	static int gd5[] = { 1, 0, 2, 3, 4 };
	static TypeDCoxeterSystem D5 = TypeDCoxeterSystem(5, ad5, 1, gd5);

	// D_6:  0
	//       *--
	//       *--*--*--*--*
	//       1  2  3  4  5
	//
	// Automorphisms:  interchange 0<->1, fix all others.
	static int ad6[] = { 0,2,3,2,2,2,
			     2,0,3,2,2,2,
			     3,3,0,3,2,2,
			     2,2,3,0,3,2,
			     2,2,2,3,0,3,
			     2,2,2,2,3,0 };
	static int gd6[] = { 1, 0, 2, 3, 4, 5 };
	static TypeDCoxeterSystem D6 = TypeDCoxeterSystem(6, ad6, 1, gd6);

	// D_7:  0
	//       *--
	//       *--*--*--*--*--*
	//       1  2  3  4  5  6
	//
	// Automorphisms:  interchange 0<->1, fix all others.
	// Parabolic subgroups:  must use 6 (or get D_6).
	static int ad7[] = { 0,2,3,2,2,2,2,
			     2,0,3,2,2,2,2,
			     3,3,0,3,2,2,2,
			     2,2,3,0,3,2,2,
			     2,2,2,3,0,3,2,
			     2,2,2,2,3,0,3,
			     2,2,2,2,2,3,0 };
	static int gd7[] = { 1, 0, 2, 3, 4, 5, 6 };
	static TypeDCoxeterSystem D7 = TypeDCoxeterSystem(7, ad7, 1, gd7);

	// D_8:  0
	//       *--
	//       *--*--*--*--*--*--*
	//       1  2  3  4  5  6  7
	//
	// Automorphisms:  interchange 0<->1, fix all others.
	// Parabolic subgroups:  must use 0 (or get A_7), must use 7 (or get D_7), must use 6/7 (or get D_6).
	static int ad8[] = { 0,2,3,2,2,2,2,2,
			     2,0,3,2,2,2,2,2,
			     3,3,0,3,2,2,2,2,
			     2,2,3,0,3,2,2,2,
			     2,2,2,3,0,3,2,2,
			     2,2,2,2,3,0,3,2,
			     2,2,2,2,2,3,0,3,
			     2,2,2,2,2,2,3,0};
	static int gd8[] = { 1, 0, 2, 3, 4, 5, 6, 7 };
	static TypeDCoxeterSystem D8 = TypeDCoxeterSystem(8, ad8, 1, gd8);

	static int ad9[] = { 0,2,3,2,2,2,2,2,2,
			     2,0,3,2,2,2,2,2,2,
			     3,3,0,3,2,2,2,2,2,
			     2,2,3,0,3,2,2,2,2,
			     2,2,2,3,0,3,2,2,2,
			     2,2,2,2,3,0,3,2,2,
			     2,2,2,2,2,3,0,3,2,
			     2,2,2,2,2,2,3,0,3,
			     2,2,2,2,2,2,2,3,0 };
	static int gd9[] = { 1, 0, 2, 3, 4, 5, 6, 7, 8 };
	static TypeDCoxeterSystem D9 = TypeDCoxeterSystem(9, ad9, 1, gd9);

	static int ad10[] = { 0,2,3,2,2,2,2,2,2,2,
			     2,0,3,2,2,2,2,2,2,2,
			     3,3,0,3,2,2,2,2,2,2,
			     2,2,3,0,3,2,2,2,2,2,
			     2,2,2,3,0,3,2,2,2,2,
			     2,2,2,2,3,0,3,2,2,2,
			     2,2,2,2,2,3,0,3,2,2,
			     2,2,2,2,2,2,3,0,3,2,
			     2,2,2,2,2,2,2,3,0,3,
			     2,2,2,2,2,2,2,2,3,0 };
	static int gd10[] = { 1, 0, 2, 3, 4, 5, 6, 7, 8, 9 };
	static TypeDCoxeterSystem D10 = TypeDCoxeterSystem(10, ad10, 1, gd10);

	// E_6:        5
	//             *
	//       *--*--*--*--*
	//       0  1  2  3  4
	//
	// Automorphisms:  interchange 0<->4, 1<->3, fix 2 and 5.
	static int ae6[] = { 0,3,2,2,2,2,
			     3,0,3,2,2,2,
			     2,3,0,3,2,3,
			     2,2,3,0,3,2,
			     2,2,2,3,0,2,
			     2,2,3,2,2,0 };
	static int ge6[] = { 4, 3, 2, 1, 0, 5 };
	static CoxeterSystem E6 = CoxeterSystem(6, ae6, 1, ge6);

	// TEMPORARY:  to compare with coxeter 1.0.
	// E_6:        2
	//             *
	//       *--*--*--*--*
	//     0 1  3  4  5  6
	//
	// Automorphisms:  interchange 1<->6, 3<->5, fix 2 and 4.
	//static int ae6[] = { 0,2,2,2,2,2,2,
	//		     2,0,2,3,2,2,2,
	//		     2,2,0,2,3,2,2,
	//		     2,3,2,0,3,2,2,
	//		     2,2,3,3,0,3,2,
	//		     2,2,2,2,3,0,3,
	//		     2,2,2,2,2,3,0 };
	//static int ge6[] = { 0, 6, 2, 5, 4, 3, 1 };
	//static CoxeterSystem E6 = CoxeterSystem(7, ae6, 1, ge6);

	// E_7:        5
	//             *
	//       *--*--*--*--*--*
	//       0  1  2  3  4  6
	//
	// Parabolic subgroups:  must use 0 (or get D_6), 6 (or get E_6)
	static int ae7[] = { 0,3,2,2,2,2,2,
			     3,0,3,2,2,2,2,
			     2,3,0,3,2,3,2,
			     2,2,3,0,3,2,2,
			     2,2,2,3,0,2,3,
			     2,2,3,2,2,0,2,
			     2,2,2,2,3,2,0 };
	static CoxeterSystem E7 = CoxeterSystem(7, ae7, 0, NULL);

	// E_8:        5
	//             *
	//       *--*--*--*--*--*--*
	//       0  1  2  3  4  6  7
	//
	// Parabolic subgroups:  must use 0 (or get D_7), 6/7 (or get E_6), 7 (or get E_7), 5 (or get A_7)
	static int ae8[] = { 0,3,2,2,2,2,2,2,
			     3,0,3,2,2,2,2,2,
			     2,3,0,3,2,3,2,2,
			     2,2,3,0,3,2,2,2,
			     2,2,2,3,0,2,3,2,
			     2,2,3,2,2,0,2,2,
			     2,2,2,2,3,2,0,3,
			     2,2,2,2,2,2,3,0 };
	static CoxeterSystem E8 = CoxeterSystem(8, ae8, 0, NULL);

class CoxeterNames
{
	public:
        	map<string, CoxeterSystem*> names;

		CoxeterNames::CoxeterNames()
		{
        		names["A2"] = &A2;
		        names["A3"] = &A3;
		        names["A4"] = &A4;
		        names["A5"] = &A5;
		        names["A6"] = &A6;
		        names["A7"] = &A7;
		        names["A8"] = &A8;
		        names["A9"] = &A9;
		        names["A10"] = &A10;
		        names["A11"] = &A11;
		        names["A12"] = &A12;
		        names["B3"] = &B3;
		        names["B4"] = &B4;
		        names["B5"] = &B5;
		        names["B6"] = &B6;
		        names["B7"] = &B7;
		        names["D3"] = &D3;
		        names["D4"] = &D4;
		        names["D5"] = &D5;
		        names["D6"] = &D6;
		        names["D7"] = &D7;
		        names["D8"] = &D8;
		        names["D9"] = &D9;
		        names["D10"] = &D10;
		        names["E6"] = &E6;
		        names["E7"] = &E7;
		        names["E8"] = &E8;
		        names["F4"] = &F4;
		        names["G2"] = &G2;
		}
};


