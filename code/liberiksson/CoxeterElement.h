#include "CoxeterSystem.h"


//////////////////////////////////////////////////////////////////////
//
// This is an implementation of the numbers game by Kimmo Eriksson, 
// as described by Bjorner/Brenti in Combinatorics of Coxeter Groups.
//
//////////////////////////////////////////////////////////////////////


class CoxeterElement
{
	public:
		CoxeterSystem* coxeter_system;  // coxeter_system is not allocated in this class.
		int size;  // number of generators in the coxeter matrix from coxeter_system->size.

		int* word;  // this is generalized 1-line notation from numbers game...
		int* one_line;  // this is the usual 1-line notation (only guarenteed to make sense in types A and D)...
		int length; // WARNING:  users should not update word or one_line directly, without adjusting length.

		CoxeterElement::CoxeterElement(CoxeterSystem* cs);  // construct identity element.
		CoxeterElement::CoxeterElement(CoxeterSystem* cs, int w[], int ol[]);  // construct element with given word and 1-line notation.
		CoxeterElement::CoxeterElement(CoxeterSystem* cs, int red[], int len);  // construct element with given reduced word.
		CoxeterElement::CoxeterElement(const CoxeterElement& cp);   // copy-constructor: for passing by reference in functions and creating copies on the heap.
		//CoxeterElement::operator=(const CoxeterElement& cp);   // = operator copy-constructor: for passing by reference in functions and creating copies on the heap.
		CoxeterElement::~CoxeterElement();

		int CoxeterElement::equals(CoxeterElement* cp);

		void CoxeterElement::print();
		void CoxeterElement::print_reduced_expression();

		void CoxeterElement::right_multiply(int s);
		void CoxeterElement::left_multiply(int s);  // WARNING:  this is not speed-optimized...
		int CoxeterElement::get_length();
		void CoxeterElement::get_reduced_expression(int reduced[]);  // requires an allocated int array of length at least this->length.
		void CoxeterElement::get_reduced_expression(int reduced[], int& count_moves);
		int CoxeterElement::get_rank(int reduced[]);

		void CoxeterElement::mask(CoxeterElement* r, long long m1, long long m2, int reduced[]);
		void CoxeterElement::mask(CoxeterElement* r, long long m1, long long m2, int imax, int reduced[]);
		int CoxeterElement::defect(long long m1, long long m2, int reduced[]);
		int CoxeterElement::deodhar();

		int CoxeterElement::contains_one_line_pattern(CoxeterElement& pattern);

		void CoxeterElement::print_heap();

};


