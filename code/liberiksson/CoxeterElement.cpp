#include "CoxeterElement.h"


	// construct identity element
	CoxeterElement::CoxeterElement(CoxeterSystem* cs)
	{
		coxeter_system = cs;
		size = coxeter_system->size;

		word = new int[size];
		one_line = new int[size+1];
		for (int i = 0; i < size; i++)
		{
			word[i] = 1;
			one_line[i] = i+1;
		}
		one_line[size] = size+1;

		length = 0;
	}

	// from 1-line.
	CoxeterElement::CoxeterElement(CoxeterSystem* cs, int w[], int ol[])
	{
		coxeter_system = cs;
		size = coxeter_system->size;
		word = new int[size];
		one_line = new int[size+1];
		for (int i = 0; i < size; i++)
		{ word[i] = w[i]; one_line[i] = ol[i]; }
		one_line[size] = ol[size];
		length = get_length();
	}

	// from reduced expression
	CoxeterElement::CoxeterElement(CoxeterSystem* cs, int red[], int len)
	{
		coxeter_system = cs;
		size = coxeter_system->size;
		word = new int[size];
		one_line = new int[size+1];
		for (int i = 0; i < size; i++)
		{
			word[i] = 1;
			one_line[i] = i+1;
		}
		one_line[size] = size+1;
		length = 0;

		for (int i = 0; i < len; i++)
		{ right_multiply(red[i]); }
	}

	// copy constructor
	CoxeterElement::CoxeterElement(const CoxeterElement& cp)
	{
		coxeter_system = cp.coxeter_system;
		size = coxeter_system->size;

		word = new int[size];
		one_line = new int[size+1];
		for (int i = 0; i < size; i++)
		{
			word[i] = cp.word[i];
			one_line[i] = cp.one_line[i];
		}
		one_line[size] = cp.one_line[size];

		length = cp.length;
	}

	//CoxeterElement::operator=(const CoxeterElement& cp);   // = operator copy-constructor: for passing by reference in functions and creating copies on the heap.

	int CoxeterElement::equals(CoxeterElement* cp)
	{
		if (cp->size != this->size) { return 0; }
		int ret = 1;
		for (int i = 0; i < this->size; i++)
		{
			if (cp->word[i] != this->word[i]) { ret = 0; }
		}
		return ret;
	}

	void CoxeterElement::print()
	{
		cout << "[ ";
		for (int i = 0 ; i < size; i++ ) { cout << word[i] << " "; }
		cout << "]";

		cout << " { ";
		for (int i = 0 ; i < size+1; i++ ) { cout << one_line[i] << " "; }
		cout << "}";
	}

	void CoxeterElement::print_reduced_expression()
	{
		int* r = new int[length];
		get_reduced_expression(r);
		cout << "( ";
		for (int i = 0; i < length; i++)
		{ cout << r[i] << " "; }
		cout << ")";
		delete[] r;
	}

	void CoxeterElement::right_multiply(int s)
	{
		for (int i = 0; i < size; i++)
		{
			if (coxeter_system->get_exponent(s,i) == 0) { word[i] = word[i]; }  // do nothing.
			else if (coxeter_system->get_exponent(s,i) == 2) { word[i] = word[i]; }  // commutes so do nothing.
			else if (coxeter_system->get_exponent(s,i) == 3) { word[i] = word[s] + word[i]; }
			else if (coxeter_system->get_exponent(s,i) == 4) { word[i] = (2*word[s]) + word[i]; }
			else if (coxeter_system->get_exponent(s,i) == 6) { word[i] = (3*word[s]) + word[i]; }
			else if (coxeter_system->get_exponent(s,i) == -1) { word[i] = (2*word[s]) + word[i]; }  // meaning infinity
			else { word[i] = 0; cout << "ERROR:  coxeter matrix entry " << coxeter_system->get_exponent(s,i) << " not supported." << endl; } 
		}

		if ( word[s] > 0 ) { length++; } else { length--; }

		word[s] = 0 - word[s];

		// perform type-dependent multiplication on one_line.
		coxeter_system->right_multiply(one_line, s);
	}

	// WARNING:  left multiplication is not optimized.
	void CoxeterElement::left_multiply(int s)
	{
		int* reduced = new int[ length ];
		get_reduced_expression(reduced);

		CoxeterElement t = CoxeterElement(coxeter_system);
		int count_moves = 1;
		t.right_multiply(s);
		for (int i = 0; i < length; i++)
		{
			if (t.word[ reduced[i] ] > 0) { count_moves++; } else { count_moves--; }
			t.right_multiply(reduced[i]);
		}

		for (int i = 0; i < size; i++)
		{
			word[i] = t.word[i];
			one_line[i] = t.one_line[i];
		}
		one_line[size] = t.one_line[size];
		length = count_moves;

		delete [] reduced;
	}

	int CoxeterElement::get_length()
	{
		int count_moves = 0;
		get_reduced_expression(NULL, count_moves);
		return count_moves;
	}

	void CoxeterElement::get_reduced_expression(int reduced[])
	{
		int c = 0;
		get_reduced_expression(reduced, c);
	}

	// requires an allocated array of length at least this->length.
	void CoxeterElement::get_reduced_expression(int reduced[], int& count_moves)
	{
		count_moves = 0;

		CoxeterElement t = *this;  // using copy-constructor, t is on the heap...
		
		while (1==1)
		{
			int move = 0;
			for (move = 0; move < size; move++)
			{
				if (t.word[move] < 0) { break; }
			}

			if (move >= size) { break; }

			if (reduced != NULL) { reduced[length-1-count_moves] = move; }
			count_moves++;
			t.right_multiply(move);
		}
	}

	int CoxeterElement::get_rank(int reduced[])
	{
		vector<int> support;
		
		for (int i = 0; i < length; i++)
		{
			int contains = 0;
			for (int j = 0; j < support.size(); j++)
			{ 
				if (support[j] == reduced[i]) { contains = 1; } 
			}
			if (contains == 0) { support.push_back( reduced[i] ); }
		}

		return support.size();
	}

	void CoxeterElement::mask(CoxeterElement* r, long long m1, long long m2, int reduced[])
	{
		mask(r, m1, m2, length, reduced);
	}

	void CoxeterElement::mask(CoxeterElement* q, long long m1, long long m2, int imax, int reduced[])
	{
		for (int i = 0; i <= imax; i++)
		{
		  if ( i < DATATYPE_BITLENGTH-1 )
		  {
			if (((m1 >> i) & 1) == 1)
			{
				if (DEBUG_VERBOSE) { cout << "mask(); r mult by " << reduced[i] << " (" << i << ")" << endl; }
				q->right_multiply( reduced[i] );
			}
		  }
		  else
		  {
			if (( (m2 >> (i - (DATATYPE_BITLENGTH-1)) ) & 1) == 1)
			{
				if (DEBUG_VERBOSE) { cout << "mask(); r mult by " << reduced[i] << " (" << i << ")" << endl; }
				q->right_multiply( reduced[i] );
			}
		  }
		}
	}

	int CoxeterElement::defect(long long m1, long long m2, int reduced[])
	{
		int def = 0;

		for (int i = 0; i < length-1; i++)  // NOTE:  only up to length-1...
		{
			CoxeterElement* q = new CoxeterElement(coxeter_system);
			mask(q, m1, m2, i, reduced);
			int l1 = q->get_length();
			q->right_multiply( reduced[i+1] );
			int l2 = q->get_length();
			if (l2 < l1) { def++; }
			delete q;
		}

		return def;
	}

	int CoxeterElement::deodhar()
	{
		int ret = 1;
		int* reduced = new int[ length ];
		
		get_reduced_expression(reduced);

		// NOTE:  for each mask, we only need STRICT subwords.  The Deodhar statistic is not satisfied for w=x.
		
		long long target_m1;
		long long target_m2;

		if ( length < DATATYPE_BITLENGTH-1 )
		{ 
			target_m1 = (long long) pow(2.0, (double) length)-1; 
			target_m2 = (long long) 0;
		}
		else
		{
			target_m1 = (long long) pow(2.0, (double) (DATATYPE_BITLENGTH-1))-1;
			target_m2 = (long long) pow(2.0, (double) (length-DATATYPE_BITLENGTH+1))-1;
		}

		if (DEBUG_VERBOSE) { cout << "  deodhar():  length " << length << " target m1 " << target_m1 << " target m2 " << target_m2 << endl; }

		for (long long m2 = 0; m2 <= target_m2; m2++)
		{
			for (long long m1 = 0; (m2 < target_m2 && m1 < pow(2.0, (double) (DATATYPE_BITLENGTH-1))) || (m2 == target_m2 && m1 < target_m1); m1++)
			{

		  		if (DEBUG_VERBOSE) { 
		    		cout << "(mask: ";
		    		for (int ii = 0; ii < DATATYPE_BITLENGTH-1; ii++) { cout << ((m1 >> ii) & 1); }
		    		cout << " ";
		    		for (int ii = 0; ii < DATATYPE_BITLENGTH-1; ii++) { cout << ((m2 >> ii) & 1); }
		    		cout << ")" << endl;
				}

		  		CoxeterElement* q = new CoxeterElement(coxeter_system);  // was size
		  		mask(q, m1, m2, reduced);
		  		int ds = (length - q->get_length() - 1) - (2*defect(m1, m2, reduced));
		  		// want ds >= 0 in order to satisfy Thm 3.
		  		if (DEBUG_VERBOSE) { cout << "d stat: " << ds << " = " << length << "-" << q->get_length() << "-1-2*" << defect(m1, m2, reduced) << endl; }

		  		if (ds < 0) 
				{ 
		    			if (DEBUG_VERBOSE) { 
						cout << "(deodhar test failed on: ";
		  				for (int ii = 0; ii < DATATYPE_BITLENGTH-1; ii++) { cout << ((m1 >> ii) & 1); }
		  				cout << " ";
			  			for (int ii = 0; ii < DATATYPE_BITLENGTH-1; ii++) { cout << ((m2 >> ii) & 1); }
		  				cout << ")" << endl;
		    			}
		  			ret = 0; 
		  			delete q; 
		  			break;   // comment out the break to see _all_ masks.
		  		}
		  		else
		  		{
		    			delete q;
		  		}
			}
		}

		delete[] reduced;
		return ret;
	}

	CoxeterElement::~CoxeterElement()
	{
		delete[] word;
		delete[] one_line;
	}

	void CoxeterElement::print_heap()
	{
		int* reduced = new int[ length ];
		get_reduced_expression(reduced);

		int level[length];
		for (int i = 0; i < length; i++) { level[i] = 0; }

		int heap[size];
		for (int i = 0; i < size; i++) { heap[i] = 0; }

		cout << "length : " << length << endl;
		for (int i = 0; i < length; i++) { cout << reduced[i] << " "; }
		cout << endl;

		for (int i = 0; i < length; i++)
		{
			for ( int m = 0 ; m < size; m++ )
			{
				if (this->coxeter_system->get_exponent(m, reduced[i]) >= 3) 
				{ heap[m] = heap[reduced[i]]+1; }
			}

			level[ i ] = heap[reduced[i]];
		}

		cout << endl;

		for (int le = length-1; le >= 0; le--)
		{
			for (int seek = 0 ; seek < size; seek++)
			{
				int found = 0;
				for (int i = 0; i < length; i++)
				{
				if (level[i] == le && reduced[i] == seek) { cout <<  "* ";  found = 1; }
				}

				if (found == 0) { cout << "  "; }
	 		}
			cout << endl;
		}

		// add coalescing code?  For each entry, l->r if there is nothing immediately above (but there is eventually), then move the elt up.

		cout << endl;
	}


	// This checks type A, B, D-style one_line pattern containment:  i.e. the bars must be in the same position, and the digits flatten.
	// WARNING:  Currently only used/extensively tested with type D.  There is a size difference in the one_line array for type A.
	// WARNING:  To optimize, we should really only calculate the flattening of w once, and check the entire _list_ of patterns inside the loop.  However, the patterns are not all the same size.
	int CoxeterElement::contains_one_line_pattern(CoxeterElement& pattern)
	{
		// look at all subwords of size pattern->size in this.
		// see if they flatten to the first pattern->size entries of pattern, with the same bar pattern.
		int any_match = 0;
		
		int k = pattern.size;

		long lower_swm = 0;
		for (int i = 0; i < k; i++ ) { lower_swm = lower_swm + (long) pow(2.0, (double) i); }
		long upper_swm = 0;
		for (int i = 1; i <= k; i++ ) { upper_swm = upper_swm + (long) pow(2.0, (double) ((this->size)-i)); }

		for (long swm = lower_swm; swm <= upper_swm; swm++)
		{
			// make sure there are exactly k positions.
			int on = 0;
			for (int ii = 0; ii < 32-1; ii++) { if ( ((swm >> ii) & 1) == 1 ) { on++; } }
			if ( on != k ) { continue; }

			// map to indicies to get subword.
			int* sub_one_line = new int[ k ];
			int iii=0;
			for (int ii = 0; ii < 32-1; ii++) { if ( ((swm >> ii) & 1) == 1 ) { sub_one_line[iii] = ii; iii++; } }

			// flatten and test.
			int* flattened = new int[k];
			for (int i = 0; i < k; i++) { flattened[i] = 0; }
			int next_digit = 1;
			for (int search_for = 1; search_for <= this->size; search_for++)
			{
			  for (int i = 0; i < k; i++)
			  {
				  if (abs(this->one_line[ sub_one_line[i] ]) == search_for)
				  {
					  flattened[i] = next_digit;
					  if (this->one_line[ sub_one_line[i] ] < 0) { flattened[i] = 0 - flattened[i]; }
					  next_digit++;
				  }
			  }
			}

			int matches = 1;
			for (int i = 0; i < k; i++)
			{
				if (flattened[i] != pattern.one_line[i]) { matches = 0; }
			}

			// test output.
			if (matches == 1) { any_match = 1;  if (DEBUG_VERBOSE) { cout << " (one_line pattern match on "; pattern.print(); cout << ") " << endl; }  break; }
			 
			if (DEBUG_VERBOSE)
			{	
			  cout << "n = " << this->size << " k = " << k << " : " << lower_swm << " to " << upper_swm << ": ";
			  for (int ii = 0; ii < 32-1; ii++) { cout << ((swm >> ii) & 1); }
			  cout << " : ";
			  for (int i = 0; i < k; i++) { cout << sub_one_line[i]; }
			  cout << " : ";
			  for (int i = 0; i < k; i++) { cout << flattened[i]; }
			  cout << "  Match on "; pattern.print(); cout << " : " << matches << endl;
			}
		}

		return any_match;
	}


//////////////////////////////////////////////////////////////////////
//
//  End of eriksson library code.
//
//////////////////////////////////////////////////////////////////////


