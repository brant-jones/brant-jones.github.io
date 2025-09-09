#include "CoxeterElement.h"
#include "Masks.h"

/////////////////////////////////////////////////////////////////////
//
//  This is code which classifies the minimally non-Deodhar elements
//  of various Coxeter groups.
//  
//  We also verify the "mu in {0,1}" property for 
//  Deodhar elements of these Coxeter groups.
//
/////////////////////////////////////////////////////////////////////

static const int DEBUG_VERBOSE_GEN = 0;
static const int VERBOSE = 0;

// generates all the elements above t in the 2-weak order, puts them in PROCESSED_ELEMENTS.
int generate_up_ideal(CoxeterElement& t, int max_length, vector<CoxeterElement>& PROCESSED_ELEMENTS)
{
	list<CoxeterElement> toproc;
	toproc.push_back(t);

	while (!toproc.empty())
	{
		CoxeterElement current = toproc.front();
		toproc.pop_front();

		if (max_length > 0 && current.length > max_length) { return 0; }

		int* reduced = new int[current.length];
		current.get_reduced_expression(reduced);

		// see if current elt has already been processed
		int match_exists = 0;
		for (int i = 0; i < PROCESSED_ELEMENTS.size(); i++)
		{
			if ( current.equals(&(PROCESSED_ELEMENTS[i])) == 1 )
			{ match_exists = 1; }
		}

		if ( match_exists != 0 )
		{
			delete reduced;
			continue;
		}

		if (DEBUG_VERBOSE_GEN)
		{ cout << "adding "; current.print(); current.print_reduced_expression(); }

		// add current elt
		PROCESSED_ELEMENTS.push_back( current );

		// hit current on the right with all possible generators, s.t. it is short-braid-avoiding
		for (int i = 0; i < current.size; i++)
		{
		  if (current.word[i] > 0) 
		  { 
			int non_comms = 0;
			
			// search backwards for generator i:
			int k = 0;
			for ( k = current.length-1 ; k >= 0; k-- ) { if (reduced[k] == i) { break; } }

			if ( k >= 0 )
			{
			  // search forwards for non commuting gens...
			  for ( int m = k+1 ; m < current.length; m++ )
			  {
				if (current.coxeter_system->get_exponent(i, reduced[m]) >= 3) { non_comms++; }
			  }
			}
			else { non_comms = 2; }
			
			if ( non_comms >=2 )
			{
		    	  if (DEBUG_VERBOSE_GEN) { cout << "multiplying by " << i << " and adding to list." << endl; }
			  CoxeterElement v = CoxeterElement(current.coxeter_system, current.word, current.one_line);
			  v.right_multiply(i);
			  toproc.push_back(v);
			}
			else
			{
		    	  if (DEBUG_VERBOSE_GEN) { cout << "multiplying by " << i << " would create short-braid." << endl; }
			}
		  }
		}

		// hit current on the left with all possible generators, s.t. it is short-braid-avoiding
		for (int i = 0; i < current.size; i++)
		{
			int non_comms = 0;
			
			// search forwards for generator i:
			int k = 0;
			for ( k = 0; k < current.length; k++ ) { if (reduced[k] == i) { break; } }

			if ( k < current.length )
			{
			  // search back for non commuting gens...
			  for ( int m = k ; m >= 0; m-- )
			  {
				if (current.coxeter_system->get_exponent(i, reduced[m]) >= 3) { non_comms++; }
			  }
			}
			else { non_comms = 2; }
			
			if ( non_comms >=2 )
			{
		    	  if (DEBUG_VERBOSE_GEN) { cout << "left multiplying by " << i << " and adding to list." << endl; }
			  CoxeterElement v = CoxeterElement(current.coxeter_system, current.word, current.one_line);
			  v.left_multiply(i);
			  if (v.get_length() > current.length) { toproc.push_back(v); }
			}
			else
			{
		    	  if (DEBUG_VERBOSE_GEN) { cout << "left multiplying by " << i << " would create short-braid." << endl; }
			}
		}

		delete reduced;
	} // end while there are elements yet to process
}

int generate_all_elements_breadth_first(CoxeterSystem* coxeter_system, int max_length)
{
	int deodhar_count = 0;  // Keep enumeration for futher directions section of paper.
	int total_count = 0;    // This won't match the total # elts in group since we cut the recursion when we find a bad pattern.

	vector<CoxeterElement> PROCESSED_ELEMENTS;
	vector<CoxeterElement> NON_DEODHAR_PATTERNS;
	list<CoxeterElement> toproc;

	// initialization of bad D8 1-line pattern:
	int D8_PATTERN_OL[9] = { -1, 6, 7, 8, -5, 2, 3, 4, 9}; int D8_PATTERN_W[8] = {5, 5, 1, 1, -11, 5, 1, 1};
	CoxeterElement D8_PATTERN = CoxeterElement(&D8, D8_PATTERN_W, D8_PATTERN_OL); 

	CoxeterElement t = CoxeterElement(coxeter_system);  // create identity elt.
	toproc.push_back(t);

	int current_length = -1;

	while (!toproc.empty())
	{
		CoxeterElement current = toproc.front();

		toproc.pop_front();

		if (max_length > 0 && current.length > max_length) { return 0; }
		if (current_length < current.length) 
		{ 
			current_length = current.length;  
			cout << "  (evaluating length " << current_length << " elements, with " << toproc.size() << " elements left to process...) " << endl; 
		}

		int* reduced = new int[current.length];
		current.get_reduced_expression(reduced);

		// see if current elt has already been processed
		int match_exists = 0;
		for (int i = 0; i < PROCESSED_ELEMENTS.size(); i++)
		{
			if ( current.equals(&(PROCESSED_ELEMENTS[i])) == 1 )
			{ match_exists = 1; }
		}

		if ( match_exists == 1 )
		{
			delete reduced;
			continue;
		}

		if (DEBUG_VERBOSE_GEN)
		{ cout << "adding "; current.print(); current.print_reduced_expression(); }

		PROCESSED_ELEMENTS.push_back( current );
		total_count++;

	map<string, int> mus;
	int DEFAULT_MAP_VALUE = mus["undef"];
	int dt = 1;
        Masks masks(current.length);
        while (masks.exhausted() == 0)
        {
		if (!masks.proper()) { masks.next(); continue; }

                // calculate the contribution of this mask to the Kazhdan-Lusztig basis element.
                int defect_count = 0;
                CoxeterElement tc(coxeter_system);
                for (int i = 0; i < current.length; i++)
                {
                        // build indexing element using 1-entries through this position i.
                        if (masks.get_value(i) == 1)
                        { tc.right_multiply(reduced[i]); }

                        // see if this position is a defect:  is reduced_expression[i+1] a right descent for t?
                        if (i < current.length-1 && tc.has_right_descent(reduced[i+1]))
                        { defect_count++; }
                }

		// Check the Deodhar statistic:  1 = mu mask, 0 = not Deodhar.
		int deodhar_statistic = (current.length - tc.get_length()) - (2*defect_count);
		if (deodhar_statistic <= 0 && masks.proper())
		{
			dt = 0;
			break;
		}

		// Check for mu in {0, 1}.
		if (deodhar_statistic == 1)
		{
				  // The current mask is a mu-mask, so update mu values.
				  string tc_reduced;
				  tc.sprint_reduced_expression(tc_reduced);
				  
				  if (mus[ tc_reduced ] != DEFAULT_MAP_VALUE)
				  {
				     if(VERBOSE)
				     {
						cout << "WARNING:  0-1 mu property failed on: ";
						string mm;
						masks.sprint(mm);
						cout << mm << endl;
					//	print_reduced_expression();
						cout << " (but this may be a non-Deodhar elt.)." << endl;
				     }	
				  }
				  else
				  { mus[ tc_reduced ] = 0; }

				  mus[ tc_reduced ]++;
		}

                masks.next();
        }

	  	if ( dt == 1 ) 
		{ 
			deodhar_count++; 

			// Print mu values.
			if (VERBOSE) 
			{
    			cout << "w = ";
			current.print_reduced_expression();
			cout << endl;
			}
			map<string, int>::iterator iter;   
  			for( iter = mus.begin(); iter != mus.end(); iter++ ) 
			{
			  if (VERBOSE)
			  {
    				cout << "mu = " << iter->second << " for x = " << iter->first << endl;
			  }
			  if (iter->second > 1) 
			  { 
			    cout << "ERROR:  Found NON-01 MU VALUE:  ";  current.print();  current.print_reduced_expression();
    			    cout << "    mu = " << iter->second << " for x = " << iter->first << endl;
			  }
  			}

			// check special 1-line pattern for any system _containing_ D8.  
			if ( (current.coxeter_system->equals(&D8) == 1) || (current.coxeter_system->equals(&D9) == 1) || (current.coxeter_system->equals(&D10) == 1) )
			{
				if ( current.contains_one_line_pattern(D8_PATTERN) )
				{ cout << "ERROR:  cannot use D8 1-line pattern for Deodhar characterization:  "; current.print(); current.print_reduced_expression(); cout << endl; }
			}
		}

	  	if (dt == 0)
	  	{
			if (VERBOSE) 
			{ cout << "  Found non-Deodhar element:  ";  current.print();  current.print_reduced_expression(); }

			int contains_pattern = 0;

			// see if any patterns in NON_DEODHAR_PATTERN list are _equal_ to current.
			for (int i = 0; i < NON_DEODHAR_PATTERNS.size(); i++)
			{
				if (current.equals(&(NON_DEODHAR_PATTERNS[i])))
				{
					contains_pattern = 1;
					break;
				}
			}

			// If not, add current, and add up-ideals in 2-weak order generated by current and it's Coxeter embeddings to NON_DEODHAR_PATTERN list.
			if (contains_pattern == 0)
			{
				// check special 1-line pattern for any system _containing_ D8.  
				// Would be nice if there were a more generic way to do this.
				if ( (current.coxeter_system->equals(&D8) == 1) || (current.coxeter_system->equals(&D9) == 1) || (current.coxeter_system->equals(&D10) == 1) )
				{
					if ( current.contains_one_line_pattern(D8_PATTERN) )
					{ contains_pattern = 1; }
				}

				if (contains_pattern == 1)
				{
				  cout << "  (D8 1-line pattern found in "; current.print(); current.print_reduced_expression(); cout << " of rank " << current.get_rank(reduced) << ")" << endl;
				}
				else
				{
				  cout << "Found MINIMAL PATTERN of rank " << current.get_rank(reduced) << ":  ";  current.print();  cout << " "; current.print_reduced_expression();  cout << endl;
				}

				int rv = generate_up_ideal(current, 0, NON_DEODHAR_PATTERNS);
		
				for (int i = 0; i < current.coxeter_system->automorphism_group_size; i++)
				{
			  		CoxeterElement v = CoxeterElement(current.coxeter_system);
					if (VERBOSE) { cout << " adding related graph automorphic element:  "; }
					for (int j = 0; j < current.length; j++)
					{
			  			v.right_multiply(current.coxeter_system->get_automorphism(i, reduced[j]));
					}
					if (VERBOSE) { v.print(); cout << endl; }

					rv = generate_up_ideal(v, 0, NON_DEODHAR_PATTERNS);
				}
			}

			if (contains_pattern == 1 && VERBOSE)
			{
				cout << "  (contains previous pattern, so breaking)..." << endl;
			}

			// In any event, by lemma:  once we find a non-Deodhar pattern, there's no need to recurse further in 2-sided weak order...
			delete reduced;
			continue;
		}

		// Now, if the current element is Deodhar, extend it in all possible (short-braid-avoiding) ways, and add these extensions to the list of elements to process.
		for (int i = 0; i < current.size; i++)
		{
		  if (current.word[i] > 0) 
		  { 
			int non_comms = 0;
			
			// search backwards for generator i:
			int k = 0;
			for ( k = current.length-1 ; k >= 0; k-- ) { if (reduced[k] == i) { break; } }

			if ( k >= 0 )
			{
			  // search forwards for non commuting gens...
			  for ( int m = k+1 ; m < current.length; m++ )
			  {
				if (current.coxeter_system->get_exponent(i, reduced[m]) >= 3) { non_comms++; }
			  }
			}
			else { non_comms = 2; }
			
			if ( non_comms >=2 )
			{

		    	  if (DEBUG_VERBOSE_GEN) { cout << "multiplying by " << i << " and adding to list." << endl; }
			  CoxeterElement v = CoxeterElement(current.coxeter_system, current.word, current.one_line);
			  v.right_multiply(i);
			  toproc.push_back(v);
			}
			else
			{
		    	  if (DEBUG_VERBOSE_GEN) { cout << "multiplying by " << i << " would not be fc." << endl; }
			}
		  }
		}

		delete reduced;
	} // end while there are elements yet to process

	cout << "Finished:  found " << deodhar_count << " Deodhar elements (out of " << total_count << " short-braid-avoiding elements processed, used " << NON_DEODHAR_PATTERNS.size() << " non-Deodhar patterns).  Also verified mu in {0,1} property for these Deodhar elements (check stdout for errors)." << endl;

        if (VERBOSE)
        { for (int i = 0; i < NON_DEODHAR_PATTERNS.size(); i++) { NON_DEODHAR_PATTERNS[i].print(); NON_DEODHAR_PATTERNS[i].print_reduced_expression(); cout << endl; } }

	cout << "  (Consistency:  check that #short-braid-avoiding elts = " << (NON_DEODHAR_PATTERNS.size() + deodhar_count) <<  " = total non-Deodhar elts + Deodhar elts.)" << endl;
	cout << endl;
}

int main(int argc, char* argv[])
{
	/////////////////////////////////////////////////////////////
	//
	// Here are # fully-commutative elts in ranks 2 - 14:  (See Stembridge for details.)
	// Type A:  {5,14,42,132,429,1430,4862,16796,58786,208012,742900,2674440,9694845}
	// Type B:  {7,24,83,293,1055,3860,14299,53481,201551,764217,2912167,11143499,42791039}
	// Type D:  {4,14,48,167,593,2144,7864,29171,109173,411501,1560089,5943199,22732739}
	// Type E:
	//   E6:  662
	//   E7:  2670
	//   E8:  10846
	//
	/////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////
	//  Finite exceptional types.
	/////////////////////////////////////////////////////////////

	cout << "Type G2: " << endl;
	G2.print_matrix(); cout << endl;
	generate_all_elements_breadth_first(&G2, 0);
	cout << "G2 has 5 short-braid-avoiding elements." << endl;
	cout << endl << endl;

	cout << "Type F4: " << endl;
	F4.print_matrix(); cout << endl;
	generate_all_elements_breadth_first(&F4, 0);
	cout << "F4 has 42 short-braid-avoiding elements." << endl;
	cout << endl << endl;

	/////////////////////////////////////////////////////////////
	//  Minimally non-Deodhar families.
	/////////////////////////////////////////////////////////////

	cout << endl << "Minimal non-Deodhar embedded factor patterns: " << endl;;

	cout << "Generating minimal patterns for A7 (includes linear types BC, F4, G2, H3, H4): " << endl;
	A7.print_matrix(); cout << endl;
	generate_all_elements_breadth_first(&A7, 0);
	cout << "A7 has 1430 short-braid-avoiding elements." << endl;
	cout << endl << endl;

	cout << "Generating minimal patterns for D8 (excluding 1-line pattern {-1, 6, 7, 8, -5, 2, 3, 4, 9}): " << endl;
	D8.print_matrix(); cout << endl;
	generate_all_elements_breadth_first(&D8, 0);
	cout << "D8 has 7864 short-braid-avoiding elements." << endl;
	cout << endl << endl;

	cout << "Generating minimal patterns for E7: " << endl;
	E7.print_matrix(); cout << endl;
	generate_all_elements_breadth_first(&E7, 0);
	cout << "E7 has 2670 short-braid-avoiding elements." << endl;
	cout << endl << endl;

	/////////////////////////////////////////////////////////////
	//
	//  Some (optional) checks for other types
	//
	/////////////////////////////////////////////////////////////

/*
	cout << endl << "Enumerative data and consistency checks: " << endl;;

        cout << "Type E6: " << endl;
        E6.print_matrix(); cout << endl;
        generate_all_elements_breadth_first(&E6, 0);
        cout << "E6 has 24 short-braid-avoiding elements." << endl;
        cout << endl << endl;

	cout << "Type E8: " << endl;
	E8.print_matrix(); cout << endl;
	generate_all_elements_breadth_first(&E8, 0);
	cout << "E8 has 10846 short-braid-avoiding elements." << endl;
	cout << endl << endl;

        cout << "Type D3: " << endl;
        D3.print_matrix(); cout << endl;
        generate_all_elements_breadth_first(&D3, 0);
        cout << "D3 has 24 short-braid-avoiding elements." << endl;
        cout << endl << endl;

        cout << "Type D4: " << endl;
        D4.print_matrix(); cout << endl;
        generate_all_elements_breadth_first(&D4, 0);
        cout << "D4 has 83 short-braid-avoiding elements." << endl;
        cout << endl << endl;

        cout << "Type D5: " << endl;
        D5.print_matrix(); cout << endl;
        generate_all_elements_breadth_first(&D5, 0);
        cout << "D5 has 293 short-braid-avoiding elements." << endl;
        cout << endl << endl;

        cout << "Type D6: " << endl;
        D6.print_matrix(); cout << endl;
        generate_all_elements_breadth_first(&D6, 0);
        cout << "B6 has 1055 short-braid-avoiding elements." << endl;
        cout << endl << endl;

        cout << "Type D7: " << endl;
        D7.print_matrix(); cout << endl;
        generate_all_elements_breadth_first(&D7, 0);
	cout << "D7 has 2144 short-braid-avoiding elements." << endl;
        cout << endl << endl;

	cout << "Type B3: " << endl;
	B3.print_matrix(); cout << endl;
	generate_all_elements_breadth_first(&B3, 0);
	cout << "B3 has 24 short-braid-avoiding elements." << endl;
	cout << endl << endl;

	cout << "Type B4: " << endl;
	B4.print_matrix(); cout << endl;
	generate_all_elements_breadth_first(&B4, 0);
	cout << "B4 has 83 short-braid-avoiding elements." << endl;
	cout << endl << endl;

	cout << "Type B5: " << endl;
	B5.print_matrix(); cout << endl;
	generate_all_elements_breadth_first(&B5, 0);
	cout << "B5 has 293 short-braid-avoiding elements." << endl;
	cout << endl << endl;

	cout << "Type B6: " << endl;
	B6.print_matrix(); cout << endl;
	generate_all_elements_breadth_first(&B6, 0);
	cout << "B6 has 1055 short-braid-avoiding elements." << endl;
	cout << endl << endl;

	cout << "Type B7: " << endl;
	B7.print_matrix(); cout << endl;
	generate_all_elements_breadth_first(&B7, 0);
	cout << "B7 has 3860 short-braid-avoiding elements." << endl;
	cout << endl << endl;

	cout << "Type A8: " << endl;
	A8.print_matrix(); cout << endl;
	generate_all_elements_breadth_first(&A8, 0);
	cout << "A8 has 4862 short-braid-avoiding elements." << endl;
	cout << endl << endl;

	cout << "Type A9: " << endl;
	A9.print_matrix(); cout << endl;
	generate_all_elements_breadth_first(&A9, 0);
	cout << "A9 has 16796 short-braid-avoiding elements." << endl;
	cout << endl << endl;

	cout << "Type A10: " << endl;
	A10.print_matrix(); cout << endl;
	generate_all_elements_breadth_first(&A10, 0);
	cout << "A10 has 58786 short-braid-avoiding elements." << endl;
	cout << endl << endl;

	cout << "Type D9 (excluding 1-line pattern {-1, 6, 7, 8, -5, 2, 3, 4, 9}): " << endl;
	D9.print_matrix(); cout << endl;
	generate_all_elements_breadth_first(&D9, 0);
	cout << "D9 has 29171 short-braid-avoiding elements." << endl;
	cout << endl << endl;

	cout << "Type D10 (excluding 1-line pattern {-1, 6, 7, 8, -5, 2, 3, 4, 9}): " << endl;
	D10.print_matrix(); cout << endl;
	generate_all_elements_breadth_first(&D10, 0);
	cout << "D10 has 109173 short-braid-avoiding elements." << endl;
	cout << endl << endl;
*/
}

