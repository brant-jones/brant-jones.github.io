#include "CoxeterElement.h"
#include "Masks.h"

/////////////////////////////////////////////////////////////////////
//
//  This is code which tests a given element for the Deodhar 
//  property.
//
/////////////////////////////////////////////////////////////////////

static const int DEBUG_VERBOSE_GEN = 0;
static const int VERBOSE = 0;

struct strCmp {
    bool operator()( string s1, string s2 ) const {
      return ( s1 < s2 );
    }
};


/////////////////////////////////////////////////////////////////////
// Input:  Coxeter matrix (type), w redexp, x redexp, -mu flag (optional)
// Output:  Mu coefficient, P_{x,w} if x provided, or the entire C'_{x,w} if not.
/////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
	if (argc < 2) { 
		cout << "No arguments given." << endl;  
                cout << "Usage:  ./deodhar A3 -w 1021" << endl;
                cout << "Given a reduced expression for a Deodhar element," << endl;
                cout << "returns the Kazhdan-Lusztig basis element." << endl;
		cout << "Optional argument:  -x <reduced expression> only prints the polynomial P_{x,w}(q)." << endl;
		cout << "Optional argument:  -mu prints the mu-coefficients." << endl;
		cout << "Optional argument:  -masks prints all masks." << endl;
		return 0; 
	}

	string na = argv[1];
	cout << "Coxeter type " << na << " with Coxeter matrix: " << endl;

	CoxeterNames names;
	CoxeterSystem* coxeter_system = names.names[ na ];
	if (coxeter_system == NULL) { cout << "Coxeter type not supported." << endl;  return 0; }
	coxeter_system->print_matrix();

	string w_s = "";
	string x_s = "";
	int look_for_mu = 0;
	int print_all_masks = 0;
	int w_is_not_deodhar = 0;

	for (int i = 2; i < argc; i++)
	{
		string a = argv[i];
		if (a == "-w") { i++; w_s = argv[i]; }
		else if (a == "-x") { i++; x_s = argv[i]; }
		else if (a == "-mu") { look_for_mu = 1;  print_all_masks = 1; }
		else if (a == "-masks") { print_all_masks = 1; }
		//else if (a == "-conj") { test_conj = 1; }
	}

	if (w_s == "") { cout << "No w expression given." << endl;  return 0; }

        int w_red_l = w_s.length();
        int w_red[w_red_l];

        for (int i=0; i < w_red_l; i++)
        {
		stringstream strstr(w_s.substr(i,1));
                strstr >> w_red[i];
                //w_red[i] = atoi( &(w_s.c_str()[i]) );
        }

	CoxeterElement w = CoxeterElement(coxeter_system, w_red, w_red_l);
	string w_name;
	w.sprint_reduced_expression(w_name);
	cout << "w is using reduced expression " << endl;
	cout << "....:  " << w_name << "." << endl;

        int x_red_l = x_s.length();
        int x_red[w_red_l];
	string x_name;

        for (int i=0; i < x_red_l; i++)
        {
	  stringstream strstr(x_s.substr(i,1));
          strstr >> x_red[i];
        }

	CoxeterElement x = CoxeterElement(coxeter_system, x_red, x_red_l);

	if (x_s != "")
	{
	  x.sprint_reduced_expression(x_name);
	  cout << "x is using reduced expression " << x_name << "." << endl;
	}

	map< string, map<string, int> > kl_basis_element;  // T_oneline is the first index, and q^defects is the second index.

	int length = w.length;
	int reduced_expression[length];

	w.get_reduced_expression(reduced_expression);

        Masks masks(length);
        while (masks.exhausted() == 0)
        {
		// print mask to string.
		string mask_as_string;
		masks.sprint(mask_as_string);

                // calculate the contribution of this mask to the Kazhdan-Lusztig basis element.
                int defect_count = 0;
                CoxeterElement t(coxeter_system);
                for (int i = 0; i < length; i++)
                {
                        // build indexing element using 1-entries through this position i.
                        if (masks.get_value(i) == 1)
                        { t.right_multiply(reduced_expression[i]); }

                        // see if this position is a defect:  is reduced_expression[i+1] a right descent for t?
                        if (i < length-1 && t.has_right_descent(reduced_expression[i+1]))
                        { 
				defect_count++; 
				mask_as_string.erase(2*(i+2)+1, 1);
				mask_as_string.insert(2*(i+2)+1, "d");
			}
                }

                string t_basis;
		t.sprint_reduced_expression(t_basis);

		// Check the Deodhar statistic:  1 = mu mask, 0 = not Deodhar.
		int deodhar_statistic = (length - t.get_length()) - (2*defect_count);
		if (deodhar_statistic <= 0 && masks.proper())
		{
			if (print_all_masks == 1)
			{
			  w_is_not_deodhar = 1;
			}
			else
			{
			  cout << "The element w has non-Deodhar mask:" << endl;
			  cout << mask_as_string << endl;
			  return 0;
			}
		}

                string q_monomial;
                if (defect_count == 0) { q_monomial = "1"; }
                else if (defect_count == 1) { q_monomial = "q"; }
                else
                {
                        stringstream s;
                        s << defect_count;
                        q_monomial = "q^" + s.str();
                }

                kl_basis_element[t_basis][q_monomial]++;

		if (print_all_masks == 1)
		{
			if ((x_s == "") || (t_basis == x_name)) 
			{
			  if (look_for_mu == 0 || (look_for_mu == 1 && deodhar_statistic == 1))
			  {
			  cout << "mask:  ";
			  cout << mask_as_string << ":  " << q_monomial;

			  if (deodhar_statistic == 1)
 			  {
			    cout << " (mu mask) ";
			  }
			  if (deodhar_statistic <= 0 && masks.proper())
			  {
			    cout << " (not Deodhar) ";
			  }
			  cout << endl;
			  }
			}
		}

                //cout << "  adding " << q_monomial << " for " << t_basis << "  (" << defect_count << ")  " ;
                //masks.print();

                masks.next();
        }

	if (w_is_not_deodhar == 0)
	{
        cout << "In lex order on reduced expressions: " << endl;
        for( map<string, map<string, int>, strCmp>::iterator iter = kl_basis_element.begin(); iter != kl_basis_element.end(); iter++ )
        {
                string t_name = (*iter).first;

		// if x is defined, then restrict to P_{x,w} entry.
		if (!(x_s == "") && !(t_name == x_name)) { continue; } 

                string q_polynomial;
                int started = 0;
                for( map<string, int, strCmp>::iterator iter2 = (*iter).second.begin(); iter2 != (*iter).second.end(); iter2++ )
                {
                        if (started == 1) { q_polynomial = q_polynomial + " + "; }

                        if ((*iter2).second == 1)
                        {
                          q_polynomial = q_polynomial + (*iter2).first;
                        }
                        else
                        {
                          stringstream s;
                          s << (*iter2).second;
                          q_polynomial = q_polynomial + s.str() + "." + (*iter2).first;
                        }

                        started = 1;
                }

		        // check w_s is the internal reduced experession, so masks match up.
                string statement = "P(" + w_name + "," + t_name + ") = " + q_polynomial;
                cout << statement;

		cout << endl;
	}
	}
}

