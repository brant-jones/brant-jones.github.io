#include "CoxeterSystem.h"


CoxeterSystem::CoxeterSystem(const CoxeterSystem& cs)
{
	size = cs.size;
	coxeter_matrix = new int[size * size];
	for (int i = 0; i < size*size; i++)
	{
		coxeter_matrix[i] = cs.coxeter_matrix[i];
	}

	automorphism_group_size = cs.automorphism_group_size;
	automorphism_group = new int[automorphism_group_size * size];
	for (int i = 0; i < automorphism_group_size*size; i++)
	{
		automorphism_group[i] = cs.automorphism_group[i];
	}
}

CoxeterSystem::CoxeterSystem(int n, int cm[], int ags, int ag[])
{
	size = n;
	coxeter_matrix = new int[size * size];
	for (int i = 0; i < size*size; i++)
	{
		coxeter_matrix[i] = cm[i];
	}

	automorphism_group_size = ags;
	automorphism_group = new int[automorphism_group_size * size];
	for (int i = 0; i < ags*size; i++)
	{
		automorphism_group[i] = ag[i];
	}
}

int CoxeterSystem::equals(CoxeterSystem* cs)
{
		if (cs->size != this->size) { return 0; }
		int ret = 1;
		for (int i = 0; i < size*size; i++)
		{
			if (cs->coxeter_matrix[i] != this->coxeter_matrix[i]) { ret = 0; }
		}
		return ret;
}


void CoxeterSystem::print_matrix()
{
		for (int i = 0; i < size; i++)
		{ 
			for (int j = 0; j < size; j++)
			{
				cout << get_exponent(i,j) << " ";
			}
			cout << endl;
		}
}

int CoxeterSystem::get_exponent(int x, int y)
{
	return coxeter_matrix[ (size * x) + y ];
}

int CoxeterSystem::get_automorphism(int i, int j)
{
	return automorphism_group[ (size * i) + j ];
}

int CoxeterSystem::right_multiply(int one_line[], int i)
{
	int t = one_line[i];
	one_line[i] = one_line[i+1];
	one_line[i+1] = t;
	return 0;
}

CoxeterSystem::~CoxeterSystem()
{
	delete[] coxeter_matrix;
	delete[] automorphism_group;
}


/////////////////////////////////////////////////////////////////////

int TypeDCoxeterSystem::right_multiply(int one_line[], int i)
{
	if (i > 0)
	{
	int t = one_line[i-1];
	one_line[i-1] = one_line[i];
	one_line[i] = t;
	}
	else if (i == 0)
	{
	int t = one_line[i];
	one_line[i] = one_line[i+1];
	one_line[i+1] = t;
	one_line[i] = 0 - one_line[i];
	one_line[i+1] = 0 - one_line[i+1];
	}

	return 0;
}


