#include "Masks.h"

Masks::Masks(int length)
{
        this->length = length;
        size = (int) (length / MASK_TYPESIZE);
        size++;
        iterator = new MASK_DATATYPE[size];
        for (int i = 0; i < size; i++) { iterator[i] = (MASK_DATATYPE) 0; }
        pointer = 0;
        exhausted_flag = 0;

        //cout << "Mask initialized:  length " << length << ", size " << size << endl;
}

Masks::~Masks()
{
	delete[] iterator;
}

int Masks::exhausted()
{
        int done = 1;
	if (length == 0) { return 1; }
        for (int i = 0; i < length; i++)
        {
                if (get_value(i) == 0) { done = 0; }
        }

        if (done == 1)
        {
                exhausted_flag = 1;
                return 0;
        }

        return exhausted_flag;
}

int Masks::next()
{
        int i = 0;
        while (iterator[i] == ((MASK_DATATYPE) MASK_FLIP_AT))
        {
                iterator[i] = 0;
                i++;
        }
        iterator[i]++;

        return 0;
}

int Masks::proper()
{
	for (int i = 0; i < length; i++)
	{ if (get_value(i) == 0) { return 1; } }
	return 0;
}

int Masks::get_value(int position)
{
        int p = ((int) (position/MASK_TYPESIZE));
        int q = position - (p * (int) MASK_TYPESIZE);
        //cout << "get_value(" << position << "): p " << p << ", q " << q;
        //cout << "  iterator[0] = " << iterator[0] << ":: " << ((iterator[p] >> q) & 1) << endl;
        return ((iterator[p] >> q) & 1);
        return 0;
}

void Masks::print()
{
        for (int i = 0; i < length; i++)
        { cout << get_value(i); }

	//cout << " : " << iterator[0] << ", " << iterator[1]; // << ", " << iterator[2];

        cout << endl;
}

void Masks::sprint(string& s)
{
	s.append("( ");
        for (int i = 0; i < length; i++)
        { 
		char str[255];
		sprintf(str, "%d", get_value(i));
		s.append(str);
		s.append(" ");
	}
	s.append(")");
}

void Masks::test()
{
        cout << "abcdefghijklmnopqrstuvwxyz1234567890" << endl;
        //iterator[0] = (unsigned int) 4294967290U;  // 2^32
        //iterator[0] = (unsigned int)   2147483640U;  // 2^31
        //iterator[0] = (unsigned int)   1073741820U;  // 2^30
        //iterator[1] = (unsigned int)   15U;

	int c=0;
        while (exhausted() == 0)
        {
                print();
		c++;
                next();
        }
	cout << "Counted " << c << " masks." << endl;
}

