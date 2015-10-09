#include <iostream>
using namespace std;

#include "Polynom.h"

// c'tors
Polynom::Polynom ()
{
    // default c'tor: create zero polynom
    m_count = 1;
    m_coefficients = new double[1];
    m_coefficients[0] = 0.0;
}

Polynom::Polynom (const Polynom& p)
{
    // copy c'tor
    m_count = p.m_count;
    m_coefficients = new double[p.m_count];
    for (int i = 0; i < m_count; i++)
        m_coefficients[i] = p.m_coefficients[i];
}

Polynom::Polynom (double coefficients[], int count)
{
    // remove leading zeros
    int max = count - 1;
    while (max > 0 && coefficients[max] == 0)
        max --;

    m_count = max + 1;
    m_coefficients = new double[m_count];

    for (int i = 0; i < m_count; i ++)
        m_coefficients[i] = coefficients[i];
}

// d'tor
Polynom::~Polynom ()
{
    delete[] m_coefficients;
}

// getter
int Polynom::Rank () const
{
    return m_count - 1;
}

bool Polynom::IsZero () const
{
    return m_count == 1 && m_coefficients[0] == 0;
}

// unary mathematical operators + and -
Polynom operator+ (const Polynom& p)
{
    return Polynom (p);
}

Polynom operator- (const Polynom& p)
{
    Polynom tmp (p);
    for (int i = 0; i < tmp.m_count; i++)
        tmp.m_coefficients[i] *= -1.0;
    return tmp;
}

// binary mathematical operators +, -, *, / and %
Polynom operator+ (const Polynom& p1, const Polynom& p2)
{
    int count = (p1.m_count <= p2.m_count) ? p2.m_count : p1.m_count;

    // create array of coefficients
    double* coefficients = new double[count];
    for (int i = count - 1; i >= 0; i --)
    {
        double coeff = 0.0;
        if (i < p1.m_count)
            coeff += p1.m_coefficients[i];
        if (i < p2.m_count)
            coeff += p2.m_coefficients[i];
        coefficients[i] = coeff;
    }

    // create result polynom
    Polynom tmp (coefficients, count);

    // release temporary coefficients array
    delete[] coefficients;

    return tmp;
}

Polynom operator- (const Polynom& p1, const Polynom& p2)
{
    Polynom tmp = - p2;
    tmp = p1 + tmp;
    return tmp;
}

Polynom operator* (const Polynom &p1, const Polynom &p2)
{
    // create array of coefficients
    int count = p1.m_count + p2.m_count -1;
    double* coefficients = new double[count];

    // clear coefficients array
    for (int i = 0; i < count; i ++)
        coefficients[i] = 0.0;

    // compute coefficients of polynom product
    for (int i = p1.m_count - 1; i >= 0; i--)
        for (int j = p2.m_count - 1; j >= 0; j--)
            coefficients[i + j] +=
                p1.m_coefficients[i] * p2.m_coefficients[j];

    // create result polynom
    Polynom result (coefficients, count);

    // delete temporary array of coefficients
    delete[] coefficients;

    return result;
}

Polynom operator/ (const Polynom& p1, const Polynom& p2)
{
    if (p1.m_count < p2.m_count)  // degree of numerator polynom is less than
        return Polynom();         // degree of denominator polynom

    // need copies of arguments
    Polynom tmp1 = p1;
    Polynom tmp2 = p2;

    // create coefficients array of result polynom
    int count = p1.m_count - p2.m_count + 1;
    double* rescoeff = new double[count];

    // clear coefficients array
    for (int i = 0; i < count; i ++)
        rescoeff[i] = 0.0;

    // apply algorithm of polynom division
    for (int i = count - 1; i >= 0 ; i--)
    {
        // premature end of division reached (comparing degrees)
        if (tmp1.m_count < p2.m_count)
            break;

        // calculate next coefficient of result polynom
        double coeff =
            tmp1.m_coefficients[tmp1.m_count - 1] /
            tmp2.m_coefficients[tmp2.m_count - 1];

        // multiply denominator polynom with coefficient
        tmp2 = tmp2 * coeff;

        // calculate difference of ranks
        int diffRank = tmp1.m_count - p2.m_count;

        // multiply denominator polynom with one ore more 'x'
        tmp2.MultiplyX(diffRank);

        // subtract denominator polynom from numerator polynom
        tmp1 = tmp1 - tmp2;

        // poke calculated coefficient into result polynom
        rescoeff[diffRank] = coeff;

        // restore denominator polynom
        tmp2 = p2;
    }

    return Polynom(rescoeff, count); 
}

Polynom operator% (const Polynom& p1, const Polynom& p2)
{
    return p1 - (p1 / p2) * p2;
}

// binary mathematical assignment operators +=, -=, *=, /= and %=
Polynom& operator+= (Polynom& p1, const Polynom& p2)
{
    p1 = p1 + p2;
    return p1;
}

Polynom& operator-= (Polynom& p1, const Polynom& p2)
{
    p1 = p1 - p2;
    return p1;
}

Polynom& operator*= (Polynom& p1, const Polynom& p2)
{
    p1 = p1 * p2;
    return p1;
}

Polynom& operator/= (Polynom& p1, const Polynom& p2)
{
    p1 = p1 / p2;
    return p1;
}

Polynom& operator%= (Polynom& p1, const Polynom& p2)
{
    p1 = p1 % p2;
    return p1;
}




// horner scheme
double Polynom::Compute (double x)
{
    double y = m_coefficients[m_count - 1];        
    for (int i = m_count - 2; i >= 0; i--) 
        y = m_coefficients[i] + y * x;
    return y;
}

// apply horner scheme, using array subscripting operator
double Polynom::operator[](double x)
{
    return Compute (x);
}

// apply horner scheme, using function call operator
double Polynom::operator() (double x)
{
    return Compute (x);
}

// comparison operators
bool operator== (const Polynom& p1, const Polynom& p2)
{
    if (p1.m_count != p2.m_count)
        return false;

    for (int i = 0; i < p1.m_count; i ++)
        if (p1.m_coefficients[i] != p2.m_coefficients[i])
            return false;

    return true;
}

bool operator!= (const Polynom& p1, const Polynom& p2)
{
    return ! (p1 == p2);
}

bool operator<  (const Polynom& p1, const Polynom& p2)
{
    if (p1.m_count < p2.m_count)
        return true;

    if (p1.m_count > p2.m_count)
        return false;

    for (int i = p1.m_count - 1; i >= 0; i --)
    {
        if (p1.m_coefficients[i] < p2.m_coefficients[i])
            return true;
        if (p1.m_coefficients[i] > p2.m_coefficients[i])
            return false;
    }

    return false;
}

bool operator<= (const Polynom& p1, const Polynom& p2)
{
    return (p1 < p2) || (p1 == p2);
}

bool operator>  (const Polynom& p1, const Polynom& p2)
{
    return ! (p1 <= p2);
}

bool operator>= (const Polynom& p1, const Polynom& p2)
{
    return ! (p1 < p2);
}

// output
ostream& operator<< (ostream& os, const Polynom& p)
{
    for (int i = p.m_count - 1; i >= 0; i --)
    {
        if (p.m_coefficients[i] == 0.0 && p.m_count > 1)
            continue;

        if (i < p.m_count - 1)
            if (p.m_coefficients[i] > 0.0)
                os << '+';
        
        if (p.m_coefficients[i] != 1 && p.m_coefficients[i] != -1)
        {
            os << p.m_coefficients[i];
        }
        else if (i > 0)
        {
            if (p.m_coefficients[i] == -1)
                os << '-';
        }
        else
            os << p.m_coefficients[i];

        if (i > 0)
            os << 'x' << '^' << i;
    }

    return os;
}

// assignment operator
Polynom& Polynom::operator= (const Polynom& p)
{
    if (this == &p)
        return *this;

    // release left side
    delete[] m_coefficients;

    // copy right side
    m_count = p.m_count;
    m_coefficients = new double[p.m_count];
    for (int i = 0; i < m_count; i++)
        m_coefficients[i] = p.m_coefficients[i];

    return *this;
}

// private helper operators
Polynom operator* (double d, const Polynom& p)
{
    Polynom q = p;
    for (int i = 0; i < p.m_count; i ++)
        q.m_coefficients[i] *= d;
    return q;
}

Polynom operator* (const Polynom& p, double d)
{
    return d * p;
}

// private helper methods
void Polynom::MultiplyX ()
{
    // create new array of coefficients
    double* tmp = new double[m_count + 1];

    // compute new coefficients
    tmp[0] = 0.0;
    for (int i = 1; i <= m_count; i ++)
        tmp[i] = m_coefficients[i - 1];

    // switch coefficients buffer
    delete[] m_coefficients;
    m_coefficients = tmp;
    m_count ++;
}

void Polynom::MultiplyX (int k)
{
    for (int i = 0; i < k; i ++) 
        MultiplyX ();
}

