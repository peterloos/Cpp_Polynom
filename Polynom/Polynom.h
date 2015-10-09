class Polynom
{
private:
    double* m_coefficients;
    int     m_count;

public:
    // c'tors / d'tor
    Polynom  ();
    Polynom  (double coefficients[], int count);
    Polynom  (const Polynom&);
    ~Polynom ();

public:
    // getter
    int Rank () const;
    bool IsZero () const;

    // assignment operator
    Polynom& operator= (const Polynom& p);

    // unary mathematical operators
    friend Polynom operator+ (const Polynom& p);
    friend Polynom operator- (const Polynom& p);

    // binary mathematical operators
    friend Polynom operator+ (const Polynom& p1, const Polynom& p2);
    friend Polynom operator- (const Polynom& p1, const Polynom& p2);
    friend Polynom operator* (const Polynom& p1, const Polynom& p2);
    friend Polynom operator/ (const Polynom& p1, const Polynom& p2);
    friend Polynom operator% (const Polynom& p1, const Polynom& p2);

    // binary mathematical assignment operators
    friend Polynom& operator+= (Polynom& p1, const Polynom& p2);
    friend Polynom& operator-= (Polynom& p1, const Polynom& p2);
    friend Polynom& operator*= (Polynom& p1, const Polynom& p2);
    friend Polynom& operator/= (Polynom& p1, const Polynom& p2);
    friend Polynom& operator%= (Polynom& p1, const Polynom& p2);

    // comparison operators
    friend bool operator== (const Polynom& p1, const Polynom& p2);
    friend bool operator!= (const Polynom& p1, const Polynom& p2);
    friend bool operator<  (const Polynom& p1, const Polynom& p2);
    friend bool operator<= (const Polynom& p1, const Polynom& p2);
    friend bool operator>  (const Polynom& p1, const Polynom& p2);
    friend bool operator>= (const Polynom& p1, const Polynom& p2);

    // index operator
    double operator[] (double x);

    // function call operator
    double operator() (double x);

    // output
    friend ostream& operator<< (ostream& os, const Polynom& p);

private:
    // private helper operators
    friend Polynom operator* (const Polynom& p, double d);
    friend Polynom operator* (double d, const Polynom& p);

    // private helper methods
    void MultiplyX ();
    void MultiplyX (int k);

    // horner scheme
    double Compute (double x);
};