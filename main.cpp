/// ------------------------------------------
/// @file main.cpp
///
/// @brief Main starting point for Matrix testing
/// ------------------------------------------

#include <iostream>
#include <time.h>
#include <chrono>

#include "inc/Matrix.h"
#include "inc/Complex.h"
#include "inc/Poly.h"

using std::cout;
using std::endl;


class Timer
{
private:
	// Type aliases to make accessing nested type easier
	using Clock = std::chrono::steady_clock;
	using Second = std::chrono::duration<double, std::ratio<1> >;

	std::chrono::time_point<Clock> m_beg { Clock::now() };

public:
	void reset()
	{
		m_beg = Clock::now();
	}

	double elapsed() const
	{
		return std::chrono::duration_cast<Second>(Clock::now() - m_beg).count();
	}
};

// Generate a random nxn matrix
Matrix<double> gen_random_mat(size_t len, double lower, double upper);

int main()
{
    Matrix<double> test = gen_random_mat(20, 0, 10);

    Timer t;

    auto QR_set = test.qr_decompose();
    cout << "-----A------" << endl;
    cout << test << endl;

    cout << "-----Q------" << endl;
    cout << QR_set.first << endl;

    cout << "-----R------" << endl;
    cout << QR_set.second << endl;

    t.reset();
    cout << "-----QR-------" << endl;
    cout << test.inverse() << endl;
    cout << t.elapsed() << " sec" << endl;

    /*
    vec = {2,4,3};
    test.setCol(1, vec);
    col = test.getCol(1);
    cout << "Col: ";
    for (auto num : col)
    {
        cout << num << ", ";
    }
    cout << endl;

    cout << test << endl;

    cout << "----------------------" << endl;

    Vector<double> vect({2,3,4});
    Vector<double> vect2({5,6,7});

    cout << vect << endl;
    cout << vect*vect2 << endl;
    cout << vect.crossR3(vect2) << endl;

    cout << vect.magnitude() << endl;
    cout << vect2.magnitude() << endl;

    cout << vect.cosineAng(vect2) << endl;

    vect2 = vect;
    cout << vect2 << endl;
    */

    return 0;
}

Matrix<double> gen_random_mat(size_t len, double lower, double upper)
{
    Matrix<double> outMat(len, len);

    srandom(time(NULL));
    const long max_rand = 1000000L;

    for (size_t i = 0; i < len; i++)
    {
        for(size_t j = 0; j < len; j++)
        {
            double random_double = lower + (upper - lower) * (random() % max_rand) / upper;

            outMat.set(i,j, random_double);
        }
    }

    return outMat;
}