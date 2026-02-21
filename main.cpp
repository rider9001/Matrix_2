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
    Matrix<double> test = gen_random_mat(2, 0, 100);

    Timer t;

    cout << "-----A------" << endl;
    cout << test << endl;

    // cout << "-----Q------" << endl;
    // cout << QR_set.first << endl;

    // cout << "-----R------" << endl;
    // cout << QR_set.second << endl;

    // cout << "-----R^-1-----" << endl;
    // cout << QR_set.second.inverse() << endl;

    // t.reset();
    // cout << "-----QR-------" << endl;
    // auto inv = test.inverse();
    // cout << inv << endl;
    // cout << t.elapsed() * 1e6 << " micros" << endl;

    // cout << "-----proof-----" << endl;
    // cout << test % inv << endl;

    cout << "Egienvalues via QR convergence:" << endl;
    t.reset();
    auto eigenvalues = test.eigenvalues();
    for(auto val : eigenvalues)
    {
        cout << val << ", ";
    }
    cout << endl;

    cout << "Eigenvectors:" << endl;
    auto eigenvectors = test.eigenvectors();
    for (size_t i = 0; i < eigenvectors.size(); i++)
    {
        auto Vec = eigenvectors.at(i).internal_norm();
        for (size_t j = 0; j < Vec.size(); j++)
        {
            cout << Vec.get(j) << ", ";
        }
        cout << endl;
    }
    cout << t.elapsed() * 1e6 << " micros" << endl;

    return 0;
}

Matrix<double> gen_random_mat(size_t len, double lower, double upper)
{
    Matrix<double> outMat(len, len);

    srandom(time(NULL));
    const long max_rand = upper;

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