/// ------------------------------------------
/// @file main.cpp
///
/// @brief Main starting point for Matrix testing
/// ------------------------------------------

#include <iostream>
#include <time.h>
#include <chrono>

#include "inc/Matrix.h"
#include "inc/Vector.h"
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
Matrix<double> gen_random_mat(const size_t& len, const double& lower, const double& upper);

// Generate random n length vector
Vector<double> gen_random_vec(const size_t& len, const double& lower, const double& upper);

// Generate random number, between lower and upper
long gen_random(const double& lower, const double& upper);

int main()
{
    srandom(time(NULL));

    const size_t mat_size = 2;
    const double lower = 1;
    const double upper = 10;

    Matrix<double> test = gen_random_mat(mat_size, lower, upper);

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

    cout << "Eigenvalues via QR convergence:" << endl;
    t.reset();
    auto eigenvalues = test.eigenvalues();
    for(auto val : eigenvalues)
    {
        cout << val << ", ";
    }
    cout << endl;

    cout << "Eigenvectors: " << endl;
    std::vector<Vector<double>> e_vecs = test.eigenvectors();
    for(Vector<double> vec : e_vecs)
    {
        cout << vec << endl;
    }
    cout << t.elapsed() * 1e6 << " micros" << endl;

    return EXIT_SUCCESS;
}

long gen_random(const double& lower, const double& upper)
{
    const long max_rand = upper;
    return lower + (upper - lower) * (random() % max_rand) / upper;
}

Matrix<double> gen_random_mat(const size_t& len, const double& lower, const double& upper)
{
    Matrix<double> outMat(len, len);

    for (size_t i = 0; i < len; i++)
    {
        for(size_t j = 0; j < len; j++)
        {
            outMat.set(i,j, gen_random(lower, upper));
        }
    }

    return outMat;
}

Vector<double> gen_random_vec(const size_t& len, const double& lower, const double& upper)
{
    Vector<double> outVec(len);

    for (size_t i = 0; i < len; i++)
    {
        outVec.set(i, gen_random(lower, upper));
    }

    return outVec;
}