/// ------------------------------------------
/// @file main.cpp
///
/// @brief Main starting point for Matrix testing
/// ------------------------------------------

#include <iostream>

#include "inc/Matrix.h"

using std::cout;
using std::endl;

template<typename T>
void print(Matrix<T> mat)
{
    for (size_t i = 0; i < mat.getRowCount(); i++)
    {
        std::stringstream row;
        for (size_t j = 0; j < mat.getColCount(); j++)
        {
            row << mat.get(i,j);
            if (j+1 != mat.getColCount())
            {
                row << ", ";
            }
        }

        cout << row.str() << endl;
    }
}

int main()
{
    auto mat = Matrix<double>({{5,6,9},{2,1,6},{1,2,3}});

    print(mat);
    cout << "-----------------" << endl;
    print(mat.adjoint());

    return 0;
}