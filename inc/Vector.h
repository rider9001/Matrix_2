/// ------------------------------------------
/// @file Vector.h
///
/// @brief Header/Source file for vector object
///
/// @note Must implement all functions upon definition due to template format
/// ------------------------------------------
#pragma once

#include <stdexcept>
#include <sstream>
#include <string>
#include <cstring>
#include <cmath>
#include <vector>

#include "Complex_C.h"
#include "Complex_P.h"

/// @brief Templated class for storing, acsessing and performing operations on a vector of values
/// Vectors are fixed length, defined upon creation
template <typename T>
class Vector
{
    public:
        ///--------------------------------------------------------
        /// @brief Constructor for a vector object
        ///
        /// @param len length of the vector in items of the given type T, must be <1
        ///
        /// @throws std::invalid_argument if length < 1
        Vector(const size_t len)
        {
            if (len < 1)
            {
                throw std::invalid_argument("Cols/Rows of a matrix must be above 0");
            }

            m_length = len;
            vec_data = new T[m_length];
        }

        ///--------------------------------------------------------
        /// @brief Constructor using an initializer list to populate vector
        ///
        /// @param vecData array of data to populate vector with
        ///
        /// @throws std::invalid_argument if initializer list is empty
        Vector(const std::initializer_list<T>& vecData)
        {
            if(vecData.size() == 0)
            {
                throw std::invalid_argument("Initializer list empty");
            }

            m_length = vecData.size();
            vec_data = new T[m_length];

            size_t i = 0;
            for (auto num : vecData)
            {
                vec_data[i++] = num;
            }
        }

        ///--------------------------------------------------------
        /// @brief Constructor using a populated std vector
        ///
        /// @param inpVec input std vector to create math vector from
        Vector(const std::vector<T> inpVec)
        {
            if(inpVec.size() == 0)
            {
                throw std::invalid_argument("Initializer vector empty");
            }

            m_length = inpVec.size();
            vec_data = new T[m_length];

            for (size_t i = 0; i < inpVec.size(); i++)
            {
                vec_data[i] = inpVec[i];
            }
        }

        ///--------------------------------------------------------
        /// @brief Copy constructor
        ///
        /// @tparam T type stored by vector
        ///
        /// @param mat reference to copied vector
        Vector<T>(Vector<T> const& vec)
        {
            m_length = vec.size();
            vec_data = new T[m_length];

            memcpy(vec_data, vec.get_data(), m_length * sizeof(T));
        }

        ///--------------------------------------------------------
        /// @brief Destructor
        ~Vector()
        {
            delete vec_data;
        }

        ///--------------------------------------------------------
        /// @brief Operator overload of +, implements vector-vector addition
        ///
        /// @param vec reference to rval matrix
        ///
        /// @return result of summed vectors
        Vector<T> operator+(Vector<T> const& vec)
        {
            if (vec.size() != m_length)
            {
                throw std::invalid_argument("Vector addition requires matricies of same dimensions");
            }

            Vector<T> outVec(m_length);

            for (size_t i = 0; i < m_length; i++)
            {
                outVec.get_data[i] = vec_data[i] + vec.get_data[i];
            }

            return outVec;
        }

        ///--------------------------------------------------------
        /// @brief Operator overload of +, implements vector-vector subtraction
        ///
        /// @param vec reference to rval matrix
        ///
        /// @return result of summed vectors
        Vector<T> operator-(Vector<T> const& vec)
        {
            if (vec.size() != m_length)
            {
                throw std::invalid_argument("Vector subtraction requires matricies of same dimensions");
            }

            Vector<T> outVec(m_length);

            for (size_t i = 0; i < m_length; i++)
            {
                outVec.get_data()[i] = vec_data[i] - vec.get_data()[i];
            }

            return outVec;
        }

        ///--------------------------------------------------------
        /// @brief Operator overload of *, implements vector dot product
        ///
        /// @param vec vector to dot product with
        ///
        /// @return dot product of the two vectors
        T operator*(Vector<T> const& vec)
        {
            if (vec.size() != m_length)
            {
                throw std::invalid_argument("Dot product requires vectors of same size");
            }

            T sum = 0;
            for (size_t i = 0; i < m_length; i++)
            {
                sum += vec_data[i] * vec.get_data()[i];
            }

            return sum;
        }

        ///--------------------------------------------------------
        /// @brief Cross product of two vectors in R3 space
        ///
        /// @param vec R3 vector to cross product
        ///
        /// @return cross product in R3
        Vector<T> crossR3(Vector<T> const& vec)
        {
            if (m_length != 3 or vec.size() != 3)
            {
                throw std::invalid_argument("R3 cross product vectors must both be 3 elements");
            }

            Vector<T> outVec(3);
            outVec.set(0, get(1)*vec.get(2) - get(2)*vec.get(1));
            outVec.set(1, get(2)*vec.get(0) - get(0)*vec.get(2));
            outVec.set(2, get(0)*vec.get(1) - get(1)*vec.get(1));

            return outVec;
        }

        ///--------------------------------------------------------
        /// @brief Operator overload of *, implements scalar multiplication of vector
        ///
        /// @param num scalar to multiply by
        ///
        /// @return multiplied vector
        Vector<T> operator*(T const& num)
        {
            Vector<T> outVec(m_length);

            for (size_t i = 0; i < m_length; i++)
            {
                outVec.get_data()[i] = vec_data[i] * num;
            }

            return outVec;
        }

        ///--------------------------------------------------------
        /// @brief Operator overload of /, implements scalar division of vector
        ///
        /// @param num scalar to divide by
        ///
        /// @return multiplied vector
        Vector<T> operator/(T const& num)
        {
            Vector<T> outVec(m_length);

            for (size_t i = 0; i < m_length; i++)
            {
                outVec.get_data()[i] = vec_data[i] / num;
            }

            return outVec;
        }

        ///--------------------------------------------------------
        /// @brief Operator overload for ==, tests equivelance between two vectors
        ///
        /// @param vec vector to compare to
        ///
        /// @return are two vectors exactly identical?
        bool operator==(Vector<T> const& vec)
        {
            if (m_length != vec.size())
            {
                return false;
            }

            for(size_t i = 0; i < m_length; i++)
            {
                if (vec_data[i] != vec.get_data()[i])
                {
                    return false;
                }
            }

            return true;
        }

        ///--------------------------------------------------------
        /// @brief Operator overload for !=, tests non-equivelance between two vectors
        ///
        /// @param vec vector to compare to
        ///
        /// @return are two vectors not exactly identical?
        bool operator!=(Vector<T> const& vec)
        {
            return !(*this == vec);
        }

        ///--------------------------------------------------------
        /// @brief finds magnitude for the vector
        /// e.g: sqrt(a^2 + b^2 + c^2...)
        ///
        /// @return magnitude of the vector
        T magnitude() const
        {
            T sum = 0;
            for (size_t i = 0; i < m_length; i++)
            {
                // do not replace with pow, need to account for custom types
                sum += get_data()[i] * get_data()[i];
            }

            // Dynamically chooses which routine to implement depending on complex type
            if constexpr (std::is_same_v<T, Complex_C_t>)
            {
                sum = powReal(sum, 0.5);
            }
            else if constexpr (std::is_same_v<T, Complex_P_t>)
            {
                sum = powReal(sum, 0.5);
            }
            else
            {
                sum = std::sqrt(sum);
            }

            return sum;
        }

        ///--------------------------------------------------------
        /// @brief Calculates the cos of angle between two vectors
        ///
        /// @param vec vector to calculate angle between
        ///
        /// @return cosine of angle to vec
        T cosineAng(Vector<T> const& vec)
        {
            return (*this * vec) / (magnitude() * vec.magnitude());
        }

        ///--------------------------------------------------------
        /// @brief finds the scalar distance from another vector
        ///
        /// @param vec direction to find scalar distance from
        ///
        /// @return distance of vector from input direction
        T scalar_in_direction(const Vector<T>& vec)
        {
            return (*this * vec) / vec.magnitude();
        }

        ///--------------------------------------------------------
        /// @brief Calculates the normailzed vector
        ///
        /// @return normalized vector
        Vector<T> normalize() const
        {
            return *this / magnitude();
        }

        ///--------------------------------------------------------
        /// @brief returns value at the index given
        ///
        /// @param index index to read
        ///
        /// @return value at index
        ///
        /// @throws invalid_argument if index is out of bounds
        T get(const size_t& index) const
        {
            if (!(index < m_length))
            {
                std::stringstream err;
                err << "Index " << index << ", is not in bounds (" << m_length - 1 << ")";
                throw std::invalid_argument(err.str());
            }

            return vec_data[index];
        }

        ///--------------------------------------------------------
        /// @brief sets the value at index to val
        ///
        /// @param index index to set value
        /// @param val value to set
        ///
        /// @throws invalid_argument if index is out of bounds
        void set(const size_t& index, const T& val)
        {
            if (!(index < m_length))
            {
                std::stringstream err;
                err << "Index " << index << ", is not in bounds (" << m_length - 1 << ")";
                throw std::invalid_argument(err.str());
            }

            vec_data[index] = val;
        }

        ///--------------------------------------------------------
        /// @brief Returns the number of elements in the vector
        ///
        /// @return vector length
        size_t size() const
        {
            return m_length;
        }

        ///--------------------------------------------------------
        /// @brief Returns a reference of the internal array for the matrix
        ///
        /// @tparam T type of matrix
        ///
        /// @returns reference of internal array
        T* get_data() const
        {
            return vec_data;
        };

    private:
        /// @brief Stores matrix values
        T* vec_data;

        /// @brief Length of the vector in items
        size_t m_length;
};

///--------------------------------------------------------
/// @brief Overload of <<, used to convert vector into an output stream
///
/// @param os output stream
/// @param mat vector to push to output stream
///
/// @return output stream
template <typename T>
std::ostream& operator<<(std::ostream& os, const Vector<T>& vec)
{
    for (size_t i = 0; i < vec.size(); i++)
    {
        os << vec.get(i);

        if (i+1 != vec.size())
        {
            os << ", ";
        }
    }

    return os;
}