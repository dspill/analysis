#ifndef _VECTOR_HPP
#define _VECTOR_HPP

#include <vector>
#include <numeric>
#include <functional>
#include <iomanip>
#include <algorithm>


// arithmetic operations on vectors of T
// multiplication by scalar
    template<typename T>
inline std::vector<T>& operator*=(std::vector<T>& v, const double a)
{
    //std::transform(v.begin(), v.end(), v.begin(),
            //std::bind(std::multiplies<T>(), std::placeholders::_1, a));
    std::transform(v.begin(), v.end(), v.begin(), [a](auto& c){return c*a;});
    return v;
}
    template<typename T>
inline std::vector<T> operator*(std::vector<T> v, const double a)
{

    v *= a;
    return v;
}

// division by scalar
    template<typename T>
inline std::vector<T>& operator/=(std::vector<T>& v, const double a)
{
    v *= 1./a;
    return v;
}

    template<typename T>
inline std::vector<T> operator/(std::vector<T> v, const double a)
{
    v /= a;
    return v;
}

// scalar multiplication
    template<typename T>
inline double operator*(std::vector<T> v_1, const std::vector<T>& v_2)
{

    return std::inner_product(v_1.cbegin(), v_1.cend(), v_2.cbegin(), .0);
}

// addition
    template<typename T>
inline std::vector<T>& operator+=(std::vector<T>& v_1, const std::vector<T>& v_2)
{
    std::transform (v_1.begin(), v_1.end(), v_2.cbegin(), v_1.begin(), std::plus<T>());
    return v_1;
}
    template<typename T>
inline std::vector<T> operator+(std::vector<T> v_1, const std::vector<T>& v_2)
{
    v_1 += v_2;
    return v_1;
}

// substraction
    template<typename T>
inline std::vector<T>& operator-=(std::vector<T>& v_1, const std::vector<T>& v_2)
{
    std::transform (v_1.begin(), v_1.end(), v_2.cbegin(), v_1.begin(), std::minus<T>());
    return v_1;
}
    template<typename T>
inline std::vector<T> operator-(std::vector<T> v_1, const std::vector<T>& v_2)
{
    v_1 -= v_2;
    return v_1;
}

    template<typename T>
inline std::ostream & operator<<(std::ostream & os, const std::vector<T> & v)
{
    for(auto it = v.cbegin(); it != v.cend(); ++it)
        os << *it << ' ';
    return os;
}


// arithmetic operations on arrays of T
// multiplication by scalar
    template<typename T, size_t s>
inline std::array<T, s>& operator*=(std::array<T, s>& v, const double a)
{
    //std::transform(v.begin(), v.end(), v.begin(),
            //std::bind(std::multiplies<T>(), std::placeholders::_1, a));
    std::transform(v.begin(), v.end(), v.begin(), [a](auto& c){return c*a;});
    return v;
}
    template<typename T, size_t s>
inline std::array<T, s> operator*(std::array<T, s> v, const double a)
{

    v *= a;
    return v;
}

// division by scalar
    template<typename T, size_t s>
inline std::array<T, s>& operator/=(std::array<T, s>& v, const double a)
{
    v *= 1./a;
    return v;
}

    template<typename T, size_t s>
inline std::array<T, s> operator/(std::array<T, s> v, const double a)
{
    v /= a;
    return v;
}

// scalar multiplication
    template<typename T, size_t s>
inline double operator*(std::array<T, s> v_1, const std::array<T, s>& v_2)
{

    return std::inner_product(v_1.cbegin(), v_1.cend(), v_2.cbegin(), .0);
}

// addition
    template<typename T, size_t s>
inline std::array<T, s>& operator+=(std::array<T, s>& v_1, const std::array<T, s>& v_2)
{
    std::transform (v_1.begin(), v_1.end(), v_2.cbegin(), v_1.begin(), std::plus<T>());
    return v_1;
}
    template<typename T, size_t s>
inline std::array<T, s> operator+(std::array<T, s> v_1, const std::array<T, s>& v_2)
{
    v_1 += v_2;
    return v_1;
}

// substraction
    template<typename T, size_t s>
inline std::array<T, s>& operator-=(std::array<T, s>& v_1, const std::array<T, s>& v_2)
{
    std::transform (v_1.begin(), v_1.end(), v_2.cbegin(), v_1.begin(), std::minus<T>());
    return v_1;
}
    template<typename T, size_t s>
inline std::array<T, s> operator-(std::array<T, s> v_1, const std::array<T, s>& v_2)
{
    v_1 -= v_2;
    return v_1;
}

    template<typename T, size_t s>
inline std::ostream & operator<<(std::ostream & os, const std::array<T, s> & v)
{
    for(auto it = v.cbegin(); it != v.cend(); ++it)
        os << *it << ' ';
    return os;
}
#endif
