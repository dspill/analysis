/*
   Copyright (C) 2012,2013
   Max Planck Institute for Polymer Research
   Copyright (C) 2008,2009,2010,2011
   Max-Planck-Institute for Polymer Research & Fraunhofer SCAI

   This file is part of ESPResSo++.

   ESPResSo++ is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   ESPResSo++ is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>. 
   */

#ifndef _REAL3D_HPP
#define _REAL3D_HPP

#include <cmath>
#include <ostream>
#include <assert.h>

/** @file Real3D.hpp */

/** This class is based on the Real3D class provided by ESPResSo++.
 * It models a 3-dimensional vector and includes among others some useful
 * operations such as addition and scalar product between vectors.
 */
class Real3D {
    double data[3];

    private:
    //friend class boost::serialization::access;
    //template<class Archive>
    //void serialize(Archive & ar, const unsigned int version) {
    //for(int i = 0; i < 3; ++i) ar & data[i];      
    //}
    //TODO noexcept move constructor

    public:

    typedef double* iterator;

    Real3D(std::initializer_list<double> il);
    Real3D() noexcept;
    Real3D(double v) noexcept; 
    Real3D(double x, double y, double z) noexcept;
    Real3D(const Real3D& v) noexcept;
    Real3D(const double v[3]) noexcept;

    // assignment is not the same as initialization
    Real3D& operator=(const Real3D& v) noexcept;

    double& operator[](int i);
    const double& operator[] (int i) const;

    double& at(int i);
    const double& at(int i) const;

    void setItem(int i, double v);
    double getItem(int i) const;

    // unary operators
    Real3D& operator+=(const Real3D& v);
    Real3D& operator-=(const Real3D& v);
    Real3D& operator*=(const double v);
    Real3D& operator/=(const double v);

    // bool operators
    bool operator==(const Real3D& v) const;
    bool operator!=(const Real3D& v) const;

    // elementwise binary operators
    Real3D operator+ (const Real3D &v) const;
    Real3D operator- (const Real3D &v) const;
    Real3D operator* (double v) const;
    Real3D operator/ (double v) const;
    /** Cross product of two Real3D. */
    Real3D cross(const Real3D& v) const;

    // binary dot product
    double operator* (const Real3D& v) const;

    double sqr() const;
    double abs() const;

    // STL iterator interface
    iterator begin();
    iterator end();

    const double* get() const { return data; }
    double* get() { return data; }

    static void registerPython();
};

//////////////////////////////////////////////////
// Global operators
Real3D operator*(double s, const Real3D& v);
std::ostream &operator<<(std::ostream &out, const Real3D& v);

//////////////////////////////////////////////////
// INLINE IMPLEMENTATION
//////////////////////////////////////////////////

//////////////////////////////////////////////////
// Real3D
inline Real3D::Real3D() noexcept {}

inline Real3D::Real3D(double v) noexcept
{ data[0] = data[1] = data[2] = v; }

inline Real3D::Real3D(double x, double y, double z) noexcept {
    data[0] = x;
    data[1] = y;
    data[2] = z;
}

inline Real3D::Real3D(const Real3D &v) noexcept {
    for (int i = 0; i < 3; i++)
        data[i] = v[i];
}

inline Real3D::Real3D(const double v[3]) noexcept {
    for (int i = 0; i < 3; i++)
        data[i] = v[i];
}

inline Real3D::Real3D(std::initializer_list<double> il)
{
    assert(il.end() - il.begin() == 3);
    std::copy(il.begin(), il.end(), data);
}

inline Real3D &Real3D::operator=(const Real3D &v) noexcept {
    data[0] = v[0];
    data[1] = v[1];
    data[2] = v[2];
    return *this;
}

inline double &Real3D::operator[](int i) 
{ return data[i]; }    

inline const double &Real3D::operator[](int i) const
{ return data[i]; }    

inline double &Real3D::at(int i) {
    assert(i < 0 || i > 2);
    return (*this)[i];
}

inline const double &Real3D::at(int i) const {
    assert(i < 0 || i > 2);
    return (*this)[i];
}

inline void Real3D::setItem(int i, double v)
{ this->at(i) = v; }

inline double Real3D::getItem(int i) const
{ return this->at(i); }

// unary operators
inline Real3D& Real3D::operator+=(const Real3D &v)
{ for (int i = 0; i < 3; i++) data[i] += v.data[i]; return *this; }

inline Real3D& Real3D::operator-=(const Real3D &v)
{ for (int i = 0; i < 3; i++) data[i] -= v.data[i]; return *this; }

inline Real3D& Real3D::operator*=(const double v)
{ for (int i = 0; i < 3; i++) data[i] *= v; return *this; }

inline Real3D& Real3D::operator/=(const double v) { 
    double v_1 = 1.0/v;
    for (int i = 0; i < 3; i++) 
        data[i] *= v_1; 
    return *this;
}

// bool operators
inline bool Real3D::operator==(const Real3D &v) const {
    return 
        (data[0] == v.data[0]) &&
        (data[1] == v.data[1]) &&
        (data[2] == v.data[2]);
}

inline bool Real3D::operator!=(const Real3D &v) const 
{ return ! (*this == v); }

// elementwise binary operators
inline Real3D Real3D::operator+ (const Real3D &v) const
{ return Real3D(*this) += v; }

inline Real3D Real3D::operator- (const Real3D &v) const
{ return Real3D(*this) -= v; }

inline Real3D Real3D::operator* (double v) const
{ return Real3D(*this) *= v; }

inline Real3D Real3D::operator/ (double v) const
{ return Real3D(*this) /= v; }

// binary dot product
inline double Real3D::operator* (const Real3D& v) const
{ return data[0]*v.data[0] + data[1]*v.data[1] + data[2]*v.data[2]; }

/** Cross product of two Real3D. */
inline Real3D Real3D::cross(const Real3D& v) const {
    return Real3D(data[1]*v[2] - data[2]*v[1],
            data[2]*v[0] - data[0]*v[2],
            data[0]*v[1] - data[1]*v[0]);
}

inline double Real3D::sqr() const
{ return data[0]*data[0] + data[1]*data[1] + data[2]*data[2]; }

inline double Real3D::abs() const
{ return sqrt(sqr()); }

inline Real3D::iterator Real3D::begin() { return data; }
inline Real3D::iterator Real3D::end() { return data+3; }

//////////////////////////////////////////////////
// Global operators
inline Real3D operator*(double s, const Real3D &v) 
{ return Real3D(v)*s; }

inline std::ostream &operator<<(std::ostream &out, 
        const Real3D &v) {
    return out << v[0] << ' ' << v[1] << ' ' << v[2];
}

//}
#endif
