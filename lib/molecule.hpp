#ifndef _MOLECULE_HPP
#define _MOLECULE_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <omp.h>
#include <assert.h>
#include <random>
#include "timeseries.hpp"
#include "Real3D.hpp"
#include "couf.hpp"

class Molecule
{
    using r3dcit = std::vector<Real3D>::const_iterator;

    private:

    const r3dcit _coordinates_begin;
    const r3dcit _coordinates_end;

    const r3dcit _velocities_begin;
    const r3dcit _velocities_end;

    public:
    //Molecule(const Molecule&) = delete; // disable copy constructor

    Molecule(r3dcit coordinates_begin, r3dcit coordinates_end, 
            r3dcit velocities_begin, r3dcit velocities_end) :
        _coordinates_begin(coordinates_begin), 
        _coordinates_end(coordinates_end),
        _velocities_begin(velocities_begin),
        _velocities_end(velocities_end) 
    {
        //std::cout << "Created " << *this << '\n';
        assert(coordinates_begin <= coordinates_end);
        assert(velocities_begin <= velocities_end);
    }

    Molecule(r3dcit coordinates_begin, r3dcit coordinates_end) :
        _coordinates_begin(coordinates_begin), 
        _coordinates_end(coordinates_end) 
    {
        //std::cout << "Created " << *this << '\n';
        assert(coordinates_begin <= coordinates_end);
    }

    size_t size() const
    {
        return std::distance(_coordinates_begin, _coordinates_end);
    }

    bool has_velocities() const
    {
        return std::distance(_velocities_begin, _velocities_end);
    }

    bool is_null() const {return !size();}

    Real3D coordinate(size_t index) const
    {
        assert(index < size());
        return *(_coordinates_begin + index);
    }

    Real3D velocity(size_t index) const
    {
        assert(index < size());
        return *(_velocities_begin + index);
    }

    Real3D bond(const int index) const
    {
        assert(index < static_cast<int>(size()) - 1);
        return *(_coordinates_begin + index + 1)
            - *(_coordinates_begin + index);
    }

    double bond_length(const int index) const
    {
        return bond(index).abs();
    }

    // analysis 
    Real3D center_of_mass() 
    {
        Real3D r_cm(0.);
        for(r3dcit cit = _coordinates_begin; cit != _coordinates_end; ++cit)
        {
            r_cm += *cit;
        }
        return r_cm / size();
    }

    double radius_of_gyration_squared() 
    {
        Real3D r_cm = center_of_mass();
        double r_gyr_sq = 0.;

        for(r3dcit cit = _coordinates_begin; cit != _coordinates_end; ++cit)
        {
            r_gyr_sq += (*cit - r_cm).sqr();
        }
        return r_gyr_sq / size();
    }

    Real3D end_to_end()  {return *(_coordinates_end - 1) - *_coordinates_begin;}

    double end_to_end_squared() {return end_to_end().sqr();}

    double mean_squared_displacement(const Molecule & earlier_molecule) const
    {
        /* computes square displacement relative to earlier_molecule */
        assert(size() == earlier_molecule.size());
        double msd = 0.;

        for(size_t ipart = 0; ipart < size(); ++ipart)
        {
            msd += (coordinate(ipart) -
                    earlier_molecule.coordinate(ipart)).sqr();
        }
        return msd / size();
    }

    double mean_bond_length() const
    {
        double mbl = 0.;
        for(size_t i = 0; i < size() - 1; ++i) mbl += bond_length(i);

        return mbl/(size() - 1);
    }

    double max_bond_length() const
    {
        double largest = 0.;
        for(size_t i = 0; i < size() - 1; ++i)
        {
            if(bond_length(i) > largest) largest = bond_length(i);
        }

        return largest;
    }

    // operators
    Molecule & operator=(const Molecule&) = delete; // disable assign

    friend std::ostream & operator<<(std::ostream & os, 
            const Molecule & molecule)
    {
        if(molecule.is_null()) 
        {
            os << "Empty molecule\n";
        }
        else
        {
            os << "Molecule with " << molecule.size() << " particles.\n";
        }
        return os;  
    }
};
#endif
