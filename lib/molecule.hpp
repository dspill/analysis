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

/** @file molecule.hpp 
 * The Molecule class allows for the modelling of chain molecules.
 * For a Molecule object only the reference to the first and last bead of the
 * molecule is stored in form of iterators.
 * It is based on a full vector of coordinates (and optionally a vector of
 * velocities) that may correspond to multiple molecules and is stored at
 * another place (typically in a Frame object).
 */
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

    /** Initialize Molecule with references to the start/end coordinates and
     * velocities.
     */
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

    /** Initialize Molecule with references to the start/end coordinates.
     */
    Molecule(r3dcit coordinates_begin, r3dcit coordinates_end) :
        _coordinates_begin(coordinates_begin), 
        _coordinates_end(coordinates_end) 
    {
        //std::cout << "Created " << *this << '\n';
        assert(coordinates_begin <= coordinates_end);
    }

    /** @return The number of beads in the Molecule. */
    size_t size() const
    {
        return std::distance(_coordinates_begin, _coordinates_end);
    }

    /** @returns Whether or not the Molecule stores velocities. */
    bool has_velocities() const
    {
        return std::distance(_velocities_begin, _velocities_end);
    }

    /** @return Whether or not there are any particles in the Molecule. */
    bool is_null() const {return !size();}

    /** 
     * @param[in] index Index within the Molecule of the particle in question.
     * @return Coordinate vector of the particle with the given index. */
    Real3D coordinate(size_t index) const
    {
        assert(index < size());
        return *(_coordinates_begin + index);
    }

    /** 
     * @param[in] index Index within the Molecule of the particle in question.
     * @return Velocity vector of the particle with the given index. */
    Real3D velocity(size_t index) const
    {
        assert(has_velocities());
        assert(index < size());
        return *(_velocities_begin + index);
    }

    /** 
     * @param[in] index Index within the Molecule of the bond in question.
     * @return Bond vector of the particle with the given index. */
    Real3D bond(const int index) const
    {
        assert(index < static_cast<int>(size()) - 1);
        return *(_coordinates_begin + index + 1)
            - *(_coordinates_begin + index);
    }

    /** 
     * @param[in] index Index within the Molecule of the bond in question.
     * @return Bond length of the particle with the given index. */
    double bond_length(const int index) const
    {
        //assert(bond(index).abs() < 2.);
        return bond(index).abs();
    }

    /** Test for any bonds that might be too large.  */
    bool consistent() const
    {
        bool c = true;
        for(size_t i = 0; i < size() - 1; ++i)
        {
            if(bond_length(i) > 1.733)
            {
                //std::cerr << "large bond (" << bond_length(i) << ") detected\n";
                c = false;
            }
        }
        return c;
    }

    // analysis 
    /** @return The center of mass vector. */
    Real3D center_of_mass() 
    {
        Real3D r_cm(0.);
        for(r3dcit cit = _coordinates_begin; cit != _coordinates_end; ++cit)
        {
            r_cm += *cit;
        }
        return r_cm / size();
    }

    // redundand
    Real3D rouse_mode_0() 
    {
        Real3D rm(0.);
        for(r3dcit cit = _coordinates_begin; cit != _coordinates_end; ++cit)
        {
            rm += *cit;
        }

        return rm / sqrt(size());
    }

    /** 
     * @param[in] p Index of the Rouse mode.
     * @return The corresponding Rouse mode. */
    Real3D rouse_mode(const size_t p) 
    {
        Real3D rm(0.);
        if(p == 0)
        {
            double pref = 1./sqrt(size());
            for(r3dcit cit = _coordinates_begin; cit != _coordinates_end; ++cit)
            {
                rm += pref * *cit;
            }
        }
        else
        {
            size_t i = 1;
            double pref_1 = sqrt(2./ size());
            double pref_2 = p * M_PI / size();
            for(r3dcit cit = _coordinates_begin; cit != _coordinates_end; ++cit)
            {
                rm += pref_1 * cos(pref_2 * (i - .5));
            }
        }
        return rm;
    }

    /** @return The squared radius of gyration. */
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

    /** @return The end-to-end vector. */
    Real3D end_to_end(){
        return *(_coordinates_end - 1) - *_coordinates_begin;
    }

    /** @return The squared end-to-end vector. */
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

    /** @return The average bond length. */
    double mean_bond_length() const
    {
        double mbl = 0.;
        for(size_t i = 0; i < size() - 1; ++i) mbl += bond_length(i);

        return mbl/(size() - 1);
    }

    /** @return The largest bond length. */
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
    /** Disable the assignment operator. */
    Molecule & operator=(const Molecule&) = delete; // disable assign

    /** Print Molecule info. */
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
