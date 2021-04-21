#ifndef _LB_FRAME_HPP
#define _LB_FRAME_HPP

#include <iostream>
#include <iterator>
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

#include "Real3D.hpp"
#include "vector.hpp"
/** @file lb_frame.hpp */

/** Struct storing an index triple corresponding to some lattice site. */
struct Index
{
    size_t x;
    size_t y;
    size_t z;
};

/** Class that stores a D3Q19 Lattice-Bolzmann configuration.
 * For each lattice site only the density and the flow velocity is stored.
 */
class LB_Frame 
{
    private:
        std::vector<double> _densities;
        std::vector<Real3D> _velocities;
        std::vector<size_t> _box{0, 0, 0};

        /**
         * Read LB configuration from file stream.
         */
        void pure_read(std::ifstream & stream)
        {
            size_t index;
            Index i;
            double pop;
            Real3D vel;
            std::string str;
            while(std::getline(stream, str))
            {
                // get lattice point
                std::string buf;            // Have a buffer string
                std::stringstream ss(str);  // Insert the string into a stream

                ss >> buf;
                i.x = stoi(buf);
                ss >> buf;
                i.y = stoi(buf);
                ss >> buf;
                i.z = stoi(buf);

                index = get_index(i);
                assert(_densities.size() > index);
                assert(_velocities.size() > index);

                ss >> buf;
                pop = stof(buf);
                ss >> buf;
                vel[0] = stof(buf);
                ss >> buf;
                vel[1] = stof(buf);
                ss >> buf;
                vel[2] = stof(buf);

                _densities[index] = pop;
                _velocities[index] = vel;
            }
        }

        /**
         * Read LB configuration from file.
         * The file needs the following column structure:
         * i_x i_y i_z rho v_x v_y vz,
         * where the i are the lattice indexes, rho is the density at the
         * lattice site and v_j are the velocity components at the lattice
         * site.
         */
        void pure_read(const std::string & filename)
        {
            std::ifstream stream(filename);
            if (!stream)
            {
                throw std::runtime_error("Could not open input file\n");
            }
            pure_read(stream);
            stream.close();
        }

    public:
        /* Constructors */
        LB_Frame(){}

        /**
         * Construct LB_Frame from either file or file stream.
         */
        template <typename T, typename U>
            LB_Frame(T input, U box)
            {
                set_box(box);
                pure_read(input);
            }

        /**
         * Construct empty LB_Frame with given box size.
         */
        template <typename T>
            LB_Frame(T box)
            {
                set_box(box);
            }

        /** @returns Box dimensions. */
        std::vector<size_t> box() const
        {
            return _box;
        }

        /** @returns Box length in given direction. */
        size_t box(size_t i) const
        {
            assert(i < 3);
            return _box[i];
        }

        /** @returns Box length in given direction. */
        template <typename T>
            void add(T input)
            {
                if(is_null()) 
                {
                    throw std::runtime_error("You can not add to NULL frame");
                }
                else pure_read(input);
            }

        /** Set cubic box. */
        void set_box(size_t L)
        {
            _box[0] = L;
            _box[1] = L;
            _box[2] = L;
            _densities.resize(size());
            _velocities.resize(size());
        }

        /** Set box. */
        void set_box(std::vector<size_t> box)
        {
            assert(box.size() == 3);
            _box = box;
            _densities.resize(size());
            _velocities.resize(size());
        }

        size_t get_index(size_t index) const
        {
            return index;
        }

        /**
         * Get linear lattice index from Cartesian index triplet.
         * @param[in] i_x x-coordinate.
         * @param[in] i_y y-coordinate.
         * @param[in] i_z z-coordinate.
         * @return linear index.
         */
        size_t get_index(size_t i_x, size_t i_y, size_t i_z) const
        {
            return i_z + box(2) * (i_y + box(1) * i_x);
        }

        size_t get_index(Index i) const
        {
            return get_index(i.x, i.y, i.z);
        }

        /**
         * Compute Cartesian indexes from linear index via modulo operations.
         */
        Index to_Index(size_t index)
        {
            size_t z = index % box(2);
            index -= z;
            size_t y = (index / box(2)) % box(1);
            index -=  y * box(2);
            size_t x = index / (box(1) * box(2));

            return Index{x, y, z};
        }

        /** @return Number of lattice sites. */
        size_t size() const
        {
            return box(0) * box(1) * box(2);
        }

        /** @return true if box is zero. */
        bool is_null() const
        {
            return size() == 0;
        }

        /* analysis methods */
        /** @return Density at given lattice index. */
        template <typename T>
            double density(T index) const
            {
                return _densities[get_index(index)];
            }

        /** @return Velocity at given lattice index. */
        template <typename T>
            Real3D velocity(T index) const
            {
                return _velocities[get_index(index)];
            }

        /** @return Momentum density at given lattice index. */
        template <typename T>
            Real3D momentum_density(T index) const
            {
                return density(index) * velocity(index);
            }

            Real3D total_momentum_density() const
            {
                Real3D total{0.};
                for(size_t i = 0; i < size(); ++i)
                    total += momentum_density(i);

                return total;
            }


        /** @return Total density. */
        double total_density() const
        {
            return std::accumulate(_densities.cbegin(), _densities.cend(), 0.0);
        }

        /** @return Average density. */
        double mean_density() const
        {
            return total_density() / size();
        }

        /** @return Print configuration info. */
        friend std::ostream & operator<<(std::ostream & os, const LB_Frame & f)
        {
            if(f.is_null())
            {
                os << "Null frame\n";
            }
            else
            {
                os << "Frame with " << f.size() << " sites. Box: "; 
                os << f.box(0) << ' ';
                os << f.box(1) << ' ';
                os << f.box(2);
                os << '\n';
            }
            return os;
        }
};

#endif
