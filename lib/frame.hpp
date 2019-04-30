#ifndef _FRAME_HPP
#define _FRAME_HPP
#include "molecule.hpp"
#include "fftw3.h"
#include <typeinfo>
#include <memory>
#include <random>
#include <set>
#include <iomanip>

/** @file frame.hpp */

/**
 * Class that stores particle coordinates and velocities as well as box size.
 * If _particles_per_molecule == 0 we will consider all particles to be in one
 * molecule (zero molecules if system is empty).
 */
class Frame
{
    private:
        std::vector<Real3D> _coordinates; //!< vector of bead coordinates
        std::vector<Real3D> _velocities; //!< vector of bead velocities
        std::vector<int> _types; //!< vector of bead types

        Real3D _box{Real3D(0.)}; //!< 3d vector of box lengths
        size_t _particles_per_molecule{0}; //!< number of particles in one molecule

        auto c_begin() const {return _coordinates.begin();}
        auto c_end() const {return _coordinates.end();}
        auto v_begin() const {return _velocities.begin();}
        auto v_end() const {return _velocities.end();}
        auto t_begin() const {return _types.begin();}
        auto t_end() const {return _types.end();}

        /* folding */
        double fold(const double coordinate, const double box_length) const
        {
            if(box_length == 0.) return 0.;

            double folded_coordinate;

            folded_coordinate = fmod(coordinate, box_length);
            folded_coordinate += (coordinate < 0.) * box_length;

            return folded_coordinate;
        }

        Real3D fold(const Real3D coordinate) const
        {
            double x, y, z;

            x = fold(coordinate[0], _box[0]);
            y = fold(coordinate[1], _box[1]);
            z = fold(coordinate[2], _box[2]);

            return Real3D(x, y, z);
        }

        Real3D fold_2d(const Real3D coordinate) const
        {
            double x, y, z;
            x = fold(coordinate[0], _box[0]);
            y = fold(coordinate[1], _box[1]);
            z = coordinate[2];

            return Real3D(x, y, z);
        }

        double minimal_distance(const double c1, const double c2,
                const double box_length) const
        {
            /* compute the minimal distance between two particles along one
             * dimension while considering the periodic images of the box */
            double fc1 = fold(c1, box_length);
            double fc2 = fold(c2, box_length);

            double dist = fc2 - fc1;
            dist = dist - box_length * round(dist/box_length);

            return dist;
        }

        Real3D minimal_distance(const Real3D c1, const Real3D c2, 
                const Real3D box) const
        {
            double x, y, z;
            x = minimal_distance(c1[0], c2[0], box[0]);
            y = minimal_distance(c1[1], c2[1], box[1]);
            z = minimal_distance(c1[2], c2[2], box[2]);

            //return sqrt(x*x + y*y + z*z);
            return Real3D(x, y, z);
        }


    public:
        /* Constructors */
        Frame(){}

        Frame(Real3D box) : _box(box) {}

        Frame(std::vector<Real3D> coordinates, std::vector<Real3D> velocities,
                Real3D box, size_t particles_per_molecule = 1)
            : _coordinates(coordinates), _velocities(velocities),
            _box(box), _particles_per_molecule(particles_per_molecule) {}

        Frame(std::ifstream & stream, const size_t particles_per_molecule = 0,
                const char format = 'x')
        {
            if(!stream) 
            {
                throw std::runtime_error(std::string("file input error"));
            }
            switch(format) {
                case 'x'  :
                    read_xyz(stream, particles_per_molecule);
                    break;
                default :
                    throw std::runtime_error("Unknown file type.");
            }
        }

        Frame(const char * filename, const size_t particles_per_molecule = 0,
                const char format = 'x') 
            : _particles_per_molecule(particles_per_molecule)
        {
            std::ifstream stream(filename);
            if(!stream) 
            {
                throw std::runtime_error(std::string("file ") + filename 
                        + std::string(" could not be found"));
            }
            switch(format) {
                case 'x'  :
                    read_xyz(stream, particles_per_molecule);
                    break;
                default :
                    throw std::runtime_error("Unknown file type.");
            }
        }

        // operators
        friend std::ostream & operator<<(std::ostream & os, const Frame & frame)
        {
            if(frame.is_null())
            {
                os << "Null frame\n";
            }
            else
            {
                os << "Frame with " << frame.size() << " particles in "
                    << frame.number_of_molecules() << " molecules "
                    << "of length " << frame.particles_per_molecule() <<
                    ", Box: "
                    << frame.box() << '\n';
            }
            return os;
        }

        /**
         * This evaluates true if lhs and rhs are the same by value
         */
        friend bool operator==(const Frame& lhs, const Frame& rhs)
        {
            return lhs._coordinates == rhs._coordinates
                && lhs._velocities == rhs._velocities
                && lhs._box == rhs._box
                && lhs._particles_per_molecule == rhs._particles_per_molecule;
        }

        // member functions
        void read_xyz(const std::string & filename,
                const size_t particles_per_molecule = 0)
        {
            std::ifstream stream(filename);
            read_xyz(stream, particles_per_molecule);
        }

        void read_xyz(std::ifstream & stream,
                const size_t particles_per_molecule = 0)
        {
            std::string str;
            std::getline(stream, str);
            if(stream.eof())
            {
                return;
            }

            // get number of particles
            size_t number_of_particles = std::stoi(str);
            _coordinates.clear();
            _coordinates.reserve(number_of_particles);

            // iterate box dimensions
            std::getline(stream, str);
            std::string buf;            // Have a buffer string
            std::stringstream ss0(str);  // Insert the string into a stream
            ss0 >> buf;
            _box[0] = stod(buf);
            ss0 >> buf;
            _box[1] = stod(buf);
            ss0 >> buf;
            _box[2] = stod(buf);

            // iterate coordinates and eventually velocities
            Real3D vec;
            for(size_t i_line = 0; i_line < number_of_particles; ++i_line)
            {
                // TODO reset instead?
                std::getline(stream, str);
                std::stringstream ss(str);

                // read types, coordinates and eventually velocities into vec
                ss >> buf;
                _types.push_back(stoi(buf));
                ss >> buf;
                vec[0] = stod(buf);
                ss >> buf;
                vec[1] = stod(buf);
                ss >> buf;
                vec[2] = stod(buf);
                _coordinates.push_back(vec);

                // also add velocities if there are any
                if(ss >> buf)
                {
                    _velocities.reserve(number_of_particles);
                    vec[0] = stod(buf);
                    ss >> buf;
                    vec[1] = stod(buf);
                    ss >> buf;
                    vec[2] = stod(buf);
                    _velocities.push_back(vec);
                }
            }

            // sanity check
            assert( _coordinates.size() == number_of_particles);
            assert( _types.size() == number_of_particles);
            assert( _velocities.size() == 0
                    || _velocities.size() == number_of_particles );

            _velocities.shrink_to_fit();
            _particles_per_molecule = particles_per_molecule;
        }

        // getter
        bool has_velocities() const {return _velocities.size() > 0;};

        bool is_null() const {return _coordinates.size() == 0;};

        size_t size() const {return _coordinates.size();}

        Real3D box() const {return _box;}

        double box(size_t i) const
        {
            assert(i < 3);
            return _box[i];
        }

        bool is_square() const 
        {
            return box(0) == box(1)
                && box(1) == box(2)
                && box(2) == box(0);
        }

        size_t particles_per_molecule() const
        {
            if(_particles_per_molecule == 0)
                return size();
            else
                return _particles_per_molecule;
        }

        size_t number_of_molecules() const
        {
            if(size() == 0)
                return 0;
            else
                return size() / particles_per_molecule();
        }

        auto c_cbegin() const {return _coordinates.cbegin();}
        auto c_cend() const {return _coordinates.cend();}
        auto v_cbegin() const {return _velocities.cbegin();}
        auto v_cend() const {return _velocities.cend();}
        auto t_cbegin() const {return _types.cbegin();}
        auto t_cend() const {return _types.cend();}

        Real3D coordinate(const size_t index) const
        {
            assert(index < size());
            return _coordinates[index];
        }

        Real3D folded_coordinate(const size_t index) const
        {
            return fold(coordinate(index));
        }

        Real3D velocity(const size_t index) const
        {
            assert(index < size());
            return _velocities[index];
        }

        int type(const size_t index) const
        {
            assert(index < size());
            if(_types.empty()) return (int) index / particles_per_molecule();
            else return _types[index];
        }

        Molecule molecule(const size_t index) const
        {
            if(index >= number_of_molecules())
            {
                std::cerr << *this;
                throw std::out_of_range("Molecule index out of range.");
            }
            size_t ppm = particles_per_molecule();
            size_t start_index = index * ppm;

            auto c_begin = _coordinates.begin() + start_index;
            auto c_end = c_begin + ppm;

            if(!has_velocities())
            {
                return Molecule(c_begin, c_end);
            }
            else
            {
                auto v_begin = _velocities.begin() + start_index;
                auto v_end = v_begin + ppm;
                return Molecule(c_begin, c_end, v_begin, v_end);
            }
        }

        // manipulation
        void set_box(const Real3D box){_box = box;}

        void set_particles_per_molecule(const size_t particles_per_molecule)
        {
            _particles_per_molecule = particles_per_molecule;
        }

        void set_number_of_molecules(const size_t number_of_molecules)
        {
            assert(number_of_molecules >= 1);
            if(size() % number_of_molecules != 0)
            {
                std::cerr << *this;
                throw std::runtime_error("Number of molecules is not "
                        "commensurable with number of particles.");
            }
            set_particles_per_molecule(size() / number_of_molecules);
        }

        void set_types(size_t max_type)
        {
            for(size_t i_part = 0; i_part < size(); ++i_part)
            {
                _types[i_part] = (i_part / particles_per_molecule()) %
                    max_type;
            }
        }

        void set_types()
        {
            for(size_t i_part = 0; i_part < size(); ++i_part)
            {
                _types[i_part] = i_part / particles_per_molecule();
            }
        }

        void reduce_types(size_t max_type)
        {
            for(int type : _types)
            {
                type %= max_type;
            }
        }

        void add_particle(Real3D coordinate, int type = 0)
        {
            _coordinates.push_back(coordinate);
            _types.push_back(type);
            if(has_velocities()) _velocities.push_back(Real3D(0.));
        }

        void add_particle(Real3D coordinate, Real3D velocity, int type = 0)
        {
            _coordinates.push_back(coordinate);
            _velocities.push_back(velocity);
            _types.push_back(type);
        }

        void remove_particle(const size_t index)
        {
            _coordinates.erase(_coordinates.begin() + index);
            if(has_velocities())
                _velocities.erase(_velocities.begin() + index);
        }

        void clear()
        {
            _box = Real3D(0.);
            _particles_per_molecule = 0;
            _coordinates.clear();
            if(has_velocities())
                _velocities.clear();
        }

        void fold()
        {
            for(std::vector<Real3D>::iterator cit = _coordinates.begin();
                    cit != _coordinates.end(); ++cit)
            {
                *cit = fold(*cit);
            }
            return;
        }

        void fold_2d()
        {
            for(std::vector<Real3D>::iterator cit = _coordinates.begin();
                    cit != _coordinates.end(); ++cit)
            {
                *cit = fold_2d(*cit);
            }
            return;
        }

        void shift_coordinates(const Real3D shift)
        {
            for(std::vector<Real3D>::iterator cit = _coordinates.begin();
                    cit != _coordinates.end(); ++cit)
            {
                (*cit) += shift;
            }
            return;
        }

        void scale_box(const double factor)
        {
            for(std::vector<Real3D>::iterator cit = _coordinates.begin();
                    cit != _coordinates.end(); ++cit)
            {
                (*cit) *= factor;
            }
            _box *= factor;
            return;
        }

        void crop_box(const double xmin, const double xmax,
                const double ymin, const double ymax,
                const double zmin, const double zmax)
        {
            // crop box down to specified size. This will break connectivity.
            assert(xmin > xmax && ymin > ymax && zmin > zmax);
            Real3D c;
            std::vector<Real3D>::iterator it = _coordinates.begin();
            while(it != _coordinates.end())
            {
                c = fold(*it);

                if(c[0] < xmin && c[0] >= xmax) _coordinates.erase(it);
                else if(c[1] < ymin && c[1] >= ymax) _coordinates.erase(it);
                else if(c[2] < zmin && c[2] >= zmax) _coordinates.erase(it);
                else ++it;
            }

            _box = Real3D(xmax - xmin, ymax - ymin, zmax - zmin);
            return;
        }

        void crop_box(const double lx, const double ly,
                const double lz)
        {
            crop_box(0., lx, 0., ly, 0., lz);
            return;
        }

        void rotate_x(const double angle)
        {
            double c = cos(angle);
            double s = sin(angle);

            for(auto it = _coordinates.begin(); it != _coordinates.end(); ++it)
            {
                auto crd = *it;
                *it = Real3D{
                    crd * Real3D{ 1,  0,  0},
                    crd * Real3D{ 0,  c, -s},
                    crd * Real3D{ 0,  s,  c}
                };
            }
        }

        void rotate_y(const double angle)
        {
            double c = cos(angle);
            double s = sin(angle);

            for(auto it = _coordinates.begin(); it != _coordinates.end(); ++it)
            {
                auto crd = *it;
                *it = Real3D{
                    crd * Real3D{ c,  0,  s},
                    crd * Real3D{ 0,  1,  0},
                    crd * Real3D{-s,  0,  c}
                };
            }
        }

        void rotate_z(const double angle)
        {
            double c = cos(angle);
            double s = sin(angle);

            for(auto it = _coordinates.begin(); it != _coordinates.end(); ++it)
            {
                auto crd = *it;
                *it = Real3D{
                    crd * Real3D{ c, -s,  0},
                    crd * Real3D{ s,  c,  0},
                    crd * Real3D{ 0,  0,  1}
                };
            }
        }

        void rotate(const Real3D vector, const double angle)
        {
            const double x = vector[0];
            const double y = vector[1];
            const double z = vector[2];

            const double c = cos(angle);
            const double omc = 1. - c;
            const double s = sin(angle);

            for(auto it = _coordinates.begin(); it != _coordinates.end(); ++it)
            {
                auto crd = *it;
                *it = Real3D{
                    crd * Real3D{x*x*omc+c,   x*y*omc-z*s, x*z*omc+y*s},
                    crd * Real3D{x*y*omc+z*s, y*y*omc+c,   y*z*omc-x*s},
                    crd * Real3D{x*z*omc-y*s, y*z*omc+x*s, z*z*omc+c}
                };
            }
        }

        std::vector<Frame> divide_box(const int nx, const int ny = 1,
                const int nz = 1) const
        {
            // divide box into nx x ny x nz equal-sized subboxes
            Real3D new_box(box(0)/nx, box(1)/ny, box(2)/nz);

            const int nf = nx * ny * nz;
            const double dx = box(0) / nx;
            const double dy = box(1) / ny;
            const double dz = box(2) / nz;

            int ix, iy, iz, index;

            std::vector<Frame> frames(nf);
            for(int i = 0; i < nf; ++i)
            {
                frames[i].set_box(new_box);
            }

            Real3D c;
            for(auto it = c_begin(); it != c_end(); ++it)
            {
                c = fold(*it);

                ix = (int) (c[0] / dx);
                iy = (int) (c[1] / dy);
                iz = (int) (c[2] / dz);

                index = iz + nz*(iy + ny*ix);

                frames[index].add_particle(fold(c));
            }

            return frames;
        }

        /**
         * Takes a slice of thickness thickness along the plane given by an 
         * initial point x0 and a normal vector norm
         * @param[in] thickness thickness of the slice
         * @param[in] x0 initial point
         * @param[in] norm normal vector
         * @param[out] sliced Frame
         */
        Frame slice(const double thickness, 
                const Real3D x0 = Real3D(0.), 
                const Real3D norm = Real3D(-1., -1., 2.)) const
        {
            double dist;
            Real3D n = norm / norm.abs();
            Real3D coordinate;
            Frame frame{box()};

            for(size_t i_part = 0; i_part < size(); ++i_part)
            {
                coordinate = folded_coordinate(i_part);
                dist = (coordinate - x0) * n;

                if(fabs(dist) <= thickness / 2.)
                {
                    frame.add_particle(coordinate, velocity(i_part), type(i_part));
                }

            }
            return frame;
        }


        Frame slice_square(const double thickness, const double height = 0.) const
        {
            return slice(thickness, Real3D(0., 0., height), 
                    Real3D(0., 0., 1.));
        }



        Frame slice_rectangle(const double thickness) const
        {
            return slice(thickness, Real3D(0., 0., 0.), 
                    Real3D(-1., 0., 2.));
        }


        /**
         * Multiplies configuration. Works only for unfolded
         */
        Frame multiply(const size_t mx=2, const size_t my=2,
                const size_t mz=2) const
        {
            Frame new_frame;
            const Real3D b = box();
            new_frame.set_box(2 * b);

            for(size_t dx = 0; dx < mx; ++dx)
            {
                for(size_t dy = 0; dy < my; ++dy)
                {
                    for(size_t dz = 0; dz < mz; ++dz)
                    {
                        auto cit = c_begin();
                        auto vit = v_begin();
                        while(cit != c_end())
                        {
                            if(has_velocities())
                            {
                                new_frame.add_particle(
                                        *cit + Real3D(dx*b[0],dy*b[1],dz*b[2]),
                                        *vit);
                            }
                            else
                            {
                                new_frame.add_particle(
                                        *cit + Real3D(dx*b[0],dy*b[1],dz*b[2]));
                            }
                            ++cit;
                            ++vit;
                        }
                    }
                }
            }
            return new_frame;
        }

        /* analysis */
        double mean_distance() const
        {
            Real3D d;
            long double dist = 0.;

            for(auto i = c_begin(); i != c_end(); ++i)
            {
                for(auto j = c_begin(); j != i; ++j)
                {
                    dist += (minimal_distance(*j, *i, box())).abs();
                }
            }
            dist *= pow(size(), -2.);
            return dist;
        }

        double smallest_distance() const
        {
            double smallest = box(0);
            double dist;

            for(auto i = c_begin(); i != c_end(); ++i)
            {
                for(auto j = c_begin(); j != i; ++j)
                {
                    dist = minimal_distance(*j, *i, _box).abs();
                    if (dist < smallest)
                    {
                        smallest = dist;
                    }
                }
            }
            return smallest;
        }

        /**
         * Calculates mean value of function f for all molecules that are
         * contained in the system.
         * @param[in] f function that is to applied to the molecules.
         */
        template <typename T>
            T mean(T (Molecule::*f)(void)) const
            {
                T m = 0.;
                for(size_t im = 0; im < number_of_molecules(); ++im)
                {
                    m += (molecule(im).*f)();
                }
                return m / number_of_molecules();
            }

        /**
         * Returns the value of function f for each molecule in the system in
         * form of a vector.
         * @param[in] f function that is to be applied to the molecules.
         */
        template <typename T>
            std::vector<T> vector(T (Molecule::*f)(void)) const
            {
                std::vector<T> vec;
                vec.reserve(number_of_molecules());
                for(size_t im = 0; im < number_of_molecules(); ++im)
                {
                    vec.push_back((molecule(im).*f)());
                }
                assert(vec.size() == number_of_molecules());
                return vec;
            }

        /**
         * Return the size of the largest bond in the system.
         */
        double max_bond_length() const
        {
            double largest = 0.;
            for(size_t im = 0; im < number_of_molecules(); ++im)
            {
                if(molecule(im).max_bond_length() > largest)
                    largest = molecule(im).max_bond_length();
            }
            return largest;
        }

        /**
         * Calculate atomar mean-squared-displacement with respect to reference
         * frame.
         * @param[in] earlier_frame reference frame
         */
        double mean_squared_displacement(Frame earlier_frame)
        {
            double msd = 0.;
            for(size_t imol = 0; imol < number_of_molecules(); ++imol)
            {
                msd += molecule(imol).
                    mean_squared_displacement(earlier_frame.molecule(imol));
            }
            return msd / number_of_molecules();
        }

        /**
         * Write configuration to .xyz file.
         * @param[in] filename filename
         * @param[in] append whether or not to append to existing files
         * @param[in] Indicate if set to 0, every particle gets identity zero
         * in the first column of the file. For nonzero values, the molecule
         * index modulo the value of indicate is used as identity. This is
         * because VMD only supports a very limited range of integers as
         * identity.
         */
        void write_xyz(const char* filename, const bool append = false, const
                size_t precision = 11) const
        {
            if(is_null()) return;

            std::ofstream stream;

            if(append) stream.open(filename, std::ios_base::out | std::ios_base::app);
            else stream.open(filename);

            // write number of particles
            stream << std::fixed;
            stream << size() << '\n';

            // write box
            stream.precision(precision);
            stream << std::scientific;
            stream << box() << '\n';

            size_t number_of_digits
                = (size_t) log10(number_of_molecules()) + 1;

            for(size_t i_part = 0; i_part < size(); ++i_part)
            {
                // write type
                stream << std::fixed;
                stream << std::setfill(' ');
                stream << std::setw(number_of_digits);
                stream << type(i_part);

                // write coordinate
                stream << std::scientific;
                for(size_t i_dim = 0; i_dim < 3; ++i_dim)
                {
                    stream << std::setfill(' ');
                    stream << std::setw(precision + 8);
                    stream << coordinate(i_part)[i_dim];
                }

                // write velocity if any
                if(has_velocities())
                {
                    for(size_t i_dim = 0; i_dim < 3; ++i_dim)
                    {
                        stream << std::setfill(' ');
                        stream << std::setw(precision + 8);
                        stream << velocity(i_part)[i_dim];
                    }
                }
                stream << '\n';
            }
            stream.close();
            return;
        }

        void write_vtk(const char* filename, const bool append = false, const
                size_t precision = 11) const
        {
            if(is_null()) return;

            std::ofstream stream;
            if(append)
            {
                stream.open(
                        filename, std::ios_base::out | std::ios_base::app);
            }
            else stream.open(filename);

            stream.precision(precision);
            stream << std::scientific;

            // write header
            stream << "# vtk DataFile Version 1.0\n";
            stream << "Polymer System\n";
            stream << "ASCII\n\n";
            stream << "DATASET POLYDATA\n";
            stream << "POINTS " << size() << " float\n";

            // write particle coordinates
            for(const Real3D coordinate : _coordinates)
                stream << coordinate << '\n';

            // write bonds
            stream << "\nLINES " << number_of_molecules() << ' '
                << number_of_molecules() + size() << '\n';

            size_t i_part;
            for(size_t i_mol = 0; i_mol < number_of_molecules(); ++i_mol)
            {
                stream << number_of_molecules();

                for(size_t i_mon = 0; i_mon < particles_per_molecule(); ++i_mon)
                {
                    i_part = i_mol * particles_per_molecule() + i_mon;
                    stream << ' ' << i_part;
                }
                stream << '\n';
            }

            stream.close();
            return;
        }

        void write_pdb(const char* filename, const bool append = false) const
        {

            if(is_null()) return;

            std::ofstream stream;
            if(append) stream.open(filename, std::ios_base::out |
                    std::ios_base::app);
            else stream.open(filename);


            stream << "CRYST1";

            // write box
            stream << std::fixed;
            stream.precision(3);
            for(size_t i_dim = 0; i_dim < 3; ++i_dim)
            {
                stream << std::setfill(' ');
                stream << std::setw(9);
                stream << box(i_dim);
            }

            // write angles
            stream.precision(2);
            for(size_t i_dim = 0; i_dim < 3; ++i_dim)
            {
                stream << std::setfill(' ');
                stream << std::setw(7);
                stream << 90.;
            }
            stream << " P 1           1\n";

            if( (int) log10(number_of_molecules()) > 3)
            {
                throw std::runtime_error("too many chains");
            }

            for(size_t i_part = 0; i_part < size(); ++i_part)
            {
                stream << "ATOM  ";

                stream << std::fixed;
                stream << std::setfill(' ');
                stream << std::setw(5);
                // atom serial number
                stream << 0 << ' ';
                // atom name, alternate location indicator, residue name, 
                // chain identifier
                stream << "FE   UNX F";

                // type
                stream << std::setfill(' ');
                stream << std::setw(4);
                stream << type(i_part) << "    ";


                // coordinates
                stream << std::defaultfloat;
                stream << std::fixed;
                stream.precision(3);
                for(size_t i_dim = 0; i_dim < 3; ++i_dim)
                {
                    stream << std::setfill(' ');
                    stream << std::setw(8);
                    stream << coordinate(i_part)[i_dim];
                }

                stream.precision(2);

                stream << std::fixed;
                stream << std::setfill(' ');
                stream << std::setw(6);
                stream << "0.0";

                stream << std::fixed;
                stream << std::setfill(' ');
                stream << std::setw(6);
                stream << "0.0";

                stream << '\n';
            }
            stream << "END" << '\n';
            stream.close();
            return;
        }

        void write_lattice(const char* filename, const size_t lattice_size) const
        {
            // write lattice configuration to disk
            double * lattice = new double[lattice_size*lattice_size*lattice_size]{0.};
            read_lattice(*this, lattice, lattice_size);

            couf::write_3d_array_to_file(lattice, filename, lattice_size); 

            delete [] lattice;
        }

        // TODO
        /**
         * Write configuration to binary file
         */
        void write_binary(const char* filename, const bool append = false) const
        {
            if(is_null()) return;

            std::ofstream stream;

            if(append) stream.open(filename, std::ios_base::out
                    | std::ios_base::app | std::ios::binary);
            else stream.open(filename, std::ios::out | std::ios::binary);

            // write velocity flag
            bool velocity_flag = has_velocities();
            stream.write(reinterpret_cast<char *>(&velocity_flag),
                    sizeof(velocity_flag));

            // write number of particles
            size_t number_of_particles = size();
            stream.write(reinterpret_cast<const char *>(&number_of_particles),
                    sizeof(number_of_particles));

            // write box
            stream.write(const_cast<char *>(
                        reinterpret_cast<const char *>(&_box[0])),
                    3*sizeof(double));

            // write coordinates
            stream.write(const_cast<char *>(
                        reinterpret_cast<const char *>(&_coordinates[0])),
                    _coordinates.size() * sizeof(_coordinates[0]));

            // write velocities
            if(has_velocities())
            {
                stream.write(const_cast<char *>(
                            reinterpret_cast<const char *>(&_velocities[0])),
                        _velocities.size() * sizeof(_velocities[0]));
            }
            stream.close();
            return;
        }

        void make_sphere(const double radius)
        {
            if(2*radius > box(0) || 2*radius > box(1) || 2*radius > box(2))
            {
                throw std::runtime_error("radius is to big for box");
            }

            Real3D center = box() / 2.;
            Real3D r;

            for(size_t i_x = 0; i_x < box(0); ++i_x)
            {
                for(size_t i_y = 0; i_y < box(1); ++i_y)
                {
                    for(size_t i_z = 0; i_z < box(2); ++i_z)
                    {
                        r = Real3D((double) i_x, (double) i_y, (double) i_z);
                        if((r - center).abs() < radius)
                        {
                            add_particle(r);
                        }
                    }
                }
            }
        }

        void make_cube(const double L)
        {
            if(2*L > box(0) || 2*L > box(1) || 2*L > box(2))
            {
                throw std::runtime_error("L is to big for box");
            }

            Real3D start = (box() - Real3D(L)) / 2.;

            for(size_t i_x = 0; i_x < L; ++i_x)
            {
                for(size_t i_y = 0; i_y < L; ++i_y)
                {
                    for(size_t i_z = 0; i_z < L; ++i_z)
                    {
                        add_particle(start + Real3D(i_x, i_y, i_z));
                    }
                }
            }
        }

        friend double read_lattice(const Frame & frame, double *lattice, const
                size_t lattice_size, Real3D *velocity_lattice=nullptr);
        friend void read_lattice_2d(const Frame & frame, double *lattice, const
                int lattice_size);
        friend double read_lattice(const char * filename, double *lattice,
                const size_t lattice_size);
        friend double read_lattice_2d(const char * filename, double *lattice,
                const size_t lattice_size);
};

/**
 * Projects the 3d configuration onto a lattice using 2nd order
 * extrapolation. Box must be square. Repeated application adds
 * the new configuration on top of the old.
 * @param[in] frame input configuration
 * @param[out] lattice array that stores the lattice
 * @param[in] lattice_size linear size of the lattice
 * @param[out] velocity_lattice array that stores the velocity lattice
 */
double read_lattice(const Frame & frame, double *lattice, const size_t
        lattice_size, Real3D *velocity_lattice)
{
    const int original_lattice_size = round(frame.box(0));
    if(original_lattice_size == 0.) throw
        std::runtime_error("You have to set a nonzero box");

    const double lattice_constant = (double) lattice_size / original_lattice_size;
    std::vector<Real3D>::const_iterator vcit;

    int i, j, k, new_i, new_j, new_k;
    long int index;
    double weight, x, y, z, dx, dy, dz;
    double norm = 0;

    if(velocity_lattice != nullptr)
    {
        vcit = frame.v_cbegin();
    }
    for(auto ccit = frame.c_cbegin(); ccit != frame.c_cend(); ++ccit)
    {
        Real3D coordinate = *ccit;

        // coordinates in lattice units
        x = coordinate[0] * lattice_constant;
        y = coordinate[1] * lattice_constant;
        z = coordinate[2] * lattice_constant;

        x = couf::fold(x, lattice_size);
        y = couf::fold(y, lattice_size);
        z = couf::fold(z, lattice_size);

        // indices of the lattice site to the bottom left:
        i = (int) x;
        j = (int) y;
        k = (int) z;

        // loop over surrounding lattice sites and compute weights
        // according to distance
        for(int di = 0; di < 2; ++di)
        {
            dx = fabs(i + di - x);

            for(int dj = 0; dj < 2; ++dj)
            {
                dy = fabs(j + dj - y);

                for(int dk = 0; dk < 2; ++dk)
                {
                    dz = fabs(k + dk - z);

                    // at the borders, again periodic boundary conditions have
                    // to be taken into account
                    new_i = (i + di) % lattice_size;
                    new_j = (j + dj) % lattice_size;
                    new_k = (k + dk) % lattice_size;

                    index = new_k + lattice_size
                        * (new_j + lattice_size * new_i);

                    assert(index >= 0);
                    assert(index < pow(lattice_size, 3));

                    weight = (1. - dx)*(1. - dy)*(1. - dz);
                    norm += weight;

                    assert(weight >= 0. && weight <= 1.);

                    lattice[index] += weight;

                    if(velocity_lattice != nullptr)
                    {
                        velocity_lattice[index]
                            += weight*(*vcit);
                    }
                }
            }
        }
        if(velocity_lattice != nullptr)
        {
            ++vcit;
        }
    }
    assert(norm - frame.size() < 10e-8);

    return norm;
}

void read_lattice_2d(const Frame & frame, double *lattice, const int lattice_size)
{
    const int original_lattice_size = round(frame.box(0));
    const double lattice_constant = (double) lattice_size / original_lattice_size;

    int i, j, new_i, new_j = 0;
    double x, y, dx, dy = 0.;

    for(auto ccit = frame.c_cbegin(); ccit != frame.c_cend(); ++ccit)
    {
        Real3D coordinate = *ccit;
        /* coordinates in lattice units */
        x = coordinate[0] * lattice_constant;
        y = coordinate[1] * lattice_constant;

        x = couf::fold(x, lattice_size);
        y = couf::fold(y, lattice_size);

        /* indices of the lattice site to the bottom left: */
        i = (int) x;
        j = (int) y;

        /* loop over surrounding lattice sites and compute weights
         * according to distance */
        for(int di = 0; di < 2; ++di)
        {
            dx = fabs(i + di - x);

            for(int dj = 0; dj < 2; ++dj)
            {
                dy = fabs(j + dj - y);

                /* at the borders, again periodic boundary conditions have
                 * to be taken into account */
                new_i = (i + di) % lattice_size;
                new_j = (j + dj) % lattice_size;

                lattice[new_j + lattice_size * new_i] += (1 - dx)*(1 - dy);
            }
        }
    }
    return;
}

/**
 * Read 3d lattice from file. Lattice must be cubic.
 * @param[in] filename input file name
 * @param[out] lattice array that stores the lattice
 * @param[in] lattice_size linear size of the lattice
 */
double read_lattice(const char * filename, double *lattice, const size_t lattice_size)
{
    /* read lattice from datafile with format:
     * i_x i_y i_z density */
    std::ifstream instream(filename);
    if (!instream)
    {
        throw std::runtime_error("Could not open input file\n");
    }

    /* read the lattice from file */
    int i, j, k;
    int index = 0;
    double norm = 0.;
    double val;
    std::string str;
    std::string buf;
    while(getline(instream, str))
    {
        std::stringstream ss(str);
        ss >> buf;
        i = stoi(buf);
        ss >> buf;
        j = stoi(buf);
        ss >> buf;
        k = stoi(buf);
        ss >> buf;
        val = stof(buf);
        norm += val;
        ss >> buf;

        index = k + lattice_size * j + pow(lattice_size, 2) * i;
        lattice[index] = val;
    }
    instream.close();

    if(norm <= 0.)
        std::cerr << "WARNING: System seems to be empty or has negative density\n";

    if(index != pow(lattice_size, 3) - 1)
        throw std::runtime_error("Wrong lattice size given\n");

    return norm;
}

/**
 * Read 2d lattice from file. Lattice must be square.
 * @param[in] filename input file name
 * @param[out] lattice array that stores the lattice
 * @param[in] lattice_size linear size of the lattice
 */
double read_lattice_2d(const char * filename, double *lattice, const size_t
        lattice_size)
{
    /* read lattice from datafile */
    std::ifstream instream(filename);
    if (!instream)
    {
        throw std::runtime_error("Could not open input file\n");
    }

    /* read the lattice from file */
    int index = 0;
    double norm = 0.;
    double val;
    std::string str;
    std::string buf;
    while(getline(instream, str))
    {
        std::stringstream ss(str);
        while(ss >> buf)
        {
            val = stof(buf);
            lattice[index++] = val;
            norm += val;
        }
    }
    instream.close();

    if(norm <= 0.)
        std::cerr << "WARNING: System seems to be empty or has negative density\n";

    printf("Calculated a total density of %16.9e.\n", norm);

    if(index != pow(lattice_size, 2) - 1)
        throw std::runtime_error("Wrong lattice size given\n");

    return norm;
}

/**
 * Transform radial symmetric image of a Fourier transform S(\vec{q}) to a
 * vector of pairs (q, S(q)). The wave vector q has a resolution of
 * bin_width.
 * @param[in] lattice_transformed FFTW fourier transformation of a density lattice
 * @param[in] lattice_size linear size of the lattice
 * @param[in] lattice_constant lattice constant
 * @param[in] bin_width bin width
 * @return vector of pairs (q, S(q)).
 */
std::vector<std::vector<double>> linearize_lattice(const fftw_complex
        *lattice_transformed, const int lattice_size, const double
        lattice_constant, const double bin_width = 0.05, const double norm = 1.)
{
    std::vector<std::vector<double> > result;
    result.reserve(1000); // TODO estimate

    int index, inew, jnew, knew;
    double q = 0., S = 0., re, im, abs = 0., dist;
    const double q_0 = 2. * M_PI / (lattice_size * lattice_constant);

    //const double max = 10.;
    const int threshold = (lattice_size + 1) / 2;
    const double max = q_0 * sqrt(3.) * (threshold - 1) + 1.e-10;
    const int n_bins = ceil(max / bin_width);
    int bin;

    double *histogram = (double*) calloc(n_bins, sizeof(double));
    double *abscissa  = (double*) calloc(n_bins, sizeof(double));
    int *count        = (int*)    calloc(n_bins, sizeof(int));

    /* build the histogram */
    for(int i = 0; i < lattice_size; ++i)
    {
        for(int j = 0; j < lattice_size; ++j)
        {
            for(int k = 0; k < lattice_size / 2 + 1; ++k)
            {
                /*  skip trivial mode */
                if(i == 0 && j ==0 && k == 0) continue;

                index = k + (lattice_size / 2 + 1)
                    * (j + lattice_size  * i);

                /* we look for the closest distance to a corner */
                if(i < threshold) inew = i;
                else inew = (lattice_size - 1) - i;

                if(j < threshold) jnew = j;
                else jnew = (lattice_size - 1) - j;

                if(k < threshold) knew = k;
                else knew = (lattice_size - 1) - k;

                dist = sqrt(inew*inew + jnew*jnew + knew*knew);

                re = lattice_transformed[index][0];
                im = lattice_transformed[index][1];
                abs = re*re + im*im;

                q = q_0 * dist;
                S = abs / norm;

                assert(q <= max);
                bin = (int) (q / bin_width);
                histogram[bin] += S;
                abscissa[bin]  += q;
                ++count[bin];
            }
        }
    }

    /* build the vector */
    std::vector<double> tuple(2); // TODO not efficient
    for(int i = 0; i < n_bins; ++i)
    {
        if(count[i] > 0){
            // tuple[0] = abscissa[i] / count[i]; // TODO average abscissa
            tuple[0] = bin_width * (i + 0.5);
            tuple[1] = histogram[i] / count[i];
        }
        else
        {
            tuple[0] = bin_width * (i + 0.5);
            tuple[1] = 0.;
        }
        result.push_back(tuple);
    }
    free(histogram);
    free(abscissa);
    free(count);

    return result;
}

/**
 * Calculate the structure factor of a frame.
 * The configuration is first mapped to a density map on a lattice. Then a
 * Fourier transform of the lattice is performed. This gives the structure
 * factor depending on VECTOR q. The lattice is then linearized to obtain
 * the dependence on the absolute value of q.
 * @param[in] input Input lattice configuration. This can either be an array or
 * a filename.
 * @param[in] bin_width bin width
 * @param[in] lattice_size linear lattice size
 * @param[in] lattice_constant lattice constant
 * @param[in] bin_width bin width
 * @param[in] norm normalization constant
 * @return vector of pairs (q, S(q)).
 */
template<typename T>
std::vector<std::vector<double>> structure_factor(const T input, const size_t
        lattice_size, const double lattice_constant, const double bin_width =
        0.1, const double norm = 1.)
{
    if(typeid(T) != typeid(Frame) && typeid(T) != typeid(const char *))
        throw std::runtime_error("Incompatible input type.\n");

    const size_t number_of_sites = pow(lattice_size, 3);
    const size_t reduced_number_of_sites = pow(lattice_size, 2) *
        (lattice_size / 2 + 1);

    /* allocate space for lattices */
    std::unique_ptr<fftw_complex> lattice_transformed{
        fftw_alloc_complex(reduced_number_of_sites * sizeof(fftw_complex))
    };
    for(size_t i = 0; i < reduced_number_of_sites; ++i)
    {
        lattice_transformed.get()[i][0] = 0.;
        lattice_transformed.get()[i][1] = 0.;
    }
    //memset(lattice_transformed.get(), 0.,
    //2*reduced_number_of_sites*sizeof(double)); // TODO

    std::unique_ptr<double>
        lattice{fftw_alloc_real(number_of_sites * sizeof(double))};
    std::fill(lattice.get(), lattice.get() + number_of_sites, 0.);

    /* generate fftw plan */
    fftw_plan plan = fftw_plan_dft_r2c_3d(lattice_size, lattice_size,
            lattice_size, lattice.get(), lattice_transformed.get(),
            FFTW_ESTIMATE); // TODO inplace?

    /* read lattice with function that is suitable for type of input */
    read_lattice(input, lattice.get(), lattice_size);

    /* do the transformation */
    fftw_execute(plan);

    return linearize_lattice(lattice_transformed.get(), lattice_size,
            lattice_constant, bin_width, norm);
}

/**
 * Calculate the exact structure factor for given wave vector q
 * @param[in] frame input configuration
 * @param[in] q wave vector
 * @return S(q)
 */
double structure_factor(const Frame & frame, const Real3D q)
{
    double re = 0.;
    double im = 0.;
    double sp;

    for(auto it = frame.c_cbegin(); it != frame.c_cend(); ++it)
    {
        sp = q * (*it);

        re += cos(sp);
        im += sin(sp);
    }

    double abs_squared = re*re + im*im;
    return abs_squared / frame.size();
}

/**
 * Compute the set of lattice vectors that lie within a spherical shell of
 * certain thickness.
 * @param[in] radius radius of shell
 * @param[in] thickness thickness of shell
 * @param[in] lattice_constant lattice constant
 * @return list of vectors
 */
std::vector<Real3D> lattice_vectors_inside_shell(const double radius,
        const double thickness, const double lattice_constant)
{
    std::vector<Real3D> list;
    Real3D vec;

    /* work in lattice units */
    const double r_lower = std::max(0., (radius - .5 * thickness) / lattice_constant);
    const double r_upper = (radius + .5 * thickness) / lattice_constant;
    const size_t n_items = std::ceil(
            4.*M_PI/3.*(pow(r_upper + 0.5*sqrt(3.), 3) - pow(r_lower - 0.5*sqrt(3.), 3)));
    list.reserve(n_items);

    for(int i = -std::floor(r_upper); i <= std::ceil(r_upper); ++i)
    {
        for(int j = -std::floor(r_upper); j <= std::ceil(r_upper); ++j)
        {
            for(int k = -std::floor(r_upper); k <= std::ceil(r_upper); ++k)
            {
                if(i*i + j*j + k*k >= r_lower * r_lower
                        && i*i + j*j + k*k <= r_upper * r_upper)
                {
                    vec[0] = i;
                    vec[1] = j;
                    vec[2] = k;
                    list.push_back(vec * lattice_constant);
                }
            }
        }
    }

    return list;
}

/**
 * Method for the brute-force calculation of the structure factor S(q).
 * All lattice vectors in a shell of thickness lattice_constant around q are
 * considered if there are less than n_rand. If there are more, n_rand vectors
 * are chosen at random.
 * @param[in] frame
 * @param[in] q absolute value of vave vector to consider
 * @param[in] n_rand number of random orientations to consider
 * @return vector of tuples (<q>, <S(q)>, \sigma(q))
 */
std::vector<double>  mean_structure_factor(const Frame & frame, const double q,
        const size_t n_rand = 128)
{
    std::vector<double> q_vec;
    q_vec.reserve(n_rand);
    std::vector<double> S_vec;
    S_vec.reserve(n_rand);

    const double lattice_constant = 2.*M_PI / std::min(frame.box(0),
            std::min(frame.box(1), frame.box(2)));
    auto list = lattice_vectors_inside_shell(q, lattice_constant, lattice_constant);
    if(list.size() == 0)
    {
        std::cerr << "WARNING: shell empty\n";
    }
    if(list.size() <= n_rand)
    {
        for(auto it = list.cbegin(); it != list.cend(); ++it)
        {
            S_vec.push_back(structure_factor(frame, *it));
            q_vec.push_back(it->abs());
        }
    }
    else
    {
        // pick items randomly
        std::random_device rd;
        std::mt19937 rng(rd());
        std::shuffle(list.begin(), list.end(), rng);

        for(size_t i_rand = 0; i_rand < n_rand; ++i_rand)
        {
            S_vec.push_back(structure_factor(frame, list[i_rand]));
            q_vec.push_back(list[i_rand].abs());
        }
    }
    std::vector<double> result(3);

    double mean_S = std::accumulate(S_vec.begin(), S_vec.end(), .0)/S_vec.size();
    double sq_sum_S = std::inner_product(S_vec.begin(), S_vec.end(), S_vec.begin(), 0.0);
    double stdev_S = std::sqrt(sq_sum_S / S_vec.size() - mean_S * mean_S);

    result[0] = q;
    result[1] = mean_S;
    result[2] = stdev_S;
    return result;
}


/**
 * Calculate Minkowski functionals (MFs). In 3 dimensions there are 4:
 * V_0: volume
 * V_1: area
 * V_2: mean curvature
 * V_3: Euler-Poincare characteristic
 * The configuration is first mapped onto a black and white lattice. Each
 * cell-center of the lattice has 8 neighbors. Hence there are 2^8 different
 * neighborhoods. The MFs are additive and rotationally invariant. By symmetry
 * the number of neighborhoods that are unique wrt. the MFs reduces to 22.
 * The MFs are computed in the lattice centers.
 *
 * See paper Arns, Knackstedt, Pinczewski, Mecke Phys. Rev. E 63 2001
 *
 * @param[in] input Input lattice configuration either as Frame or
 * filename const char *
 * @param[in] lattice_size linear lattice size of the interpolation lattice
 * @param[in] threshold lattice sites with density >= threshold will be
 * interpreted as 'black'.
 */
template<typename T>
std::array<double, 6> minkowski_functionals(const T input, 
        const size_t lattice_size, const double threshold = -1., 
        const char norm='n', const bool natural_units=false)
{
    if(typeid(T) != typeid(Frame) && typeid(T) != typeid(const char *))
        throw std::runtime_error("Incompatible input type.\n");

    std::array<double, 6> result{0., 0., 0., 0., 0., 0.};

    // values of the MFs for the 22 possible configurations
    static const std::array<std::array<double, 6>, 22> results{{
        {0., 0., 0., 0., 0., 0.},
            {1., 3., 3., 3., 1., 1.},
            {2., 4., 2., 2., 0., 0.},
            {2., 6., 6., 2., 2., -2.},
            {2., 6., 6., 6., 2., -6.},
            {3., 5., 1., 1., -1., -1.},
            {3., 7., 5., 1., 1., -3.},
            {3., 9., 9., -3., 3., -1.},
            {4., 4., 0., 0., 0., 0.},
            {4., 6., 0., 0., -2., -2.},
            {4., 6., 0., 0., -2., -2.},
            {4., 8., 4., -4., 0., 0.},
            {4., 8., 4., -4., 0., 0.},
            {4., 12., 12., -12., 4., 4.},
            {5., 9., 3., -9., -1., 3.},
            {5., 7., -1., -5., -3., 1.},
            {5., 5., -1., -1., -1., -1.},
            {6., 6., -6., -6., -6., 2.},
            {6., 6., -2., -6., -2., 2.},
            {6., 4., -2., -2., 0., 0.},
            {7., 3., -3., -3., 1., 1.},
            {8., 0., 0., 0., 0., 0.}
    }};

    // associate the 256 possible neighborhoods to the 22 configurations
    static const std::array<const std::array<double, 6> *, 256> pointers{
        &results[0],
            &results[1],
            &results[1],
            &results[2],
            &results[1],
            &results[2],
            &results[3],
            &results[5],
            &results[1],
            &results[3],
            &results[2],
            &results[5],
            &results[2],
            &results[5],
            &results[5],
            &results[8],
            &results[1],
            &results[2],
            &results[3],
            &results[5],
            &results[3],
            &results[5],
            &results[7],
            &results[9],
            &results[4],
            &results[6],
            &results[6],
            &results[10],
            &results[6],
            &results[10],
            &results[11],
            &results[16],
            &results[1],
            &results[3],
            &results[2],
            &results[5],
            &results[4],
            &results[6],
            &results[6],
            &results[10],
            &results[3],
            &results[7],
            &results[5],
            &results[9],
            &results[6],
            &results[11],
            &results[10],
            &results[16],
            &results[2],
            &results[5],
            &results[5],
            &results[8],
            &results[6],
            &results[10],
            &results[11],
            &results[16],
            &results[6],
            &results[11],
            &results[10],
            &results[16],
            &results[12],
            &results[15],
            &results[15],
            &results[19],
            &results[1],
            &results[3],
            &results[4],
            &results[6],
            &results[2],
            &results[5],
            &results[6],
            &results[10],
            &results[3],
            &results[7],
            &results[6],
            &results[11],
            &results[5],
            &results[9],
            &results[10],
            &results[16],
            &results[2],
            &results[5],
            &results[6],
            &results[10],
            &results[5],
            &results[8],
            &results[11],
            &results[16],
            &results[6],
            &results[11],
            &results[12],
            &results[15],
            &results[10],
            &results[16],
            &results[15],
            &results[19],
            &results[3],
            &results[7],
            &results[6],
            &results[11],
            &results[6],
            &results[11],
            &results[12],
            &results[15],
            &results[7],
            &results[13],
            &results[11],
            &results[14],
            &results[11],
            &results[14],
            &results[15],
            &results[18],
            &results[5],
            &results[9],
            &results[10],
            &results[16],
            &results[10],
            &results[16],
            &results[15],
            &results[19],
            &results[11],
            &results[14],
            &results[15],
            &results[18],
            &results[15],
            &results[18],
            &results[17],
            &results[20],
            &results[1],
            &results[4],
            &results[3],
            &results[6],
            &results[3],
            &results[6],
            &results[7],
            &results[11],
            &results[2],
            &results[6],
            &results[5],
            &results[10],
            &results[5],
            &results[10],
            &results[9],
            &results[16],
            &results[3],
            &results[6],
            &results[7],
            &results[11],
            &results[7],
            &results[11],
            &results[13],
            &results[14],
            &results[6],
            &results[12],
            &results[11],
            &results[15],
            &results[11],
            &results[15],
            &results[14],
            &results[18],
            &results[2],
            &results[6],
            &results[5],
            &results[10],
            &results[6],
            &results[12],
            &results[11],
            &results[15],
            &results[5],
            &results[11],
            &results[8],
            &results[16],
            &results[10],
            &results[15],
            &results[16],
            &results[19],
            &results[5],
            &results[10],
            &results[9],
            &results[16],
            &results[11],
            &results[15],
            &results[14],
            &results[18],
            &results[10],
            &results[15],
            &results[16],
            &results[19],
            &results[15],
            &results[17],
            &results[18],
            &results[20],
            &results[2],
            &results[6],
            &results[6],
            &results[12],
            &results[5],
            &results[10],
            &results[11],
            &results[15],
            &results[5],
            &results[11],
            &results[10],
            &results[15],
            &results[8],
            &results[16],
            &results[16],
            &results[19],
            &results[5],
            &results[10],
            &results[11],
            &results[15],
            &results[9],
            &results[16],
            &results[14],
            &results[18],
            &results[10],
            &results[15],
            &results[15],
            &results[17],
            &results[16],
            &results[19],
            &results[18],
            &results[20],
            &results[5],
            &results[11],
            &results[10],
            &results[15],
            &results[10],
            &results[15],
            &results[15],
            &results[17],
            &results[9],
            &results[14],
            &results[16],
            &results[18],
            &results[16],
            &results[18],
            &results[19],
            &results[20],
            &results[8],
            &results[16],
            &results[16],
            &results[19],
            &results[16],
            &results[19],
            &results[18],
            &results[20],
            &results[16],
            &results[18],
            &results[19],
            &results[20],
            &results[19],
            &results[20],
            &results[20],
            &results[21]
    };
    assert(pointers.size() == 256);

    const size_t number_of_sites = pow(lattice_size, 3);
    std::unique_ptr<double[]> lattice = std::make_unique<double[]>(number_of_sites);
    for(size_t i = 0; i < number_of_sites; ++i)  lattice[i] = 0.;

    /* read lattice with function that is suitable for type of input */
    const double mean_density = read_lattice(input, lattice.get(), lattice_size)
        / number_of_sites;

    double new_threshold{0.};
    if(threshold < 0.) new_threshold = mean_density;
    else new_threshold = threshold;

    // loop over lattice centers
    size_t xn, yn, zn, i_neigh, linear_neigh, config, tag, n_black{0};
    for(size_t i_x = 0; i_x < lattice_size; ++i_x)
    {
        for(size_t i_y = 0; i_y < lattice_size; ++i_y)
        {
            for(size_t i_z = 0; i_z < lattice_size; ++i_z)
            {
                // loop over 8 neighbors
                config = 0;
                i_neigh = 0;
                for(size_t d_x = 0; d_x < 2; ++d_x)
                {
                    for(size_t d_y = 0; d_y < 2; ++d_y)
                    {
                        for(size_t d_z = 0; d_z < 2; ++d_z)
                        {
                            // fold
                            xn = (i_x + d_x) % lattice_size;
                            yn = (i_y + d_y) % lattice_size;
                            zn = (i_z + d_z) % lattice_size;

                            linear_neigh
                                = zn + lattice_size * (yn + lattice_size * xn);

                            assert(linear_neigh < number_of_sites);

                            tag = (lattice[linear_neigh] >= new_threshold);
                            config += pow(2, i_neigh) * tag;
                            assert(config < 256);
                            n_black += tag;

                            ++i_neigh;
                        } // d_z
                    } // d_y
                } // d_x

                result += *pointers[config];

            } // i_z
        } // i_y
    } // i_x

    assert(result[0] == n_black);

    switch(norm) {
        case 'n':
            // none
            result[0] /= 8.;
            result[1] /= 24.;
            result[2] /= 24.;
            result[3] /= 24.;
            result[4] /= 8.;
            result[5] /= 8.;
            break;
        case 's':
            // spherical
            result[0] /= 8.;
            result[1] /= 6.;
            result[2] *= (M_PI/12.);
            result[3] *= (M_PI/12.);
            result[4] /= 8.;
            result[5] /= 8.;
            break;
        case 'c':
            // cubic
            result[0] /= 8.;
            result[1] /= 4.;
            result[2] *= (M_PI/8.);
            result[3] *= (M_PI/8.);
            result[4] /= 8.;
            result[5] /= 8.;
            break;
        default :
            throw std::runtime_error("invalid norm argument");
    }

    if(natural_units)
    {
        result[0] /= pow(lattice_size, 3);
        result[1] /= pow(lattice_size, 2);
        result[2] /= lattice_size;
        result[3] /= lattice_size;
    }


    return result;
}
#endif
