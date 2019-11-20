#ifndef _TRAJECTORY_HPP
#define _TRAJECTORY_HPP

#include "frame.hpp"
#include <list>

/*!
 * Reader class to read trajectory files. 
 * Supported formats:
 * .xyz
 * Class can be used like a forward iterator. Dereferencing returns the current
 * Frame object.
 */
class Trajectory
{
    private:
        std::ifstream _stream; //!< stream of the input file
        size_t _index{0}; //!< current frame index
        size_t _particles_per_molecule{0}; //!< number of particles in one molecule
        Frame _frame; //!< current configuration
        double _timestep{0.}; //!< timestep between consecutive frames
        mutable size_t _frames_read{0};
        std::string _last_string;

        /**
         * Advances the trajectory by one frame without reading the new frame 
         * into memory
         */
        size_t advance_one()
        {
            if(!_stream || _stream.peek() == EOF) return 0;
            // basic function that advances _stream by one frame
            std::string str;
            std::getline(_stream, str);

            // this is to catch possible double newlines at the end of the file
            if(str == "")
            {
                do {
                    std::getline(_stream, str);
                } while(str == "" && _stream.good());
                return 0;
            }

            // get number of particles
            size_t number_of_particles = std::stoi(str);

            // iterate box and coordinates
            for(size_t i_line = 0; i_line < 1 + number_of_particles; ++i_line)
            {
                _stream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                //_stream.ignore(8*1024, '\n');
            }
            return 1;
        }

    public:
        Trajectory(const Trajectory&) = delete; // disable copy constructor

        Trajectory(const std::string & filename,
                size_t particles_per_molecule = 0, size_t offset = 0):
            _stream(std::ifstream(filename)),
            _index(0),
            _particles_per_molecule(particles_per_molecule),
            _frame(Frame(_stream, particles_per_molecule)),
            _timestep(0.),
            _frames_read{1}
        {
            advance(offset);
        }

        // getter
        //bool is_null() const {return _stream.eof();}
        bool is_null() const {
            return !is_good();
        }

        bool is_good() const {
            return _stream.good();
        }

        size_t index() const {return _index;}

        size_t frames_read() const {return _frames_read;}

        const Frame & frame() const 
        {
            return _frame;
        }

        double timestep() const {return _timestep;}

        size_t particles_per_molecule() const
        {
            if(_particles_per_molecule == 0)
                return frame().size();
            else
                return _particles_per_molecule;
        }

        size_t number_of_molecules() const
        {
            if(frame().size() == 0)
                return 0;
            else
                return frame().size() / particles_per_molecule();
        }

        size_t size()
        {
            const size_t start_index = index();
            size_t index = start_index;
            while(!is_null())
            {
                advance_one();
                index++;
            }
            move_to(start_index);
            return index;
        }

        // manipulation
        void reset()
        {
            _stream.clear();                 // clear fail and eof bits
            _stream.seekg(0, std::ios::beg); // back to the start!
            _index = 0;
            _frame = Frame(_stream, _particles_per_molecule);
            _frames_read = 1;
        }

        void advance(const size_t number=1)
        {
            if(number == 0 || _stream.peek() == EOF) return;
            for(size_t i = 0; i < number - 1; ++i)
            {
                _index += advance_one();
            }
            _index += 1;
            if(!is_null())
            {
                _frame = Frame(_stream, _particles_per_molecule);
                ++_frames_read;
            }
        }

        void loop_advance(int argc, char **argv)
        {
            size_t offset = static_cast<size_t>(
                    atoi(couf::parse_arguments(argc, argv, "--offset")));
            int max = atoi(couf::parse_arguments(argc, argv, "--max", "-1"));
            size_t step = static_cast<size_t>(
                    atoi(couf::parse_arguments(argc, argv, "--step", "1")));
            double exponent = atof(couf::parse_arguments(argc, argv, "--exp"));

            size_t ndx = index();
            size_t m;
            if(max < 0) m = std::numeric_limits<int>::max();
            else m = static_cast<size_t>(max);

            if(ndx == 0 && offset) 
            {
                advance(offset);
                return;
            }

            // advance exponentially
            if(exponent)
            {
                size_t ds = 0;
                if(ndx == 0) advance();
                else 
                {
                    ds = ceil(ndx * (pow(2, exponent) - 1));
                    if(index() + ds > m)
                    {
                        _stream.setstate(std::ios_base::eofbit);
                        return;
                    }
                    advance(ds);   
                }
                //else advance(ceil(ndx * (pow(10, exponent) - 1)));
            }
            // advance linearly
            else advance(step);

            if(ndx >= m)
            {
                _stream.setstate(std::ios_base::eofbit);
                return;
            }

            return;
        }

        void move_to(const size_t number)
        {
            if(number >= index())
            {
                advance(number - index());
            }
            else
            {
                _stream.clear();                 // clear fail and eof bits
                _stream.seekg(0, std::ios::beg); // back to the start!
                _index = 0;
                advance(number + 1);
                _index -= 1;
            }
        }

        // operators
        Trajectory & operator=(const Trajectory&) = delete; // disable assign

        Frame operator*() const {return frame();} // dereference

        //const Frame * operator->() const {return &_frame;}
        Frame * operator->() {return &_frame;}

        Frame operator[](size_t index)
        {
            size_t start_index = this->index();
            move_to(index);
            Frame frame = this->frame();
            move_to(start_index);
            return frame;
        }

        Trajectory & operator++() // prefix
        {
            advance();
            return *this;
        }

        Trajectory & operator+=(const size_t number)
        {
            advance(number);
            return *this;
        }

        friend std::ostream & operator<<(std::ostream & os,
                const Trajectory & trajectory)
        {
            if(trajectory.is_null())
            {
                os << "Trajectory ended at index: "
                    << trajectory.index() << '\n'
                    << "Frames read: " << trajectory.frames_read() << '\n';
            }
            else
            {
                os << trajectory.index() << ": " << *trajectory
                    << "Frames read " << trajectory.frames_read() << '\n';
            }
            return os;
        }

        // analysis
        template <typename T>
            double mean(T (Molecule::*f)(void), size_t step = 1)
            {
                double m = 0.;
                while(!is_null())
                {
                    m += frame().mean(f);
                    advance(step);
                }
                return m / index();
            }

        template <typename T>
            Timeseries<T> timeseries(T (Frame::*f)(void), size_t step = 1,
                    size_t max = std::numeric_limits<int>::max())
            {
                std::vector<T> data;
                while(!is_null() && index() <= max)
                {
                    data.emplace_back( (frame().*f)() );
                    advance(step);
                }
                reset();
                return Timeseries<T>(data);
            }

        /**
         * @param[in] f function that produces observable of interest
         * @param[in] step only process every step'th frame
         * @param[in] max maximum frame to process
         * @return timeseries where each timestep contains a vector of
         * observables f (for each molecule)
         */
        template <typename T>
            Timeseries<std::vector<T>> timeseries(T (Molecule::*f)(void), 
                    size_t step = 1,
                    size_t max = std::numeric_limits<int>::max())
            {
                std::vector<std::vector<T>> data;
                while(!is_null() && index() <= max)
                {
                    data.emplace_back( frame().vector(f) );
                    advance(step);
                }
                reset();
                return Timeseries<std::vector<T>>(data);
            }

        /**
         * @param[in] step only process every step'th frame
         * @param[in] max maximum frame to process
         * @return set oftimeseries where each timestep contains a vector of
         * observables f (for each molecule)
         */
        std::vector<Timeseries<std::vector<Real3D>>> 
            timeseries_set(
                    std::vector<Real3D (Molecule::*)()> vofp =
                    { &Molecule::end_to_end, &Molecule::center_of_mass },
                    size_t step = 1, 
                    size_t max = std::numeric_limits<int>::max())
            {
                std::vector<Timeseries<std::vector<Real3D>>> data(vofp.size());
                while(!is_null() && index() <= max)
                {
                    Frame f = frame();
                    std::vector<std::vector<Real3D>> vec(vofp.size());
                    for(size_t im = 0; im < f.number_of_molecules(); ++im)
                    {
                        Molecule m = f.molecule(im);

                        for(size_t i_f = 0; i_f < vofp.size(); ++i_f)
                        {
                            vec[i_f].push_back((m.*vofp[i_f])());
                        }
                    }
                    for(size_t i = 0; i < data.size(); ++i)
                    {
                        data[i].push_back(vec[i]);
                    }
                    advance(step);
                }
                reset();

                return data;
            }

        /**
         * @param[in] f function that produces observable of interest
         * @param[in] step only process every step'th frame
         * @param[in] max maximum frame to process
         * @return timeseries where each timestep contains the value of f
         * averaged over all molecules
         */
        template <typename T>
            Timeseries<T> timeseries_mean(T (Molecule::*f)(void), 
                    size_t step = 1,
                    size_t max = std::numeric_limits<int>::max())
            {
                std::vector<T> data;
                while(!is_null() && index() <= max)
                {
                    data.emplace_back( frame().mean(f) );
                    advance(step);
                }
                reset();
                return Timeseries<T>(data);
            }

        template <typename T>
            Timeseries<T> timeseries_single_molecule
            (T (Molecule::*f)(void), size_t i_mol = 0, size_t step = 1,
             size_t max = std::numeric_limits<int>::max())
            {
                // make timeseries where each timestep contains a observable f
                // for molecule i_mol
                std::vector<T> data;
                while(!is_null() && index() <= max)
                {
                    data.emplace_back( (frame().molecule(i_mol).*f)() );
                    advance(step);
                }
                reset();
                return Timeseries<T>(data);
            }

        void write_xyz(const char* filename)
        {
            while(!is_null())
            {
                frame().write_xyz(filename, true);
                advance();
            }
            return;
        }
};

    template <typename T>
double correlation_function(T (Molecule::*f)(void), 
        Trajectory & trajectory, const size_t span)
{
    T val_1, val_2;
    double mean = 0.;
    size_t count = 0;

    const size_t n_mol = trajectory->number_of_molecules();
    std::list<Frame> cache;

    // fill frame cache
    while(trajectory.index() < span)
    {
        cache.emplace_back(Frame(*trajectory));
        ++trajectory;
    }
    while(!trajectory.is_null())
    {
        // TODO make molecule iterator
        Frame frame_1, frame_2;
        frame_2 = *trajectory;
        if(span == 0) frame_1 = frame_2;
        else
        {
            frame_1 = *(cache.begin());
            cache.pop_front(); 
            cache.emplace_back(frame_2);
        }

        for(size_t i_mol = 0; i_mol < n_mol; ++i_mol)
        {
            val_1 = (frame_1.molecule(i_mol).*f)();
            val_2 = (frame_2.molecule(i_mol).*f)();
            if(span == 0) assert(val_1 == val_2);

            mean += val_2 * val_1;
            ++count;
        }
        ++trajectory;
    }
    mean /= count;
    trajectory.reset();

    return mean;
}
#endif
