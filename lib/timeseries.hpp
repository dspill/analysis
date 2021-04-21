#ifndef _TIMESERIES_HPP
#define _TIMESERIES_HPP

#include "couf.hpp"
#include "vector.hpp"
#include <algorithm>
#include <sstream>
#include <iomanip>

/** @file timeseries.hpp */

/** The Timeseries class stores a time series of data points and contains some
 * functions for their analysis.
 * A data point can be a vector of points as well as the type is templated.
 * A time step, initial time and a comment can be given.
 */
template<typename T>
class Timeseries
{
    private:
        std::vector<T> _data;
        double _timestep{0.};
        double _initial_time{0.};
        std::string _comment{""};

    public:
        Timeseries() {}

        Timeseries(std::vector<T> data, 
                double timestep = 0, 
                double initial_time = 0,
                std::string comment = "") : _data(data), _timestep(timestep),
            _initial_time(initial_time), _comment(comment) {}

        Timeseries(const char * filename, size_t column = 1, size_t offset = 0)
        {
            read(filename, column, offset);
        }

        void read(const char * filename, size_t column = 1, size_t offset = 0)
        {
            std::ifstream stream(filename);
            if (!stream) throw std::runtime_error("Could not open file");
            read(stream, column, offset);
            stream.close();
        }

        void read_all(const char * filename, size_t offset = 0)
        {
            std::ifstream stream(filename);
            if (!stream) throw std::runtime_error("Could not open file");
            read_all(stream, offset);
            stream.close();
        }

        /* Read Timeseries from infile stream given a column. 
         * Lines starting with '#' are skipped.
         * @param[in] column Column of interest.
         * @param[in] offset Optional number of lines to skip.
         */
        void read(std::ifstream & stream, size_t column = 1, size_t offset = 0)
        {
            double previous = 0;
            double current = 0;
            size_t i = 0;
            std::string str;
            std::getline(stream, str);

            while(i < offset)
            {
                while(str[0] == '#') std::getline(stream, str); // skip comments
                std::getline(stream, str);
                ++i;
            }
            do {
                while(str[0] == '#') std::getline(stream, str); // skip comments
                std::string buf;            // Have a buffer string
                std::stringstream ss(str);  // Insert the string into a stream

                // get timestep TODO skip if column == 0
                ss >> buf;
                if(column == 0)
                {
                    _timestep = 1.;
                }
                else
                {
                    current = stod(buf);
                    if(_initial_time == 0.) _initial_time = current;
                    _timestep = current - previous;
                    previous = current;
                }

                // read value
                for(size_t j = 0; j < column; ++j)
                {
                    ss >> buf;
                }
                _data.push_back(stod(buf));

            } while (std::getline(stream, str));

            if(timestep() < 0)
                throw std::runtime_error("Something is wrong with the timestep");
            return;
        }

        /* Read Timeseries from all columns of an infile stream. 
         * For each time step a vector of doubles is stored.
         * Lines starting with '#' are skipped.
         * @param[in] column Column of interest.
         * @param[in] offset Optional number of lines to skip.
         */
        void read_all(std::ifstream & stream, size_t offset = 0)
        {
            std::string str;
            std::getline(stream, str);
            std::vector<double> line;

            for(size_t i = 0; i < offset; ++i)
            {
                while(str[0] == '#') std::getline(stream, str); // skip comments
                std::getline(stream, str);
            }

            do {
                while(str[0] == '#') std::getline(stream, str); // skip comments

                std::string buf;            // Have a buffer string
                std::stringstream ss(str);  // Insert the string into a stream

                // read value
                while(ss >> buf) line.push_back(stod(buf)); 
                _data.push_back(line);
                line.clear();
            } while (std::getline(stream, str));

            return;
        }

        /** @return The total number of data points. */
        size_t size() const {return _data.size();}

        /** @return The timestep that is set. */
        double timestep() const {return _timestep;}

        /** @return The comment that is set. */
        std::string comment() const {return _comment;}

        /** @return The initial time that is set. */
        double initial_time() const {return _initial_time;}

        /** @return Append new value to the end of the Timeseries. */
        void push_back(T val)
        {
            _data.push_back(val);
        }

        void emplace_back(T&& val)
        {
            _data.emplace_back(val);
        }

        /** Set the timestep. */
        void set_timestep(double timestep) {_timestep = timestep;}

        void set_comment(std::string comment)
        {
            _comment = comment;
        }

        /** Set the initial time. */
        void set_initial_time(double initial_time) {_initial_time = initial_time;}

        double time(const int step) const
        {
            return _initial_time + step * _timestep;
        }

        /** @return Constant reference to first data piont in form of an
         * iterator.  */
        typename std::vector<T>::const_iterator begin() const
        {
            return _data.begin();
        }

        /** @return Constant reference to last data piont in form of an
         * iterator.  */
        typename std::vector<T>::const_iterator end() const
        {
            return _data.end();
        }

        // analysis
        /* @return Average value of data points. */
        T mean() const
        {
            auto it = begin();
            T mean = *it;
            ++it;
            while(it != end())
            {
                mean += *it;
                ++it;
            }

            return mean / size();
        }

        /* @return Standard deviation of data points. */
        T stdev() const
        {
            double mean = this->mean();

            double sq_sum = std::inner_product(begin(), end(), begin(), 0.0);
            double stdev = std::sqrt(sq_sum / (size() - 1) - mean * mean);

            return stdev;
        }

        /* Calculates the average autocorrelation function between data points
         * that are a certain number of steps apart.
         * @param[in] Number of steps between data points.
         * @return Autocorrelation function. */
        double autocorrelation_function(size_t span) const
        {
            if(span > size())
                throw std::out_of_range("Span exceeds size of dataset.");

            size_t count = 0;
            double mean = 0.;
            T mean_i{0.};
            T mean_f{0.};

            auto it_i = begin();
            auto it_f = begin() + span;

            while(it_f != end())
            {
                mean += *it_i * *it_f;
                mean_i += *it_i;
                mean_f += *it_f;

                ++it_i;
                ++it_f;
                ++count;
            }
            return (mean - mean_i * mean_f / count) / count;
        }

        // operators
        /** Print info. */
        friend std::ostream & operator<<(std::ostream & os,
                const Timeseries & timeseries)
        {
            os << "Timeseries of lenght " << timeseries.size() << ".\n";
            return os;
        }

        /** Add a second timeseries to the current element wise. */
        Timeseries & operator+=(const Timeseries &ts)
        {
            assert(this->size() == ts.size());
            this->_data = this->_data + ts._data;
            return *this;
        }

        /** Add two Timeseries element wise. */
        Timeseries operator+(const Timeseries &ts) const
        {
            Timeseries ts2 = *this;
            ts2 += ts;
            return ts2;
        }

        /** Substract two Timeseries element wise. */
        Timeseries & operator-=(const Timeseries &ts)
        {
            assert(this->size() == ts.size());
            this->_data = this->_data - ts._data;
            return *this;
        }

        /** Substract a second Timeseries from the current element wise. */
        Timeseries operator-(const Timeseries &ts) const
        {
            Timeseries ts2 = *this;
            ts2 -= ts;
            return ts2;
        }

        /** Multiply current Timeseries element wise by a scalar. */
        Timeseries & operator*=(const double a)
        {
            this->_data *= a;
            return *this;
        }
        
        /** Multiply two Timeseries element wise by a scalar. */
        Timeseries operator*(const double a)
        {
            Timeseries ts2 = *this;
            ts2 *= a;
            return ts2;
        }

        /** Divide current Timeseries element wise by a scalar. */
        Timeseries & operator/=(const double a)
        {
            this->_data /= a;
            return *this;
        }

        /** Divide two Timeseries element wise by a scalar. */
        Timeseries operator/(const double a)
        {
            Timeseries ts2 = *this;
            ts2 /= a;
            return ts2;
        }

        // equality
        /** Test whether the data points of two Timeseries are equal by value.
         */
        bool operator==(const Timeseries &ts) const {
            assert(this->size() == ts.size());
            return std::equal(this->begin(), this->end(), ts.begin());
        }

        /** Test whether the data points of two Timeseries are not equal by
         * value.  */
        bool operator!=(const Timeseries &ts) const {
            return !(*this == ts);
        }

        void coarsen(const size_t factor) // TODO does not work
        {
            if(factor == 1) return;
            T mean;
            std::vector<T> new_data;
            new_data.reserve(size()/factor);
            auto it = begin();
            while(it != end())
            {
                if(end() - it < static_cast<int>(factor)) break;
                mean = *it;
                ++it;
                for(size_t i = 1; i < factor; ++i)
                {
                    mean += *it;
                    ++it;
                }
                mean /= static_cast<double>(factor);
                new_data.push_back(mean);
            }
            _data = new_data;
            _initial_time += (factor - 1) * timestep() / 2.;
            _timestep *= factor;
            return;
        }


        /** Write Timeseries to file. */
        void write(const char * filename) const
        {
            constexpr size_t precision = 11;
            std::ofstream stream(filename);

            if(!comment().empty())
                stream << "# " << comment() << '\n';

            size_t step = 0;
            stream << std::scientific << std::setprecision(precision);
            for(auto it = begin(); it != end(); ++it)
            {
                if(_timestep > 0.) stream << time(step);
                //for (const auto &val : *it)  TODO fix below
                //for (auto val = it->begin(); val != it->end(); ++val) 
                //stream << std::setw(precision + 7) << val;	
                //
                stream << std::setw(precision + 7) << *it;	
                stream << '\n';

                ++step;

            }
            return;
        }

        /** Clear all data from Timeseries. */
        void clear()
        {
            _data.clear();
            _timestep = 0.;
            _initial_time = 0.;
            _comment = "";
        }

        T operator [](size_t i) const {return _data[i];}
        T const& operator [](size_t i) {return _data[i];}
};
#endif
