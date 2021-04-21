#ifndef _MCTIMESERIES_HPP
#define _MCTIMESERIES_HPP

#include "couf.hpp"
#include "vector.hpp"
#include "timeseries.hpp"
#include <algorithm>
#include <sstream>
#include <iomanip>

/* @file mctimeseries.hpp 
 * Container for multi column timeseries.
 */

template<typename T>
class MCTimeseries
{
    private:
        std::vector<std::vector<T>> _data;
        double _timestep{0.};
        double _initial_time{0.};
        std::string _comment{""};

    public:
        MCTimeseries() {}

        MCTimeseries(std::vector<std::vector<T>> data) : _data(data) {}

        MCTimeseries(const char * filename, size_t offset = 0)
        {
            read_all(filename, offset);
        }

        void push_back(std::vector<T> data)
        {
            _data.push_back(data);
        }

        void emplace_back(std::vector<T>&& data)
        {
            _data.emplace_back(data);
        }

        void read_all(const char * filename, size_t offset = 0)
        {
            std::ifstream stream(filename);

            if (!stream) throw std::runtime_error("Could not open file");
            if(_data.size() != 0) throw std::runtime_error("_data not empty");

            // get file dimensions
            size_t n_col = couf::count_columns(filename);
            size_t n_lin = couf::count_lines(filename);
            // reserve space
            for(size_t i_c = 0; i_c < n_col; ++i_c)
            {
                _data.emplace_back(std::vector<T>());
                _data[i_c].reserve(n_lin);
            }
            assert(n_col == _data.size());

            // do the reading
            std::string str;
            for(size_t i = 0; i < offset; ++i)
            {
                while(str[0] == '#' || str[0] == '\0')
                    std::getline(stream, str); // skip comments
                std::getline(stream, str);
            }
            do {
                while(str[0] == '#' || str[0] == '\0')
                    std::getline(stream, str); // skip comments
                std::string buf;            // Have a buffer string
                std::stringstream ss(str);  // Insert the string into a stream

                // read value
                for(size_t i_col = 0; i_col < n_col; ++i_col)
                {
                    ss >> buf;
                    _data[i_col].push_back(stod(buf));
                }

            } while (std::getline(stream, str));

            stream.close();
            assert(consistent());
        }

        void read(const char * filename, size_t column, size_t offset = 0)
        {
            std::ifstream stream(filename);
            if (!stream) throw std::runtime_error("Could not open file");
            read(stream, column, offset);
            stream.close();
        }

        void read(std::ifstream & stream, size_t column, size_t offset = 0)
        {
            double previous = 0;
            double current = 0;
            std::string str;
            std::getline(stream, str);
            _data.emplace_back(std::vector<T>());

            for(size_t i = 0; i < offset; ++i)
            {
                while(str[0] == '#' || str[0] == '\0')
                    std::getline(stream, str); // skip comments
                std::getline(stream, str);
            }
            do {
                while(str[0] == '#' || str[0] == '\0')
                    std::getline(stream, str); // skip comments
               std::string buf;            // Have a buffer string
                std::stringstream ss(str);  // Insert the string into a stream

                // read value
                for(size_t j = 0; j < column; ++j)
                {
                    assert(ss);
                    ss >> buf;
                }
                _data[_data.size()-1].push_back(stod(buf));

            } while (std::getline(stream, str));

            if(timestep() < 0)
                throw std::runtime_error("Something is wrong with the timestep");

            assert(consistent());
            return;
        }

        bool consistent() const
        {
            for(auto d: _data)
            {
                if(d.size() != _data[0].size())
                {
                    return false;
                }
            }
            return true;
        }

        size_t size() const {return _data.size();}

        size_t number_of_columns() const {return _data.size();}

        size_t number_of_steps() const {return _data[0].size();}

        std::vector<T> data(size_t i_col) const
        {
            return _data[i_col];
        }

        double timestep() const 
        {
            assert(_data.size() > 0);
            assert(_data[0].size() > 1);

            if(_timestep == 0) return _data[0][1] - _data[0][0];
            else return _timestep;
        }

        double initial_time() const
        {
            assert(_data.size() > 0);
            assert(_data[0].size() > 0);
            return _data[0][0];
        }

        std::string comment() const 
        {
            return _comment;
        }

        void set_timestep(double timestep) {_timestep = timestep;}

        void set_initial_time(double initial_time) {_initial_time = initial_time;}

        void set_comment(std::string comment)
        {
            _comment = comment;
        }

        double time(const int step) const
        {
            return _initial_time + step * _timestep;
        }

        Timeseries<T> timeseries(size_t i_col)
        {
            return Timeseries<T> (data(i_col), timestep(), initial_time(),
                comment());
        }
};
#endif
