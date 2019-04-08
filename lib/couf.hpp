#ifndef _COUF_HPP
#define _COUF_HPP
/* Collection of Useful Functions */
#include <cmath>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstring>
#include <cfloat>
#include <vector>
#include "Real3D.hpp"

namespace couf
{
    int my_assert(bool condition)
    {
        if(condition) return 0;
        else exit(1);
    }

    /* folding */
    double fold(const double coordinate, const double box_length)
    {
        if(box_length == 0.) return 0.;

        double folded_coordinate;

        folded_coordinate = fmod(coordinate, box_length);
        folded_coordinate += (coordinate < 0.) * box_length;

        return folded_coordinate;
    }

    Real3D fold(const Real3D coordinate, const Real3D box)
    {
        double x, y, z;
        x = fold(coordinate[0], box[0]);
        y = fold(coordinate[1], box[1]);
        z = fold(coordinate[2], box[2]);

        return Real3D(x, y, z);
    }

    Real3D fold_2d(const Real3D coordinate, const Real3D box)
    {
        double x, y, z;
        x = fold(coordinate[0], box[0]);
        y = fold(coordinate[1], box[1]);
        z = coordinate[2];

        return Real3D(x, y, z);
    }

    double minimal_distance(const double c1, const double c2,
            const double box_length)
    {
        /* compute the minimal distance between two particles along one
         * dimension while considering the periodic images of the box */
        double fc1 = fold(c1, box_length);
        double fc2 = fold(c2, box_length);

        double dist = fc2 - fc1;
        dist = dist - box_length * round(dist/box_length);

        return dist;
    }

    Real3D minimal_distance(const Real3D c1, const Real3D c2, const Real3D box)
    {
        double x, y, z;
        x = minimal_distance(c1[0], c2[0], box[0]);
        y = minimal_distance(c1[1], c2[1], box[1]);
        z = minimal_distance(c1[2], c2[2], box[2]);

        //return sqrt(x*x + y*y + z*z);
        return Real3D(x, y, z);
    }

    /* write 3d table to file */
    void write_3d_array_to_file(const double *array, const char *filename,
            const size_t n_x, size_t n_y = 0, size_t n_z = 0)
    {
        if(n_y == 0) n_y = n_x;
        if(n_z == 0) n_z = n_x;

        FILE *pFile;
        pFile = fopen(filename, "w");

        size_t index = 0;

        for(size_t i_x = 0; i_x < n_x; ++i_x)
        {
            for(size_t i_y = 0; i_y < n_y; ++i_y)
            {
                for(size_t i_z = 0; i_z < n_z; ++i_z)
                {
                    fprintf(pFile, "%4ld %4ld %4ld %16.9e\n", 
                            i_x, i_y, i_z, array[index++]);
                }
            }
        }
        fclose(pFile);
    }

    void write_3d_Real3D_array_to_file(const Real3D *array, 
            const char *filename, const size_t n_x, size_t n_y = 0, size_t n_z = 0)
    {
        if(n_y == 0) n_y = n_x;
        if(n_z == 0) n_z = n_x;

        FILE *pFile;
        pFile = fopen(filename, "w");

        size_t index = 0;

        for(size_t i_x = 0; i_x < n_x; ++i_x)
        {
            for(size_t i_y = 0; i_y < n_y; ++i_y)
            {
                for(size_t i_z = 0; i_z < n_z; ++i_z)
                {
                    fprintf(pFile, "%4ld %4ld %4ld %16.9e %16.9e %16.9e\n", 
                            i_x, i_y, i_z, 
                            array[index][0], array[index][1], array[index][2]);
                    ++index;
                }
            }
        }
        fclose(pFile);
    }

    /* write table to file */
    void write_to_file(const std::vector<double> &table, const char *filename) 
    {
        FILE *pFile;
        pFile = fopen(filename, "w");

        for(std::vector<double>::const_iterator it = table.begin(); 
                it != table.end(); ++it)
        {
            fprintf(pFile, "%16.9e\n", *it);
        }
        fclose(pFile);
        return;
    }

    void write_to_file(const std::vector<std::vector<double> > &table, const char *filename) 
    {
        FILE *pFile;
        pFile = fopen(filename, "w");

        for(auto it = table.cbegin(); 
                it != table.cend(); ++it)
        {
            for(auto it2 = (*it).cbegin();
                    it2 != (*it).cend(); ++it2)
            {
                fprintf(pFile, "%16.9e ", *it2);
            }
            fprintf(pFile, "\n");
        }
        fclose(pFile);
        return;
    }

    void write_to_file(const double *array, const int n_row,
            const int n_col, const char *filename) 
    {
        // write 2d row-major array to file
        //
        //if(n_col * n_row != sizeof(array)/sizeof(*array))
        //{
        //std::cout << sizeof(array)/sizeof(double) << std::endl;
        //fprintf(stderr, "ERROR: Array size does not match supplied number of columns/rows!\n");
        //exit(1);
        //}
        FILE *pFile;
        pFile = fopen(filename, "w");

        int index = 0;

        for(int i_row = 0; i_row < n_row; ++i_row)
        {
            for(int i_col = 0; i_col < n_col; ++i_col)
            {
                fprintf(pFile, "%16.9e ", array[index++]);
            }
            fprintf(pFile, "\n");
        }
        fclose(pFile);
    }

    void write_to_file(const Real3D *array, const int n_row,
            const int n_col, const char *filename) 
    {
        // write 2d row-major array to file
        //
        //if(n_col * n_row != sizeof(array)/sizeof(*array))
        //{
        //std::cout << sizeof(array)/sizeof(double) << std::endl;
        //fprintf(stderr, "ERROR: Array size does not match supplied number of columns/rows!\n");
        //exit(1);
        //}
        FILE *pFile;
        pFile = fopen(filename, "w");

        int index = 0;

        for(int i_row = 0; i_row < n_row; ++i_row)
        {
            for(int i_col = 0; i_col < n_col; ++i_col)
            {
                Real3D val = array[index++];
                fprintf(pFile, "%16.9e %16.9e %16.9e ", val[0], val[1], val[2]);
            }
            fprintf(pFile, "\n");
        }
        fclose(pFile);
    }

    int next_lower_power_of_two(int number)
    {
        /* returns the biggest power of two, that is smaller than number */
        int test = 1;
        while(test < number)
        {
            test *= 2;
        }
        return test / 2;
    }

    // parse argv
    const char *parse_arguments(const int argc, char **argv, const char *indicator, const char *def="0")
    {
        int iarg = 1;
        while(iarg < argc)
        {
            if(strcmp(argv[iarg], indicator) == 0)
            {
                // for indicators without value return 1
                if(iarg == argc - 1 || argv[iarg + 1][0] == '-')
                {
                    return "1";
                }
                // if indicator found, return corresponding value
                return argv[iarg + 1];
            }
            ++iarg;
        }
        // if nothing was found return 0
        return def;
    }

    void print_arguments(int argc, char **argv)
    {
        int i;

        for(i = 0; i < argc; ++i)
        {
            printf("%s ", argv[i]);
        }
        printf("\n");
    }

    bool is_given(const int argc, char **argv, const char *indicator)
    {
        bool tag = false;
        int iarg = 1;
        while(iarg < argc) if(strcmp(argv[iarg++], indicator) == 0) tag = true;
        return tag;
    }

    //double minimum(double x1, double x2)
    //{
        //return x1 < x2 ? x1 : x2;
    //}

    //double maximum(double x1, double x2)
    //{
        //return x1 > x2 ? x1 : x2;
    //}

    /* compare functions for sorting with e.g. qsort()*/
    static int compare(const void * a, const void * b)
    {
        if (*(double*)a > *(double*)b) return 1;
        else if (*(double*)a < *(double*)b) return -1;
        else return 0;
    }

    //static int vcompare_full(const std::vector<double> a , const std::vector<double>  b)
    //{
    //if (a[0] > b[0]) return 1;
    //else if (a[0] < b[0]) return -1;
    //else return 0;
    //}

    // compare vectors by first element
    bool vcompare (const std::vector<double> a , const std::vector<double>  b)
    {
        return (a[0] < b[0]);
    }

    /* sum over exponentials in logarithmic representation */
    double logsumexp(double *exponents, int length)
    {
        int i;
        double delta;
        long double q = 0.;

        /* sort array in ascending order */
        qsort(exponents, length, sizeof(double), compare);

        for(i = 0; i < length - 1; ++i)
        {
            delta = exponents[i] - exponents[length - 1];
            if(delta > - DBL_MAX)
            {
                q += expl(delta);
            }
        }

        return exponents[length - 1] + log1pl( q );
    }

    /* convert rad to angle (circle with circumference length) */
    double to_angle(double length, double rad)
    {
        double angle;

        angle = 360. * rad / length;
        while(angle < 0.)
        {
            angle += 360.;
        }
        while(angle >= 360.)
        {
            angle -= 360.;
        }

        return angle;
    }

    /* round properly to next integer */
    int round_positive_nymnber(double dbl)
    {
        return (int) (dbl + .5);
    }

    int round_negative_number(double dbl)
    {
        return (int) (dbl - .5);
    }

    //int round(double dbl)
    //{
    //if(dbl >= 0) return (int) (dbl + .5);
    //else return (int) (dbl - .5);
    //}

    /* solve LSE via Gauss Algorithm */
    //const int n_points = 10;
    //void gsolve(double m[][n_points - 1], double *b, double *x, int n_points)
    //{
        //int step, i, j, k, lne;
        //double factor, temp;

        //for(step = 0; step < n_points; ++step)
        //{
            //[> move rows with leading zeros to the beginning <]
            //for(i = 0; i < n_points; ++i)
            //{
                //if(m[i][step] == 0)
                //{
                    //[> for every column... <]
                    //for(j = step; j < n_points; ++j)
                    //{
                        //temp = m[i][j];

                        //[> ...move preceding values in said column one up <]
                        //for(k = i; k > 0; --k)
                        //{
                            //m[k][j] = m[k - 1][j];
                        //}
                        //m[0][j] = temp;
                    //}
                    //temp = b[i];

                    //for(k = i; k > 0; --k)
                    //{
                        //b[k] = b[k - 1];
                    //}
                    //b[0] = temp;
                //}
            //}

            //[> produce triangular form <]
            //for(i = 0; i < n_points - step - 1; ++i)
            //{
                //if(m[i][step] != 0)
                //{
                    //factor = m[i][step] / m[i + 1][step];

                    ////m[i][step] = 0.; <- not really necessary
                    ////step + 1
                    //for(j = step; j < n_points; ++j)
                    //{
                        //m[i][j] -= factor * m[i + 1][j];
                    //}
                    //b[i] -= factor * b[i + 1];
                //}
            //}
        //}

        //[> get solutions from triangular form <]
        //for(i = n_points - 1; i >= 0; --i)
        //{
            //lne = n_points - i - 1;
            //x[i] = b[lne];

            //for(j = i + 1; j < n_points; ++j)
            //{
                //x[i] -= m[lne][j] * x[j];
            //}
            //x[i] /= m[lne][i];
        //}
    //}

    void print_vector_matrix(std::vector< std::vector<double> > A)
    {
        int n = A.size();
        for (int i=0; i<n; i++) {
            for (int j=0; j<n+1; j++) {
                std::cout << A[i][j] << "\t";
                if (j == n-1) {
                    std::cout << "| ";
                }
            }
            std::cout << "\n";
        }
        std::cout << std::endl;
    }

    std::vector<double> gauss(std::vector< std::vector<double> > A)
    {
        int n = A.size();

        for (int i=0; i<n; i++) {
            // Search for maximum in this column
            double maxEl = fabs(A[i][i]);
            int maxRow = i;
            for (int k=i+1; k<n; k++) {
                if (fabs(A[k][i]) > maxEl) {
                    maxEl = fabs(A[k][i]);
                    maxRow = k;
                }
            }

            // Swap maximum row with current row (column by column)
            for (int k=i; k<n+1;k++) {
                double tmp = A[maxRow][k];
                A[maxRow][k] = A[i][k];
                A[i][k] = tmp;
            }

            // Make all rows below this one 0 in current column
            for (int k=i+1; k<n; k++) {
                double c = -A[k][i]/A[i][i];
                for (int j=i; j<n+1; j++) {
                    if (i==j) {
                        A[k][j] = 0;
                    } else {
                        A[k][j] += c * A[i][j];
                    }
                }
            }
        }

        // Solve equation Ax=b for an upper triangular matrix A
        std::vector<double> x(n);
        for (int i=n-1; i>=0; i--) {
            x[i] = A[i][n]/A[i][i];
            for (int k=i-1;k>=0; k--) {
                A[k][n] -= A[k][i] * x[i];
            }
        }
        return x;
    }

    /* signum */
    template <typename T> int sgn(T val)
    {
        return (T(0) < val) - (val < T(0));
    }

    /* gets number of characters in buffer */
    unsigned int FileRead(std::istream &is, std::vector<char> &buff) 
    {
        is.read(&buff[0], buff.size());
        return is.gcount();
    }

    /* count lines in a buffer */
    unsigned int CountLines(const std::vector<char> &buff, int sz) 
    {
        int newlines = 0;
        const char *p = &buff[0];
        for (int i = 0; i < sz; i++) {
            if (p[i] == '\n') {
                newlines++;
            }
        }
        return newlines;
    }

    /* count lines in a file */
    int fast_count_lines(const char* filename) 
    {
        const int SZ = 1024 * 1024;

        std::vector <char> buff(SZ);

        std::ifstream ifs(filename);
        int n = 0;
        while(int cc = FileRead(ifs, buff)) {
            n += CountLines(buff, cc);
        }
        return n;
    }

    /* count noncommented, nonempty lines in a file */
    int count_lines(const char* filename)
    {
        int n = 0;
        std::ifstream ifs(filename);
        std::string str;

        while(getline(ifs, str)) 
        {
            if(str[0] == '#' || str[0] == '\0')
            {
                continue;
            }
            else
            {
                ++n;
            }
        }
        return n;
    }

    /* binomial coefficient  */
    template <class T = unsigned long>
        T binomial_coefficient(unsigned long n, unsigned long k) {
            unsigned long i;
            T b;
            if (0 == k || n == k) {
                return 1;
            }
            if (k > n) {
                return 0;
            }
            if (k > (n - k)) {
                k = n - k;
            }
            if (1 == k) {
                return n;
            }
            b = 1;
            for (i = 1; i <= k; ++i) {
                b *= (n - (k - i));
                //if (b < 0) return -1; [> Overflow <]
                b /= i;
            }
            return b;
        }

    /* average of a vector */
    template <typename T>
        T mean(std::vector<T> v)
        {
            return std::accumulate(v.begin(), v.end(), .0) / v.size();
        }

    template <typename T>
        T stddev(std::vector<T> v)
        {
            T sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), .0);
            return std::sqrt(sq_sum / v.size() - mean(v) * mean(v));
        }
}

#endif
