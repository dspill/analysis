#include <iomanip>
#include "trajectory.hpp"
#include "couf.hpp"
#include "limits.h"
#include "fftw3.h"
#include <algorithm>   // sort
#include <sys/stat.h>  // mkdir

#define QMAX 10.

using namespace std;

int main(int argc, char **argv)
{
    if(argc < 2 || strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0)
    {
        cout << "Usage: " << '\n';
        cout << argv[0] << '\n';
        cout << "  <char*> infile" << '\n';
        cout << "  --thr    <double> (opt) threshold" << '\n';
        cout << "  --cg     <int>    (opt) factor by which the lattice is "
            "coarse grained" << '\n';
        cout << "  --fg     <int>    (opt) factor by which the lattice is "
            "fine grained" << '\n';
        cout << "  --step   <int>    (opt) only evaluate every step-th frame"
            << '\n';
        cout << "  --offset <int>    (opt) number of lines to skip at the "
            "beginning" << '\n';
        cout << "  --max    <int>    (opt) highest number of frame to read" 
            << '\n';
        cout << "  --exp    <double>\n";
        cout << "  --natural_units <bool>\n";
        exit(1);
    }

    constexpr size_t precision = 10;

    // infile
    const char *infile  = argv[1];
    // outfile
    char outfile[256];

    // frames to skip at the beginning
    double threshold    = atof(couf::parse_arguments(argc, argv, "--thr", "-1"));
    int cg_factor       = atoi(couf::parse_arguments(argc, argv, "--cg", "1"));
    int fg_factor       = atoi(couf::parse_arguments(argc, argv, "--fg", "1"));
    bool natural_units  = atoi(couf::parse_arguments(argc, argv, "--natural_units", "0"));
    double bin_width    = atof(couf::parse_arguments(argc, argv, "--bw"));


    Trajectory traj{infile};
    const double q_min = 2.*M_PI / std::min(traj->box()[0], 
            std::min(traj->box()[1], traj->box()[2]));
    if(bin_width == 0.0) bin_width = q_min;
    if(bin_width < q_min)
        throw runtime_error("Bin-size is too small");
    /* ----------- input done ---------- */

    cout << "infile:    " << infile << '\n';
    cout << "threshold: " << threshold << '\n';
    cout << "cg_factor: " << cg_factor << '\n';
    cout << "fg_factor: " << fg_factor << '\n';

    /* get lattice size */
    const size_t original_lattice_size = round(traj->box()[0]);

    if(!mkdir("./dsf", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) ||
    !mkdir("./minkowski", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH))
    {
        throw std::runtime_error("could not create directory/directories");
    }

    /* loop through frames */
    double lattice_constant;
    size_t lattice_size;
    while(!traj.is_null())
    {
        cout << "reading frame " << traj.index() << '\n';
        lattice_constant = (double) cg_factor / fg_factor;
        lattice_size = round(original_lattice_size / lattice_constant);
        Frame f = *traj;

        /* compute structure factor */
        vector<vector<double>> struc_fac =
            structure_factor(f, lattice_size, lattice_constant, bin_width, 
                    traj->size());

        sprintf(outfile, "./dsf/dsf_%05zd.dat", traj.index());
        couf::write_to_file(struc_fac, outfile);

        vector<double> thresholds{.6, .7, .8, .9, 1., 1.1, 1.2, 1.3, 1.4};
        vector<double> cgs {1, 2, 4, 8};
        for(const size_t& cg : cgs)
        {
            for(const double& thr : thresholds)
            {
                lattice_constant = (double) cg;
                lattice_size = round(original_lattice_size / lattice_constant);

                sprintf(outfile,
                        "./minkowski/minkowski_functionals_fg%d_cg%zd_thr%3.2f.dat",
                        1, cg, thr);
                if(traj.index() == 0) 
                {
                    remove(outfile);
                }
                ofstream stream(outfile, ofstream::app);

                /* compute minkowski functionals */
                array<double, 6> mfs = minkowski_functionals(
                        f, lattice_size, thr, 'c', natural_units);

                /* file output */
                stream << fixed;
                stream << scientific << setw(6);
                stream << traj.index();

                stream.precision(precision);
                for(auto mf : mfs)
                {
                    stream << scientific << setw(precision + 8) << mf;
                }
                //stream << '\n';
                stream << endl;
                stream.close();
            }
        }

        /* advance to next frame */
        traj.loop_advance(argc, argv);
    }
    exit(0);
}
