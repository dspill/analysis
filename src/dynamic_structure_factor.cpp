#include <iomanip>
#include "trajectory.hpp"
#include "couf.hpp"
#include "fftw3.h"
#include "limits.h"
#include <algorithm>    // std::sort

#define QMAX 10.

using namespace std;

int main(int argc, char **argv)
{
    if(argc < 2 || strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0)
    {
        cout << "Usage: " << '\n';
        cout << argv[0] << '\n';
        cout << "  <char*> infile" << '\n';
        cout << "  --of     <char*>  (opt) outfile name" << '\n';
        cout << "  --bw     <double> (opt) bin width" << '\n';
        cout << "  --cg     <int>    (opt) factor by which the lattice is coarse grained" << '\n';
        cout << "  --fg     <int>    (opt) factor by which the lattice is fine grained" << '\n';
        cout << "  --step   <int>    (opt) only evaluate every step-th frame" << '\n';
        cout << "  --offset <int>    (opt) number of lines to skip at the beginning" << '\n';
        cout << "  --max    <int>    (opt) highest number of frame to read" << '\n';
        cout << "  --exp    <double>\n";
        cout << "  --dim    <int>    (3)   spacial dimension\n";
        exit(1);
    }

    // infile
    const char *infile  = argv[1];
    // outfile
    char outfile[256];

    // frames to skip at the beginning
    double bin_width    = atof(couf::parse_arguments(argc, argv, "--bw", "0"));
    int cg_factor       = atoi(couf::parse_arguments(argc, argv, "--cg", "1"));
    int fg_factor       = atoi(couf::parse_arguments(argc, argv, "--fg", "1"));
    size_t dim          = static_cast<size_t>(
            atoi(couf::parse_arguments(argc, argv, "--dim", "3")));

    const double lattice_constant = (double) cg_factor / fg_factor;
    /* ----------- input done ---------- */

    Trajectory traj{infile};

    /* get lattice size */
    const size_t original_side_length = round(traj->box()[0]);
    const size_t side_length = round(original_side_length / lattice_constant);

    cout << "infile                = " << infile << '\n';
    cout << "bin_width             = " << bin_width << '\n';
    cout << "cg_factor             = " << cg_factor << '\n';
    cout << "fg_factor             = " << fg_factor << '\n';
    cout << "original lattice size = " << original_side_length << '\n';
    cout << "scaled lattice size   = " << side_length << '\n';
    cout << "lattice constant      = " << lattice_constant << '\n';
    if(side_length % cg_factor != 0)
        throw runtime_error("side_length must be divisible by cg_factor.\n");

    /* loop through frames */
    while(!traj.is_null())
    {
        cout << "reading frame " << traj.index() << '\n';
        /* compute structure factor */
        vector<array<double, 2>> struc_fac =
            structure_factor(*traj, side_length, lattice_constant, bin_width, 
                    traj->size(), dim);

        /* file output */
        if(couf::is_given(argc, argv, "--of"))
            strcpy(outfile, couf::parse_arguments(argc, argv, "--of"));
        else
        {
            sprintf(outfile, "dsf%zdd_%05zd.dat", dim, traj.index());
        }
        couf::write_to_file(struc_fac, outfile);

        /* advance to next frame */
        traj.loop_advance(argc, argv);
    }

    exit(0);
}
