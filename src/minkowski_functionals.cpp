#include <iomanip>
#include "trajectory.hpp"
#include "couf.hpp"
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
        cout << "  --thr    <double> (opt) threshold" << '\n';
        cout << "  --cg     <int>    (opt) factor by which the lattice is coarse grained" << '\n';
        cout << "  --fg     <int>    (opt) factor by which the lattice is fine grained" << '\n';
        cout << "  --step   <int>    (opt) only evaluate every step-th frame" << '\n';
        cout << "  --offset <int>    (opt) number of lines to skip at the beginning" << '\n';
        cout << "  --max    <int>    (opt) highest number of frame to read" << '\n';
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
    double threshold    = atof(couf::parse_arguments(argc, argv, "--thr"));
    int cg_factor       = atoi(couf::parse_arguments(argc, argv, "--cg", "1"));
    int fg_factor       = atoi(couf::parse_arguments(argc, argv, "--fg", "1"));
    bool natural_units  = atoi(couf::parse_arguments(argc, argv, "--natural_units", "0"));

    if(couf::is_given(argc, argv, "--of"))
        strcpy(outfile, couf::parse_arguments(argc, argv, "--of"));
    else
        sprintf(outfile, "minkowski_functionals_fg%d_cg%d_thr%3.2f.dat", fg_factor, cg_factor, threshold);

    const double lattice_constant = (double) cg_factor / fg_factor;

    /* ----------- input done ---------- */

    Trajectory traj{infile};

    cout << "infile:    " << infile << '\n';
    cout << "threshold: " << threshold << '\n';
    cout << "cg_factor: " << cg_factor << '\n';
    cout << "fg_factor: " << fg_factor << '\n';

    /* get lattice size */
    const size_t original_lattice_size = round(traj->box()[0]);
    const size_t lattice_size = round(original_lattice_size / lattice_constant);
    cout << "original lattice size = " << original_lattice_size << '\n';
    cout << "scaled lattice size   = " << lattice_size << '\n';
    cout << "lattice constant      = " << lattice_constant << '\n';
    if(lattice_size % cg_factor != 0)
        throw runtime_error("Lattice_size must be divisible by cg_factor.\n");

    ofstream stream(outfile);
    stream << "# step V_0 V_1 V_2^(4) V_2^(8) V_3^(6) V_3^(26)\n";

    /* loop through frames */
    while(!traj.is_null())
    {
        cout << "reading frame " << traj.index() << '\n';
        /* compute minkowski functionals */
        array<double, 6> mfs =
            minkowski_functionals(*traj, lattice_size, threshold, 's', natural_units);

        /* file output */
        stream << std::fixed;
        stream << std::scientific << std::setw(6);
        stream << traj.index();

        stream.precision(precision);
        for(auto mf : mfs)
        {
            stream << std::scientific << std::setw(precision + 8) << mf;
        }
        stream << '\n';

        /* advance to next frame */
        traj.loop_advance(argc, argv);
    }
    stream.close();
    exit(0);
}
