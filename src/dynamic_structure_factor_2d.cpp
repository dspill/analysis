#include <iomanip>
#include <array>
#include "trajectory.hpp"
#include "couf.hpp"
#include "fftw3.h"
#include "limits.h"
#include <algorithm>    // std::sort

#define KMAX 10.

using namespace std;

int main(int argc, char **argv)
{
    //cout << Frame::minimal_distance(Real3D(1.,0,-1), Real3D(9.,0,4), Real3D(10)) << endl;
    //cout << Frame::fold(Real3D(10.,5.,-2.), Real3D(10)) << endl;
    //exit(1);

    if(argc < 2 || strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0)
    {
        cout << "Usage: " << endl;
        cout << argv[0] << endl;
        cout << "  <char*> infile" << endl;
        cout << "  --dt     <double> (opt) md timestep" << endl;
        cout << "  --bw     <double> (opt) bin width" << endl;
        cout << "  --cg     <int>    (opt) factor by which the lattice is coarse grained" << endl;
        cout << "  --fg     <int>    (opt) factor by which the lattice is fine grained" << endl;
        cout << "  --step   <int>    (opt) only evaluate every step-th frame" << endl;
        cout << "  --offset <int>    (opt) number of lines to skip at the beginning" << endl;
        cout << "  --max    <int>    (opt) highest number opf frame to read" << endl;
        exit(1);
    }

    // infile
    const char *infile  = argv[1];
    // frames to skip at the beginning
    const int offset    = atoi(couf::parse_arguments(argc, argv, "--offset"));
    int max       = atoi(couf::parse_arguments(argc, argv, "--max"));
    double timestep     = atof(couf::parse_arguments(argc, argv, "--dt"));
    double bin_width    = atof(couf::parse_arguments(argc, argv, "--bw"));
    int cg_factor       = atoi(couf::parse_arguments(argc, argv, "--cg"));
    int fg_factor       = atoi(couf::parse_arguments(argc, argv, "--fg"));
    int step            = atoi(couf::parse_arguments(argc, argv, "--step"));

    if(timestep == 0.0) timestep = 1;
    if(bin_width == 0.0) bin_width = 0.01;
    if(cg_factor == 0) cg_factor = 1;
    if(fg_factor == 0) fg_factor = 1;
    if(max == 0) max = std::numeric_limits<int>::max();;

    if(step == 0) step = 1;
    if(offset) cout << "Skipping " << offset << " frames." << endl;

    const double lattice_constant = (double) cg_factor / fg_factor;

    char outfile[256];

    Trajectory traj(infile);

    /* get lattice size and allocate arrays */
    const int original_lattice_size = round(traj->box()[0]);
    const int lattice_size = round(original_lattice_size / lattice_constant);
    const int number_of_sites = lattice_size * lattice_size;
    if(lattice_size % cg_factor != 0)
    {
        fprintf(stderr, "ERROR: lattice_size must be integer multiple of cg_factor\n");
        exit(1);
    }

    /* allocate arrays */
    fftw_complex *lattice_transformed =
        fftw_alloc_complex(lattice_size * (lattice_size / 2 + 1) * sizeof(fftw_complex));
    double *lattice  = fftw_alloc_real(number_of_sites * sizeof(double));

    /* generate fftw plan */
    fftw_plan plan = fftw_plan_dft_r2c_2d(lattice_size, lattice_size, lattice,
            lattice_transformed, FFTW_ESTIMATE);

    cout << "infile                = " << infile << '\n';
    cout << "bin_width             = " << bin_width << '\n';
    cout << "cg_factor             = " << cg_factor << '\n';
    cout << "fg_factor             = " << fg_factor << '\n';
    cout << "original lattice size = " << original_lattice_size << '\n';
    cout << "scaled lattice size   = " << lattice_size << '\n';
    cout << "lattice constant      = " << lattice_constant << '\n';

    while(!traj.is_null())
    {
        printf("Reading frame %zd\n", traj.index());

        // reset lattice
        for(int i = 0; i < number_of_sites; ++i) lattice[i] = 0.;
        //memset(lattice, 0., number_of_sites*sizeof(double));
        read_lattice(*traj, lattice, lattice_size, nullptr, 2);

        /* do the transformation */
        fftw_execute(plan);

        /* format output */
        vector<array<double, 2>> result =
        linearize_lattice(lattice_transformed, lattice_size, lattice_constant,
                bin_width, traj->size(), 2);

        /* file output */
        sprintf(outfile, "dsf_%05zd.dat", traj.index());
        couf::write_to_file(result, outfile);

        /* advance to next frame */
        traj.loop_advance(argc, argv);
    }

    fftw_destroy_plan(plan);
    fftw_free(lattice);
    fftw_free(lattice_transformed);

    exit(0);
}
