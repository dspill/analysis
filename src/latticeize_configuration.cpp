#include <iomanip>
#include "trajectory.hpp"
#include "couf.hpp"
#include "fftw3.h"
#include "limits.h"

using namespace std;

int main(int argc, char **argv)
{

    if(argc < 2 || strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0)
    {
        cout << "Usage: " << endl;
        cout << argv[0] << endl;
        cout << "  <char*> infile" << endl;
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
    size_t max          = static_cast<size_t>(atoi(couf::parse_arguments(argc, argv, "--max")));
    int cg_factor       = atoi(couf::parse_arguments(argc, argv, "--cg"));
    int fg_factor       = atoi(couf::parse_arguments(argc, argv, "--fg"));
    int step            = atoi(couf::parse_arguments(argc, argv, "--step"));

    if(cg_factor == 0) cg_factor = 1;
    if(fg_factor == 0) fg_factor = 1;
    if(step == 0) step = 1;
    if(offset) cout << "Skipping " << offset << " frames." << endl;
    if(!couf::is_given(argc, argv, "--max"))
            max = std::numeric_limits<int>::max();

    constexpr size_t dim = 3;
    const double lattice_constant = (double) cg_factor / fg_factor;
    char outfile[256];

    Trajectory traj(infile);
    traj.advance(offset);

    /* get lattice size and allocate arrays */
    const int original_lattice_size = round(traj->box()[0]);
    const int lattice_size = round(original_lattice_size / lattice_constant);
    const int number_of_sites = pow(lattice_size, dim);
    cout << "original lattice size = " << original_lattice_size << endl;
    cout << "scaled lattice size   = " << lattice_size << endl;
    cout << "lattice constant      = " << lattice_constant << endl;
    if(original_lattice_size % cg_factor != 0)
    {
        fprintf(stderr, "ERROR: lattice_size must be integer multiple of cg_factor\n");
        exit(1);
    }

    /* allocate arrays */
    double *lattice          = new double[number_of_sites]; 
    Real3D *velocity_lattice = new Real3D[number_of_sites];

    /* loop through frames */
    while(!traj.is_null() && traj.index() <= max)
    {

        printf("\rReading frame %zd ", traj.index());
        cout.flush();

        // reset lattice
        for(int i = 0; i < number_of_sites; ++i) lattice[i] = 0.;
        for(int i = 0; i < number_of_sites; ++i) velocity_lattice[i] = 0.;
        read_lattice(*traj, lattice, lattice_size, velocity_lattice);

        /* file output */
        sprintf(outfile, "lattice_configuration_%04zd.dat", traj.index());
        couf::write_3d_array_to_file(lattice, outfile, lattice_size);

        sprintf(outfile, "lattice_velocities_%04zd.dat", traj.index());
        couf::write_3d_Real3D_array_to_file(velocity_lattice, outfile, lattice_size);
        
        /* advance to next frame */
        traj += step;
    }
    int frames_read = traj.index() - offset;
    frames_read = frames_read / step + (frames_read % step != 0);

    printf("\n");
    printf("n_frames    = %zd\n", traj.index());
    printf("step        = %d\n", step);
    printf("offset      = %d\n", offset);
    printf("frames_read = %d\n", frames_read);

    delete [] lattice;
    delete [] velocity_lattice;


    exit(0);
}
