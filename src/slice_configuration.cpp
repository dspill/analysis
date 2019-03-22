#include <iomanip>
#include "trajectory.hpp"
#include "couf.hpp"

using namespace std;

int main(int argc, char **argv)
{
    if(argc < 2 || strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0)
    {
        cout << "Usage: " << '\n';
        cout << argv[0] << '\n';
        cout << "  <char*> infile" << '\n';
        cout << "  -d       <double> thickness of the slice" << '\n';
        cout << "  --of     <char*>  (opt) outfile name" << '\n';
        cout << "  --ppm    <int>    particles per molecule" << '\n';
        cout << "  --step   <int>    (opt) only evaluate every step-th frame" << '\n';
        cout << "  --offset <int>    (opt) number of lines to skip at the beginning" << '\n';
        cout << "  --max    <int>    (opt) maximum frame to read" << '\n';
        exit(1);
    }

    // infile
    const char *infile  = argv[1];
    // frames to skip at the beginning
    double thickness = atoi(couf::parse_arguments(argc, argv, "-d"));
    size_t particles_per_molecule {static_cast<size_t>(
            atoi(couf::parse_arguments(argc, argv, "--ppm")))};

    if(thickness == 0) thickness = .75;

    char outfile[256];
    sprintf(outfile, "sliced_%s", infile);
    remove(outfile);

    Trajectory traj(infile, particles_per_molecule);

    /* loop through frames */
    while(!traj.is_null())
    {
        printf("\rprocessing frame %zd ", traj.index());
        cout.flush();

        /* file output */
        //traj->set_types(128);
        traj->set_types();
        //traj->slice(thickness).write_xyz(outfile, true);
        traj->slice_square(thickness).write_xyz(outfile, true);

        /* advance to next frame */
        traj.loop_advance(argc, argv);
    }
    exit(0);
}

