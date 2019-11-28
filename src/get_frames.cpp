#include <iomanip>
#include "trajectory.hpp"
#include "couf.hpp"
#include <algorithm>  // std::sort

using namespace std;

int main(int argc, char **argv)
{
    if(argc < 2 || strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0)
    {
        cout << "Usage:\n";
        cout << argv[0] << '\n';
        cout << "  <char*> infile\n";
        cout << "  --step      <int>    (opt) only evaluate every step-th frame\n";
        cout << "  --offset    <int>    (opt) number of lines to skip at the "
            "beginning\n";
        cout << "  --max       <int>    (opt) maximum frame to read\n";
        cout << "  --append    flag     (opt) output to ONE file\n";
        cout << "  --vtk       flag     (opt) output as vtk\n";
        cout << "  --pdb       flag     (opt) output as pdb\n";
        cout << "  --ppm       <int>    particles_per_molecule (needed for vtk)\n";
        cout << "  --set_types <int>    deduct particle types from indices\n";
        cout << "  --multiply  <int>    make new system 2x2x2 times the size\n";
        cout << "  --slice     <double>\n";
        cout << "  --exp       <double>\n";
        exit(1);
    }

    // infile
    const char *infile   = argv[1];
    // frames to skip at the beginning
    const bool append    = atoi(couf::parse_arguments(argc, argv, "--append"));
    const bool write_vtk = atoi(couf::parse_arguments(argc, argv, "--vtk"));
    const bool write_pdb = atoi(couf::parse_arguments(argc, argv, "--pdb"));
    const int  set_types = atoi(couf::parse_arguments(argc, argv, "--set_types"));
    const int  multiply  = atoi(couf::parse_arguments(argc, argv, "--multiply"));
    const double  slice  = atoi(couf::parse_arguments(argc, argv, "--slice", "0"));
    size_t particles_per_molecule {static_cast<size_t>(
            atoi(couf::parse_arguments(argc, argv, "--ppm")))};

    if(write_vtk || write_pdb || set_types || multiply)
        if(particles_per_molecule == 0)
            throw runtime_error("You have to give the number of particles per "
                    "molecule.");

    char outfile[256];

    /* loop through frames */
    Trajectory traj{infile, particles_per_molecule};
    Frame frame;
    do {
        //printf("\rReading frame %zd ", traj.index());
        //cout.flush();

        /* file output */
        frame = *traj;
        if(multiply) frame = frame.multiply();
        if(set_types) frame.set_types(set_types);
        if(slice) frame = frame.slice_square(slice);
        if(write_vtk)
        {
            if(append) sprintf(outfile, "traj_filtered.vtk");
            else sprintf(outfile, "frame_%05zd.vtk", traj.index());
            frame.write_vtk(outfile, append);
        }
        else if(write_pdb)
        {
            if(append) sprintf(outfile, "traj_filtered.pdb");
            else sprintf(outfile, "frame_%05zd.pdb", traj.index());
            frame.set_types();
            frame.write_pdb(outfile, append);
        }
        else
        {
            if(append) sprintf(outfile, "traj_filtered.xyz");
            else sprintf(outfile, "frame_%05zd.xyz", traj.index());
            frame.write_xyz(outfile, append);
        }

        /* advance to next frame */
        traj.loop_advance(argc, argv);
    } while(!traj.is_null());

    exit(0);
}

