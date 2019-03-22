#include <iomanip>
#include "trajectory.hpp"
#include "couf.hpp"
#include <algorithm>  // std::sort

using namespace std;

int main(int argc, char **argv)
{
    if(argc < 2 || strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0)
    {
        cout << "Usage: " << endl;
        cout << argv[0] << endl;
        cout << "  <char*> infile" << endl;
        cout << "  --step   <int>    (opt) only evaluate every step-th frame" << endl;
        cout << "  --offset <int>    (opt) number of lines to skip at the beginning" << endl;
        cout << "  --max    <int>    (opt) maximum frame to read" << endl;
        exit(1);
    }

    // infile
    const char *infile  = argv[1];
    // frames to skip at the beginning
    const int offset    = atoi(couf::parse_arguments(argc, argv, "--offset"));
    size_t max          = static_cast<size_t>(atoi(couf::parse_arguments(argc, argv, "--max")));
    int step            = atoi(couf::parse_arguments(argc, argv, "--step"));

    if(step == 0) step = 1;
    if(offset) cout << "Skipping " << offset << " frames." << endl;
    if(!couf::is_given(argc, argv, "--max"))
        max = std::numeric_limits<int>::max();

    char outfile[256];

    Trajectory traj(infile);
    traj.advance(offset);

    /* loop through frames */
    while(!traj.is_null() && traj.index() <= max)
    {
        /* file output */
        sprintf(outfile, "multiplied_frame_%05zd.xyz", traj.index());
        traj->multiply().write_xyz(outfile);

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

    exit(0);
}

