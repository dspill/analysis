#include "trajectory.hpp"
#include "couf.hpp"
#include <iterator>

using namespace std;

int main(int argc, char **argv)
{
    if(argc < 2 || strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0)
    {
        cout << "Usage: " << '\n';
        cout << argv[0] << '\n';
        cout << "  <char*> infile" << '\n';
        cout << "  --bw     <double> (opt) bin width" << '\n';
        cout << "  --qmax   <double> (opt) largest wave-vector to consider" << '\n';
        cout << "  --nrand  <int>    (opt) amount of random numbers to use in computation of S(q)" << '\n';
        cout << "  --step   <int>    (opt) only evaluate every step-th frame" << '\n';
        cout << "  --offset <int>    (opt) number of lines to skip at the beginning" << '\n';
        cout << "  --max    <int>    (opt) highest number opf frame to read" << endl;
        exit(1);
    }

    // infile
    const char *infile  = argv[1];
    // frames to skip at the beginning
    const int offset    = atoi(couf::parse_arguments(argc, argv, "--offset"));
    size_t max          = static_cast<size_t>(atoi(couf::parse_arguments(argc, argv, "--max")));
    double bin_width    = atof(couf::parse_arguments(argc, argv, "--bw"));
    double qmax         = atof(couf::parse_arguments(argc, argv, "--qmax", "7"));
    int nrand           = atoi(couf::parse_arguments(argc, argv, "--nrand", "128"));
    int step            = atoi(couf::parse_arguments(argc, argv, "--step", "1"));

    constexpr size_t precision = 11;
    char outfile[256];

    /* ----------- input done ---------- */
    Trajectory traj{infile};
    const double q_min = 2.*M_PI / std::min(traj->box()[0], 
            std::min(traj->box()[1], traj->box()[2]));
    if(bin_width == 0.0) bin_width = q_min;
    if(bin_width < q_min)
        throw runtime_error("Bin-size is too small");

    cout << "infile:    " << infile << '\n';
    cout << "offset:    " << offset << '\n';
    cout << "max:       " << max << '\n';
    cout << "bin_width: " << bin_width << '\n';
    cout << "step:      " << step << '\n';

    traj += offset;

    /* loop through frames */
    while(!traj.is_null())
    {
        cout << "reading frame " << traj.index() << '\n';
        cout.flush();

        sprintf(outfile, "bf_dsf_%05zd.dat", traj.index());
        ofstream stream(outfile);
        stream.precision(precision);
        stream << std::scientific;

        /* compute structure factor */
        for(double q = 0.; q < qmax; q += bin_width)
        {
            cout << "computing q " << q << '\n';
            auto result =  mean_structure_factor(*traj, q, nrand);
            copy(result.begin(), result.end(), ostream_iterator<double>(stream, " "));
            stream << "\n";
        }
        stream.close();

        /* advance to next frame */
        traj.loop_advance(argc, argv);
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
