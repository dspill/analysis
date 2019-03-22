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
        cout << "  --ppm     <int> particles per molecule" << '\n';
        cout << "  --step    <int> (opt) only evaluate every step-th frame" << '\n';
        cout << "  --offset  <int> (opt) number of lines to skip at the beginning" << '\n';
        exit(1);
    }

    // infile
    const char *infile {argv[1]};
    size_t particles_per_molecule {static_cast<size_t>(
            atoi(couf::parse_arguments(argc, argv, "--ppm")))};
    size_t step        {static_cast<size_t>(
            atoi(couf::parse_arguments(argc, argv, "--step", "1")))};
    size_t max         {static_cast<size_t>(atoi(couf::parse_arguments(argc, argv, "--max")))};
    const int offset   {atoi(couf::parse_arguments(argc, argv, "--offset", "0"))};

    if(particles_per_molecule == 0) throw runtime_error("You have to give the chain length.");
    if(offset) cout << "Skipping " << offset << " frames." << '\n';
    if(!couf::is_given(argc, argv, "--max"))
        max = std::numeric_limits<int>::max();
    /* ----------- input done ---------- */

    // open infile
    Trajectory traj{infile, particles_per_molecule};
    traj += offset;

    cout << "calculating end-to-end...";
    cout.flush();
    auto ts = traj.timeseries(&Molecule::end_to_end, step, max);
    const double norm = ts.autocorrelation_function(0);
    cout << "done.\n";

    // open outfile
    const char outfile_name[256] = "autocorrelation_function.dat";
    remove(outfile_name);
    FILE *outfile;
    outfile = fopen(outfile_name, "w");
    fprintf(outfile, "# step C_ete C_rcm\n");
    for(size_t span = 0; span < ts.size(); ++span)
    {
        fprintf(outfile, "%4zd %16.9e\n",
                span * step, ts.autocorrelation_function(span)/norm);
    }
    fclose(outfile);

    exit(0);
}
