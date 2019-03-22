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
        cout << "  -o        <char*>  (opt) output file" << '\n';
        cout << "  -c        <size_t> (opt) column that is to be read" << '\n';
        cout << "  --offset  <int>    (opt) number of lines to skip at the beginning" << '\n';
        //cout << "  --step    <int>   (opt) only evaluate every step-th frame" << '\n';
        //cout << "  --dt      <double>(opt) timestep" << '\n';
        //cout << "  --max     <int>    (opt) highest number opf frame to read" << endl;
        exit(1);
    }

    // infile
    const char * infile {argv[1]};
    const char * outfile = couf::parse_arguments(argc, argv, "-o");
    size_t column        {static_cast<size_t>(
            atoi(couf::parse_arguments(argc, argv, "-c")))};
    const int offset   {atoi(couf::parse_arguments(argc, argv, "--offset"))};
    //size_t step        {static_cast<size_t>(
            //atoi(couf::parse_arguments(argc, argv, "--step")))};
    //double timestep    {atof(couf::parse_arguments(argc, argv, "--dt"))};
    //size_t max         {static_cast<size_t>(atoi(couf::parse_arguments(argc, argv, "--max")))};

    if(outfile[0] == '0') outfile = "ana_timeseries.dat";
    if(column == 0) column = 1;
    //if(step == 0) step = 1;
    //if(timestep == 0.) timestep = 1.;
    //if(offset) cout << "Skipping " << offset << " frames." << '\n';
    //if(!couf::is_given(argc, argv, "--max")) max = std::numeric_limits<size_t>::max();
    /* ----------- input done ---------- */

    // open infile
    Timeseries<double> ts(infile, column, offset);
    cout.precision(10);
    cout << scientific << ts.mean() << ' ' << ts.stdev() << '\n';

    exit(0);
}
