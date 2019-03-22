#include "timeseries.hpp"

using namespace std;

int main(int argc, char **argv)
{
    if(argc < 2 || strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0)
    {
        cout << "Usage: " << '\n';
        cout << argv[0] << '\n';
        cout << "  <char*> space separated list of input files" << '\n';
        cout << "  -o    <char*> output file" << '\n';
        cout << "  --cg  <int>   coarse-graining factor" << '\n';
        exit(1);
    }
    const char * outfile = 
        couf::parse_arguments(argc, argv, "-o", "averaged_timeseries.dat");
    const size_t cg = static_cast<size_t>(
            atoi(couf::parse_arguments(argc, argv, "--cg", "1")));

    if(cg != 1) cout << "Coarse-graining by a factor " << cg << '\n';

    Timeseries<vector<double>> ts1;
    Timeseries<vector<double>> ts2;

    int i = 1;
    cout << "reading " << argv[i] << '\n';
    ts1.read_all(argv[i++]);
    exit(1);

    size_t n_files = 1;
    while (i < argc && argv[i][0] != '-')
    {
        cout << "reading " << argv[i] << '\n';
        ts2.read_all(argv[i++]);
        if(ts1.size() != ts2.size())
        {
            throw runtime_error("cannot average over timeseries of unequal length");
        }
        //ts1 += ts2;
        //cout << ts2[0] << '\n';
        ts2.clear();
        ++n_files;
    }
    cout << "did average over " << n_files << " files\n";

    ts1 /= n_files;
    ts1.coarsen(cg);
    ts1.write(outfile);

    exit(0);
}
