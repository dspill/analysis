#include <iomanip>
#include "trajectory.hpp"
#include "couf.hpp"
#include "limits.h"
#include "fftw3.h"
#include <algorithm>   // sort
#include <sys/stat.h>  // mkdir

#define QMAX 10.

using namespace std;

int main(int argc, char **argv)
{
    if(argc < 2 || strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0)
    {
        cout << "Usage: " << '\n';
        cout << argv[0] << '\n';
        cout << " <char*> infile_1" << '\n';
        cout << " <char*> infile_2" << '\n';
        cout << " <char*> outfile" << '\n';
        exit(1);
    }

    // infile
    const char* infile_1  = argv[1];
    const char* infile_2  = argv[2];
    // outfile
    const char* outfile   = argv[3];
    /* ----------- input done ---------- */

    Trajectory t1(infile_1);
    Trajectory t2(infile_2);

    t1.write_xyz(outfile, false);
    if(*t1 == *t2)
    {
        cout << "Last frame of " << infile_1 
            << " coincides with first frame of " << infile_2 << ", skipping\n";
        ++t2;
    }
    if(t2.is_good()) t2.write_xyz(outfile, true);

    exit(0);
}
