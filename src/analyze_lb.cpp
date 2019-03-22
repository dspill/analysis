#include "lb_frame.hpp"

using namespace std;

int main(int argc, char ** argv)
{
    LB_Frame f(20);
    cout << f;

    int i = 1;
    while (i < argc && argv[i][0] != '-')
    {
        cout << "adding " << argv[i] << '\n';
        f.add(argv[i]);

        ++i;
    }
    cout << f.density(Index{7, 3, 16}) << '\n';
    cout << f.velocity(Index{7, 3, 16}) << '\n';
    cout << f.total_momentum_density()/f.size() << '\n';
    cout << f.total_density()/f.size() << '\n';

    return argc;
}
