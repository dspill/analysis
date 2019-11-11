#include <iostream>

using namespace std;

template<typename T>
void foo(T x, int y=1)
{
    cout << "generic version:\n";
    cout << x << y << '\n';
}

template<>
void foo<int>(int x, int y)
{
    cout << "integer version:\n";
    cout << x << y << '\n';
}

template<>
void foo<string>(string x, int y)
{
    cout << "string version:\n";
    cout << x << y << '\n';
}

int main()
{
    foo("Hello world", 1);
    foo(string("Hello world"), 2);
    foo(123);
}

