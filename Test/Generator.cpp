#include <bits/stdc++.h>

using namespace std;

using i64 = long long;

std::mt19937 seed(time(NULL));
template <class T> T rand(T l, T r) {
    return std::uniform_int_distribution<T>(l, r)(seed);
}


int main() {
    int n = rand(1, 100), t = rand(0, 1), k = rand(1, 100);
    cout << n << " " << k << " " << t << '\n';
    for (int i = 0; i < n; i++) {
        cout << rand(1, 10) << " ";
    }
    cout << '\n';

    return 0;
}