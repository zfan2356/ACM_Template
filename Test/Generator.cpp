#include <bits/stdc++.h>

using namespace std;

int main() {
    srand(time(NULL));
    int n = 10;
    cout << n << endl;
    for (int i = 0; i < n; i++) {
        cout << rand() % 1000 << " ";
    }
    cout << endl;


    return 0;
}