#include <bits/stdc++.h>
#include "testlib.h"

using namespace std;


int main(int argc, char *argv[]) {
    registerGen(argc, argv, 1);

    int T = 1;
    cout << T << endl;
    int n = 100000, m = 10000, s = rnd.next(1, n);
    cout << n << " " << m << " " << s << endl;
    for (int i = 0; i < n; i++) {
        cout << rnd.next(1, 1000000) << " " << rnd.next(1, 1000000) << " ";
        int l = rnd.next(1, n), r = rnd.next(1, n);
        while (l == r) {
            l = rnd.next(1, n), r = rnd.next(1, n);
        }
        if (l > r) {
            swap(l, r);
        }
        cout << l << " " << r << endl;
    }

    while (m--) {
        int l = rnd.next(1, n), r = rnd.next(1, n);
        while (l == r) {
            l = rnd.next(1, n), r = rnd.next(1, n);
        }
        if (l > r) {
            swap(l, r);
        }
        cout << l << " " << r << endl;
    }

    return 0;
}