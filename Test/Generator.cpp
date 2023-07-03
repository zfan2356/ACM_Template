#include "testlib.h"
#include <bits/stdc++.h>

using namespace std;

int main(int argc, char* argv[]) {
    registerGen(argc, argv, 1);

    int n = 7, t = 50;
    int q = rnd.next(1, 5);

    vector<int> p(n);

/* 为节点 1..n-1 设置父亲 */
    for (int i = 0; i < n; i++) {
        if (i > 0) {
            p[i] = rnd.wnext(i, t);
        }
    }

    printf("%d %d\n", n, q);

/* 打乱节点 1..n-1 */
    vector<int> perm(n);
    for (int i = 0; i < n; i++) {
        perm[i] = i;
    }
    shuffle(perm.begin() + 1, perm.end());

/* 根据打乱的节点顺序加边 */
    vector<pair<int, int> > edges;
    for (int i = 1; i < n; i++)
        if (rnd.next(2))
            edges.push_back(make_pair(perm[i], perm[p[i]]));
        else
            edges.push_back(make_pair(perm[p[i]], perm[i]));

/* 打乱边 */
    shuffle(edges.begin(), edges.end());

    for (int i = 0; i + 1 < n; i++)
        printf("%d %d %d %d\n", edges[i].first + 1, edges[i].second + 1, rnd.next(1, n - 1), rnd.next(1, 6));

    for (int i = 0; i < q; i++) {
        int u = rnd.next(1, n), v = rnd.next(1, n);
        if (u == v) {
            while (u == v) {
                v = rnd.next(1, n);
            }
        }
        if (u > v) {
            swap(u, v);
        }
        printf("%d %d %d %d\n", rnd.next(1, n), rnd.next(1, 6), u, v);
    }
    return 0;
}