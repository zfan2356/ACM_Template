```c++
char c[N][N];
int n, m;
ull pow1[N], pow2[N];
ull hash1[N][N];
void init() {
    pow1[0] = pow2[0] = 1;
    for(int i = 1; i < N; i ++ ) {
        pow1[i] = pow1[i - 1] * P1;
        pow2[i] = pow2[i - 1] * P2;
    }

    for(int i = 1; i <= n; i ++ )
        for(int j = 1; j <= m; j ++ ) {
            hash1[i][j] = hash1[i - 1][j] * P1 + hash1[i][j - 1] * P2 - hash1[i - 1][j - 1] * P1 * P2 + c[i][j];
        }
}

ull get(int lx, int ly, int rx, int ry) {
    return hash1[rx][ry] - hash1[lx - 1][ry] * pow1[rx - lx + 1] - hash1[rx][ly - 1] * pow2[ry - ly + 1] + hash1[lx - 1][ly - 1] * pow1[rx - lx + 1] * pow2[ry - ly + 1];
}
```
