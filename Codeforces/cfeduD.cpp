// CF edu D（jiangly dp）

#include <bits/stdc++.h>

#define int long long
#define endl '\n'
using namespace std;
using LL = long long;

template<typename A, typename B> inline ostream& operator<< (ostream& out, const pair<A, B>& p) { return out << "(" << p.first << ", " << p.second << ")"; }
template<typename T> inline ostream& operator<< (ostream& out, const vector<T>& a) { out << "["; for(int i = 0; i < a.size(); i ++) { if (i) out << ','; out << ' ' << a[i]; } return out << " ]"; }
template<typename T> inline ostream& operator<< (ostream& out, const set<T>& a) { return out << vector<T>(all(a)); }

const int mod = 998244353, inf = 1e18;

void solve() {
    int n, k, x; cin >> n >> k >> x;
    vector<int> a(n + 1);
    for(int i = 1; i <= n; i ++ ) {
        cin >> a[i];
    }
    vector<vector<array<int, 3>>> dp(n + 1, vector<array<int, 3>>(k + 1, array<int, 3>({-inf, -inf, -inf})));
    for(int j = 0; j <= k; j ++ ) {
        for(int t = 0; t < 3; t ++ ) {
            dp[1][j][t] = 0;
        }
    }

    for(int i = 1; i <= n; i ++ ) {
        for(int j = 0; j <= k; j ++ ) {
            for(int t = 0; t < 3; t ++ ) {
                dp[i][j][t] = max(dp[i][j][t], dp[i - 1][j][t] + (t == 1 ? a[i] - x : 0));
                if(j < k) {
                    dp[i][j + 1][t] = max(dp[i][j + 1][t], dp[i - 1][j][t] + (t == 1 ? a[i] + x : 0));
                }
            }
        }
        for(int j = 0; j <= k; j ++ ) {
            for(int t = 1; t < 3; t ++ ) {
                dp[i][j][t] = max(dp[i][j][t], dp[i][j][t - 1]);
            }
        }
    }
    cout << dp[n][k][2] << endl;
}

signed main() {
    ios::sync_with_stdio(0);
    cin.tie(0);

    int Case;
    cin >> Case;
    while(Case -- ) solve();

    return 0;
}