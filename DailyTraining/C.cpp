#include <bits/stdc++.h>

//#define int long long
#define endl '\n'
using namespace std;
using LL = long long;

template<typename A, typename B> inline ostream& operator<< (ostream& out, const pair<A, B>& p) { return out << "(" << p.first << ", " << p.second << ")"; }
template<typename T> inline ostream& operator<< (ostream& out, const vector<T>& a) { out << "["; for(int i = 0; i < a.size(); i ++) { if (i) out << ','; out << ' ' << a[i]; } return out << " ]"; }
template<typename T> inline ostream& operator<< (ostream& out, const set<T>& a) { return out << vector<T>(all(a)); }

const int mod = 998244353;


signed main() {
    ios::sync_with_stdio(0);
    cin.tie(0);

    int n, k; cin >> n >> k;
    string s; cin >> s;
    vector<pair<int, int>> seg;
    int x = 0;
    for(int i = 0; i < n; i ++ ) {
        x += (s[i] == 'X');
    }
    if (k > x) {
        k = n - k, x = n - x;
        for (int i = 0; i < n; i++) s[i] ^= 'X' ^ 'Y';
    }
    if(s == string(n, 'X')) {
        cout << max(k - 1, 0) << endl;
    } else {
        for (int i = 0; i < n; i++) {
            if (s[i] == 'X') {
                int j = i;
                while (j + 1 < n && s[j + 1] == 'X') j ++;
                if (j + 1 < n && s[j + 1] == 'Y' && i && s[i - 1] == 'Y') seg.push_back({i, j});
                i = j;
            }
        }
        sort(seg.begin(), seg.end(), [&](auto x, auto y) {
            return x.second - x.first + 1 < y.second - y.first + 1;
        });
        int ans = 0;
        int t, K;
        string S(s);
        for (t = 0, K = k; t < seg.size(); t++) {
            int v = seg[t].second - seg[t].first + 1;
            if (K - v < 0) break;
            K -= v;
            for (int i = seg[t].first; i <= seg[t].second; i++) S[i] = 'Y';
        }
        for (int i = 1; i < n; i++) {
            if (S[i] == S[i - 1] && S[i] == 'Y') ans++;
        }
        cout << ans + K << endl;
    }

    return 0;
}