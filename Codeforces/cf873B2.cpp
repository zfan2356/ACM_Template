#include <bits/stdc++.h>

//#define int long long
#define endl '\n'
using namespace std;
using ll = long long;

template<typename A, typename B>
inline ostream &operator<<(ostream &out, const pair <A, B> &p) {
    return out << "(" << p.first << ", " << p.second << ")";
}

template<typename T>
inline ostream &operator<<(ostream &out, const vector <T> &a) {
    out << "[";
    for (int i = 0; i < a.size(); i++) {
        if (i) out << ',';
        out << ' ' << a[i];
    }
    return out << " ]";
}

template<typename T>
inline ostream &operator<<(ostream &out, const set <T> &a) { return out << vector<T>(all(a)); }

const int N = 2e5 + 10, mod = 998244353;

void solve() {
    int n; cin >> n;
    vector<int> a(n), b(n);
    for (int i = 0; i < n; i ++ ) {
        cin >> a[i];
        id[i] = i + 1;
    }
    sort(b.begin(), b.end(), [&](int x, int y) {
        return a[x] < a[y];
    });
    a = b;
    ll ans = 0;
    for (int i = 0; i < n; i ++ ) {
        ans += 1ll * i * (n - i);
    }// all answer

    stack<int> stk{-1};
    vector<int> L(n), R(n);
    for (int i = 0; i < n; i ++ ) {
        while (stk.size() > 1 && a[stk.top()] > a[i]) stk.pop();
        L[i] = stk.top();
        stk.push(i);
    }
    stk = {n};
    for (int i = n - 1; i >= 0; i -- ) {
        while (stk.size() > 1 && a[stk.top()] > a[i]) stk.pop();
        R[i] = stk.top();
        stk.push(i);
    }

    set<int> s;
    s.insert(-1);
    for (int i = n - 1; i >= 0; i -- ) {
        int pos = a[i];
        int l = L[pos], r = R[pos];
        int k = l == -1 ?
    }

    cout << ans << endl;
}

signed main() {
    ios::sync_with_stdio(0);
    cin.tie(0);

#ifdef DEBUG
    freopen("D:\\mypile\\acmC++\\in.txt", "r", stdin);
    freopen("D:\\mypile\\acmC++\\out.txt", "w", stdout);
#endif

    int Case;
    cin >> Case;
    while (Case--) solve();

    return 0;
}