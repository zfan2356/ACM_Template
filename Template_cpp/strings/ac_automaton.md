### AC自动机
> 解决多模式串的匹配问题

例题: 给定若干模式串$t_i$ 以及 模板串$s$, 问$s$中每个$t_i$出现了多少次
```c++
constexpr int N = 5E5 + 10;

int trie[N][26];
int fail[N];

int tot = 1;

void solve() {
    std::string s;
    std::cin >> s;

    int q;
    std::cin >> q;

    std::vector<std::string> t(q);
    for (int i = 0; i < q; i++) {
        std::cin >> t[i];
    }

    for (int i = 0; i < 26; i++) {
        trie[0][i] = 1;
    }

    std::vector<int> end(q);
    for (int i = 0; i < q; i++) {
        int p = 1;
        for (auto c : t[i]) {
            int &q = trie[p][c - 'a'];
            if (q == 0) {
                q = ++tot;
            }
            p = q;
        }
        end[i] = p;
    }

    std::queue<int> qu;
    qu.push(1);

    std::vector<int> nodes;
    while (!qu.empty()) {
        int x = qu.front();
        nodes.push_back(x);
        qu.pop();

        for (int i = 0; i < 26; i++) {
            if (trie[x][i] == 0) {
                trie[x][i] = trie[fail[x]][i];
            } else {
                fail[trie[x][i]] = trie[fail[x]][i];
                qu.push(trie[x][i]);
            }
        }
    }

    int p = 1;
    std::vector<int> f(tot + 1);
    for (auto c : s) {
        p = trie[p][c - 'a'];
        f[p]++;
    }

    for (int i = nodes.size() - 1; i >= 0; i--) {
        int x = nodes[i];
        f[fail[x]] += f[x];
    }

    for (int i = 0; i < q; i++) {
        std::cout << f[end[i]] << '\n';
    }

}
```