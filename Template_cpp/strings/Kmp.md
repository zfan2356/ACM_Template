```c++
struct KMP {
    // 计算前缀函数, 其中S是以0开始的一个string
    std::vector<int> prefixFunction(std::string S) {
        int n = S.size();
        std::vector<int> P(n);
        for (int i = 1; i < n; i++) {
            int j = P[i - 1];
            while (j && S[i] != S[j]) {
                j = P[j - 1];
            }
            j += (S[i] == S[j]);
            P[i] = j;
        }
        return P;
    }

    // 计算Pattern序列在Text序列中的匹配次数, 返回位置的集合
    std::vector<int> work(std::string Text, std::string Pattern) {
        std::vector<int> Pos;
        std::string Splicing = Pattern + "#" + Text;
        int n = Text.size(), m = Pattern.size();
        std::vector<int> prefix = prefixFunction(Splicing);
        for (int i = m + 1; i <= n + m; i++) {
            if (prefix[i] == m) {
                Pos.push_back(i - 2 * m);
            }
        }
        return Pos;
    }
};
```


```c++
for (int i = 0, j = 0; i < m; i++) {
    while (j && s[i] != t[j]) {
         j = ne[j - 1];
        }
        j += s[i] == t[j];
        if (j == n) {
            cout << i - n + 1 << ' ';
        }
    }
```
