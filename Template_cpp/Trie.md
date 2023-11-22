## Trie

```c++
template<char FIRST = 'a'>
struct Trie {
    static constexpr int ALPHABET_SIZE = 26;
    std::vector<std::array<int, ALPHABET_SIZE>> trie;
    std::vector<int> cnt;
    int cur;

    Trie() {}
    Trie(int N_) {
        cur = 0;
        trie.assign(N_ * ALPHABET_SIZE, {});
        cnt.assign(N_ * ALPHABET_SIZE, 0);
    }

    // 插入S串
    void insert(std::string S) {
        int p = 0;
        for (auto c : S) {
            int u = c - FIRST;
            if (!trie[p][u]) {
                trie[p][u] = ++cur;
            }
            p = trie[p][u];
        }
        cnt[p]++;
    }

    // 查询S串出现的次数
    int queryFrequency(std::string S) {
        int p = 0;
        for (auto c : S) {
            int u = c - FIRST;
            if (!trie[p][u]) {
                return 0;
            }
            p = trie[p][u];
        }
        return p;
    }
};
```

### 说明
1. 传入的参数为字符串最长的长度
2. 如果字符串为大写字母, 可以传入'A'
```c++
Trie trie(MAX_LENGTH)
Trie<'A'> trie(MAX_LENGTH)
```