## Trie

```c++
struct Trie {
    static constexpr int ALPHABET_SIZE = 26;
    static constexpr char FIRST = 'a';
    std::vector<std::array<int, ALPHABET_SIZE>> trie;
    std::vector<int> cnt;
    int cur;

    Trie() {}
    Trie(int N_) {
        cur = 0;
        trie.assign(N_, {});
        cnt.assign(N_, 0);
    }

    // 插入S串
    void insert(std::string S) {
        int p = 0;
        for (auto c : S) {
            int &q = trie[p][c - FIRST];
            if (!q) {
                q = ++cur;
            }
            p = q;
        }
        cnt[p]++;
    }

    // 查询S串出现的次数
    int queryFrequency(std::string S) {
        int p = 0;
        for (auto c : S) {
            int &q = trie[p][c - FIRST];
            if (!q) {
                return 0;
            }
            p = q;
        }
        return cnt[p];
    }
};
```

### 说明
1. 传入的参数为字符串最长的长度
2. 如果字符串为大写字母, 别忘记修改`FIRST`常量
```c++
Trie trie(MAX_LENGTH)
```
3. 实际较为灵活, 模板起不到太大的作用