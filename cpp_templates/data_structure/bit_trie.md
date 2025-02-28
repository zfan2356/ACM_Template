## 01Trie

```c++
template<typename T>
struct BitTrie {
    std::vector<std::array<int, 2>> trie;
    std::vector<int> cnt;
    int cur, MAX_BIT;

    BitTrie() {}
    BitTrie(int N_) {
        cur = 0;
        MAX_BIT = std::__lg(std::numeric_limits<T>::max()) + 1;
        trie.assign(N_ * MAX_BIT, {});
        cnt.assign(N_ * MAX_BIT, 0);
    }

    void insert(T x, T v) {
        for (int i = MAX_BIT - 1, p = 0; i >= 0; i--) {
            int u = x >> i & 1;
            if (!trie[p][u]) {
                trie[p][u] = ++cur;
            }
            p = trie[p][u];
            cnt[p] += v;
        }
    }

    T queryMaxXOR(T x) {
        T res = 0;
        for (int i = MAX_BIT - 1, p = 0; i >= 0; i--) {
            int u = x >> i & 1;
            if (trie[p][u ^ 1] && cnt[trie[p][u ^ 1]]) {
                res = res * 2 + 1;
                p = trie[p][u ^ 1];
            } else {
                p = trie[p][u];
                res *= 2;
            }
        }
        return res;
    }
};


```

### 说明
初始化传入插入的规模
```c++
BitTrie<int> trie(n); // 插入n个数
```
