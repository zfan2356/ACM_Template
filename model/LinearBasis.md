```c++
template<typename T, int MAXL>
struct LinearBasis {
    int rank;
    std::vector<T> a;

    LinearBasis() {
        a.resize(MAXL);
        rank = 0;
        std::fill(a.begin(), a.end(), T());
    }

    LinearBasis(std::vector<T>& _init) {
        a.resize(MAXL);
        rank = 0;
        std::fill(a.begin(), a.end(), T());
        build(_init);
    }

    void insert(T t) {
        if (t == 0) {
            return;
        }
        for (int j = MAXL - 1; j >= 0; j--){
            if (!(t & (1ll << j))) {
                continue;
            }

            if (a[j]) {
                t ^= a[j];
            }
            else {
                for (int k = 0; k < j; k++) {
                    if (t & (1ll << k)) {
                        t ^= a[k];
                    }
                }
                for (int k = j + 1; k < MAXL; k++) {
                    if (a[k] & (1ll << j)) {
                        a[k] ^= t;
                    }
                }
                a[j] = t;
                rank++;
                break;
            }
        }
    }

    void build(std::vector<T>& _init){
        a.resize(MAXL);
        std::fill(a.begin(), a.end(), T());
        for (auto x : _init) {
            insert(x);
        }
    }

    T queryMax(){
        T res = 0;
        for (int i = 0; i < MAXL; i++) {
            res ^= a[i];
        }
        return res;
    }

    T queryKth(int k) { // 如果不满秩, k--;
        T ans = 0;
        if (k >= (1ll << (int)a.size())) {
            return -1;
        }
        for (int i = MAXL - 1; i >= 0; i--) {
            if (k & (1ll << i)) {
                ans ^= a[i];
            }
        }
        return ans;
    }

    static LinearBasis merge(const LinearBasis &a, const LinearBasis &b){
        LinearBasis res = a;
        for (int i = 0; i < MAXL; i++) {
            res.insert(b.a[i]);
        }
        return res;
    }

    friend std::ostream &operator<<(std::ostream &out, const LinearBasis<T, MAXL> &L) {
        out << "Rank: " << L.rank << '\n';
        for (auto x : L.a) {
            out << bitset<MAXL>(x) << '\n';
        }
        return out;
    }
};


```