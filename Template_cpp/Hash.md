```c++
using u64 = unsigned long long;

constexpr u64 mod = (1ull << 61) - 1;

std::mt19937_64 rnd(std::random_device{}());
std::uniform_int_distribution<u64> dist(100, mod - 1);
const u64 base = dist(rnd);

u64 power[N];

static inline u64 add(u64 a, u64 b){
    a += b;
    if (a >= mod) a -= mod;
    return a;
}

static inline u64 mul(u64 a, u64 b){
    __uint128_t c = __uint128_t(a) * b;
    return add(c >> 61, c & mod);
}

u64 merge(u64 h1, u64 h2, int len2){
    return add(mul(h1, power[len2]), h2);
}

void init(){
    power[0] = 1;
    for(int i = 1; i < N; i++)
        power[i] = mul(power[i - 1], base);
}

u64 query(const std::vector<u64> &s, int l, int r){
    return add(s[r], mod - mul(s[l - 1], power[r - l + 1]));
}

std::vector<u64> build(const std::string &s){
    int sz = s.size();
    std::vector<u64> hashed(sz + 1);
    for (int i = 0; i < sz; i++){
        hashed[i + 1] = add(mul(hashed[i], base), s[i]);
    }
    return hashed;
}

template <typename T>
std::vector<u64> build(const std::vector<T> &s){
    int sz = s.size();
    std::vector<u64> hashed(sz + 1);
    for (int i = 0; i < sz; i++){
        hashed[i + 1] = add(mul(hashed[i], base), s[i]);
    }
    return hashed;
}

int lcp(const std::vector<u64> &a, int l1, int r1, const std::vector<u64> &b, int l2, int r2){
    int len = std::min(r1 - l1 + 1, r2 - l2 + 1);
    int l = 0, r = len;
    while(l < r){
        int mid = (l + r + 1) / 2;
        if (query(a, l1, l1 + mid - 1) == query(b, l2, l2 + mid - 1)) l = mid;
        else r = mid - 1;
    }
    return r;
}
```