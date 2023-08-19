```c++
int qmi(int a, int k, int p){  // 求a^k mod p
    int res = 1 % p;
    while (k){
        if (k & 1) res = (LL)res * a % p;
        a = (LL)a * a % p;
        k >>= 1;
    }
    return res;
}


int exgcd(int a, int b, int &x, int &y){	// 扩展欧几里得算法, 求x, y，使得ax + by = gcd(a, b)
    if(!b) return x = 1, y = 0, a;
    int r = exgcd(b, a % b, x, y);
    tie(x, y) = make_tuple(y, x - (a / b) * y);
    return r;
}

int bsgs(int a, int b, int p){      //当a和p互素时
    if (1 % p == b % p) return 0;
    int k = sqrt(p) + 1;
    unordered_map<int, int> hash;
    for (int i = 0, j = b % p; i < k; i ++ ){
        hash[j] = i;
        j = (LL)j * a % p;
    }
    int ak = qmi(a, k, p);
    for (int i = 1, j = ak; i <= k; i ++ ){
        if (hash.count(j)) return i * k - hash[j];
        j = (LL)j * ak % p;
    }
    return -INF;
}

int exbsgs(int a, int b, int p){    //当a和p不互素时
    b = (b % p + p) % p;
    if (1 % p == b % p) return 0;
    int x, y;
    int d = exgcd(a, p, x, y);
    if (d > 1){
        if (b % d) return -INF;
        exgcd(a / d, p / d, x, y);
        return exbsgs(a, (LL)b / d * x % (p / d), p / d) + 1;
    }
    return bsgs(a, b, p);
}
```