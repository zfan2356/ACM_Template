#### 求逆元

```C++
vector<LL> inv(n + 1);
inv[1] = 1;
for(int i = 2; i <= n; i ++ ) inv[i] = mod - inv[mod % i] * (mod / i) % mod;
```

#### 		欧拉函数

```C++
bool st[N];
int primes[N], phi[N], cnt;
void init(int n){
    phi[1] = 1;
    for(int i = 2; i <= n; i ++ ){
        if(!st[i]){
            primes[cnt ++ ] = i;
            phi[i] = i - 1;
        }
        for(int j = 0; primes[j] * i <= n; j ++ ){
            st[primes[j] * i] = 1;
            if(i % primes[j] == 0){
                phi[i * primes[j]] = phi[i] * primes[j];
                break;
            }
            phi[i * primes[j]] = phi[i] * (primes[j] - 1);
        }
    }
}
```

小应用: 1到$n$ 中有多少对数, 满足两两互质? 

答案: $\phi(n) * 2 + 1$

#####     			2. 扩展欧拉定理

不要求$(a, m) = 1$ , 当$b \ge \phi(m)$ 时

$a^{b} \equiv a ^ {(b\ mod\ \phi(m)) + \phi(m)}\ mod\ m$ 

应用: **欧拉降幂**

处理$\phi(m),$ $\phi(\phi(m))....$  这种处理会使得$m$ 以$log$ 的速度下降到1

```c++
int n, m, w[N], q;
int phi[N], tot;
bool ok;

int get_euler(int n){		//计算单个的欧拉函数
    int res = n;
    for(int i = 2; i <= n / i; i ++ ){
        if(n % i == 0){
            res = res / i * (i - 1);
            while(n % i == 0) n /= i;
        }
    }
    if(n > 1) res = res / n * (n - 1);
    return res;
}

void init(){		//初始化phi的值
    int tmp = m;
    tot = 0;
    phi[++ tot] = tmp;
    while(tmp != 1){
        tmp = get_euler(tmp);
        phi[++ tot] = tmp;
    }
}

int qmi(int a, int k, int p){	//快速幂, 同时判断一下是否大于phi(m)
    int res = 1 % p;
    while (k){
        if (k & 1) {
            res = (LL) res * a;
            if(res >= p){
                ok = 1;
                res %= p;
            }
        }
        a = (LL)a * a;
        if(a >= p){
            ok = 1;
            a %= p;
        }
        k >>= 1;
    }
    return res;
}

int Euler_Drop(int l, int r){   //求序列l到r的幂次 mod m的值
    r = min(r, l - 1 + tot);
    int res = 1;
    for(int i = r; i >= l; i -- ){
        ok = 0;
        res = qmi(w[i], res, phi[i - l + 1]);
        if(ok) res += phi[i - l + 1];
    }
    return res % phi[1];
}
```



#### 		积性函数

![image-20220819182846898](C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20220819182846898.png)

![image-20220819182906130](C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20220819182906130.png)

![image-20220819182951598](C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20220819182951598.png)

![image-20220819183038961](C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20220819183038961.png)

例题: 求解<img src="C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20221103234046465.png" alt="image-20221103234046465" style="zoom: 80%;" />

推导如下: 

<img src="C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20221103234130677.png" alt="image-20221103234130677" style="zoom: 80%;" />

法一:($O(nlogn)$ 预处理) 

法二: 利用积性函数性质, 因为$d$ 为积性函数, $\phi$ 也为积性函数, 所以两者相乘的和也是积性函数, 所以我们设$f(i) = \sum_{d|i}{d\phi(\frac{i}{d})}$ 

那么直接筛$f$ 即可, 注意筛积性函数的几个特判, 要注意$f(p^k)$ , 这个需要手动推导一下:

$f(p^{k+1}) = \sum_{d | p^{k+1}}{d\phi(p^{k+1}/d)} = \sum_{d*p|p^{k+1}}{d*p\phi(p^k/d)} + \phi(p^{k+1}) = p * f(p^k) + p^{k}(p - 1)$ 

这样的话我们就可以筛了

```c++
void init(int n){
    f[1] = 1;
    low[1] = 1;
    for(int i = 2; i <= n; i ++ ){
        if(!st[i]){
            primes[cnt ++ ] = i;
            f[i] = i + i - 1;
            low[i] = i;
        }
        for(int j = 0; primes[j] * i <= n; j ++ ){
            st[primes[j] * i] = 1;
            if(i % primes[j] == 0){
                //此时i中最小质因子仍然是primes[j] 但是其中有
                if(low[i] == i) f[i * primes[j]] = f[i] * primes[j] + i * (primes[j] - 1);
                else f[i * primes[j]] = f[i / low[i]] * f[primes[j] * low[i]];

                low[i * primes[j]] = primes[j] * low[i];
                break;
            }
            f[i * primes[j]] = f[i] * f[primes[j]];
            // 此时i * primes[j]中, primes[j]为最小质因子, 且i中没有, 所以i与primes[j]互质, 就可以直接根据积性函数算
            low[i * primes[j]] = primes[j];
            //最小质因子自然就是primes[j]
        }
    }
}
```



#### 二次剩余

+ 定义：

如果存在一个整数$x$ , 使得 $x ^ 2 \equiv n\ (mod\ p)$ 那么就称$n$ 是模 $p$ 的二次剩余

1. 判断方法

   勒让德符号$(\frac{n}{p})$ , 如果$n$ 是模 $p$ 的二次剩余，那么$(\frac{n}{p})$ = 1；如果不是，那么$(\frac{n}{p})$ = -1；如果$p | n$ ，那么$(\frac{n}{p})$ = 0;

   其中$(\frac{n}{p}) = n^{\frac{p - 1}{2}}$ ， 其中 $p$ 为奇质数

2. 定理：

   对于方程 $x ^ 2 \equiv n\ (mod\ p)$ , 有 $\frac{p - 1}{2}$ 个 $n$ , 使得该方程有解
   
3.  $\sqrt{2}$ 在 mod $1e9 + 7$意义下的值为$59713600$ 



#### 		杜教筛

![image-20220819183131140](C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20220819183131140.png)

![image-20220819183201281](C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20220819183201281.png)



![image-20220911172751384](C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20220911172751384.png)

```c++
namespace DuJS{
    const int N = 1e7;

    unordered_map<LL, LL> Sphi;
    unordered_map<LL, int> Smu;
    int prime[N], mu[N];
    LL phi[N];
    bool check[N];

    void init() {		//先预处理
        int tot = 0;
        phi[1] = mu[1] = 1;
        for (int i = 2; i < N; i++) {
            if (!check[i]) {
                prime[++tot] = i;
                mu[i] = -1;
                phi[i] = i - 1;
            }
            for (int j = 1; j <= tot; j++) {
                if (i * prime[j] >= N) break;
                check[i * prime[j]] = true;
                if (i % prime[j] == 0) {
                    mu[i * prime[j]] = 0;
                    phi[i * prime[j]] = phi[i] * prime[j];
                    break;
                }
                mu[i * prime[j]] = -mu[i];
                phi[i * prime[j]] = phi[i] * (prime[j] - 1);
            }
        }
        for (int i = 1; i < N; i++)
            mu[i] += mu[i - 1], phi[i] += phi[i - 1];
    }

    int cal_Smu(LL n) {				//求莫比乌斯函数前缀和
        if (n < N) return mu[n];
        if (Smu[n]) return Smu[n];
        int ret = 1;
        for (LL i = 2, last; i <= n; i = last + 1) {
            last = n / (n / i);
            ret -= cal_Smu(n / i) * (last - i + 1);
        }
        return Smu[n] = ret;
    }

    LL cal_Sphi(LL n) {				//求欧拉函数前缀和
        if (n < N) return phi[n];
        if (Sphi[n]) return Sphi[n];
        LL ret = 1LL * (1 + n) * n / 2;
        for (LL i = 2, last; i <= n; i = last + 1) {
            last = n / (n / i);
            ret -= cal_Sphi(n / i) * (last - i + 1);
        }
        return Sphi[n] = ret;
    }
}
```



#### 		Min_25筛

![image-20221006210527002](C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20221006210527002.png)

![image-20221006210503955](C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20221006210503955.png)

```c++
namespace min_25 { // 假设 f(p)=p^2-p
    const int M = 1e6 + 5, inv2 = 500000004, inv3 = 333333336;
    bool isnp[M] = {1, 1};
    int p[M], num;
    int id1[M], id2[M];
    LL B, N;
    LL sp1[M], sp2[M], w[M];
    LL g1[M], g2[M];    //求所有质数的积性函数的和

    void getPrime(int n) {      //筛出1 - n内的所有的质数
        num = 0;
        for (int i = 2; i <= n; ++ i) {
            if (!isnp[i]) p[++num] = i;
            for (int j = 1; j <= num && i * p[j] <= n; j++) {
                isnp[i * p[j]] = 1;
                if (i % p[j] == 0) break;
            }
        }
        sp1[0] = 0; sp2[0] = 0;
        //根据题目要求的积性函数来处理求和, 这里f1(p) = p^2, f2(p) = p, f(p) = f1(p) - f2(p)
        for (int i = 1; i <= num; ++i) sp1[i] = (sp1[i - 1] + p[i]) % mod;
        for (int i = 1; i <= num; ++i) sp2[i] = (sp2[i - 1] + (LL) p[i] * p[i] % mod) % mod;
    }

    LL calc_g1(LL x) {  //求和1 - n的i
        x %= mod;
        return x * (x + 1) % mod * inv2 % mod;
    }

    LL calc_g2(LL x) {  //求和1 - n的i^2
        x %= mod;
        return x * (x + 1) % mod * (2 * x + 1) % mod * inv2 % mod * inv3 % mod;
    }

    LL fpk(LL x) {      //求解f(x)
        x %= mod;
        return x * (x - 1) % mod;
    }

    int getId(LL x) {       //获得值对应的id
        return x <= B ? id1[x] : id2[N / x];
    }

    LL go(LL x, int y) {        //求解h(x, y) 2 - x中, min_p(y) >= p[y];
        if (p[y] >= x) return 0;
        int k = getId(x);
        LL ans = (g2[k] - g1[k] - sp2[y] + sp1[y] + mod + mod) % mod;
        for (int i = y + 1; i <= num && (LL) p[i] * p[i] <= x; ++i) {
            LL pe = p[i];
            for (int e = 1; pe <= x; ++e, pe *= p[i])
                ans = (ans + fpk(pe) * (go(x / pe, i) + (e != 1)) % mod) % mod;
        }
        return ans;
    }

    LL solve(LL _N) {
        N = _N;
        B = sqrt(N + 0.5);
        getPrime(B);    //处理根号n下的质数
        // 求解g(n/x,0) 并记录 n/x 对应的下标
        int tot = 0;
        for (LL l = 1, v, r; l <= N; l = r + 1) {   //分块处理
            v = N / l;
            r = N / v;
            w[++tot] = v;   //记录每一块的权值
            g1[tot] = calc_g1(v) - 1;       //记录g(v, 0), 因为v的范围过大, 其实相当于离散化
            g2[tot] = calc_g2(v) - 1;
            if (v <= B) id1[v] = tot;
            else id2[N / v] = tot;  //因为v过大, 所以开两个数组存储
        }
        for (int i = 1; i <= num; ++i)
            for (int j = 1; j <= tot && (LL) p[i] * p[i] <= w[j]; j ++) {
                int k = getId(w[j] / p[i]);
                g1[j] = (g1[j] - (LL) p[i] * (g1[k] - sp1[i - 1] + mod) % mod + mod) % mod;
                g2[j] = (g2[j] - (LL) p[i] * p[i] % mod * (g2[k] - sp2[i - 1] + mod) % mod + mod) % mod;
            }
        return (go(N, 0) + 1) % mod;
    }
}

```



####  		数论分块

$$
\sum_{i=1}^{m}{[\frac{n}{i}]}的结果是多少？
$$

```c++
LL ans = 0;
LL k = n;
for(int l = 1, r; l <= min(n, m); l = r + 1){
    r = min(n / (n / l), min(n, m));
    ans += (k / l) * (r - l + 1);
}
cout << ans << endl;

```



#### 		扩展中国剩余定理

```c++
int mul(int a, int b, int mod) {	//龟速乘, 有必要可以加
    int ans = 0;
    for (; b; a = (a + a) % mod, b >>= 1)
        if (b & 1)
            ans = (ans + a) % mod;
    return ans;
}

int _mod(int a, int b) {
    return ((a % b) + b) % b;
}

int exgcd(int a, int b, int &x, int &y){	// 扩展欧几里得算法, 求x, y，使得ax + by = gcd(a, b)
    if(!b) return x = 1, y = 0, a;
    int r = exgcd(b, a % b, x, y);
    tie(x, y) = make_tuple(y, x - (a / b) * y);
    return r;
}

int exCRT(int n, int a[], int m[]){  // x = ai (mod mi)   i下标[0, n)
    for (int i = 1; i < n; i ++) {
        int k1, k2;
        int d = exgcd(m[0], m[i], k1, k2);
        int c = (a[i] % m[i] - a[0] % m[i] + m[i]) % m[i];
        if (c % d) return -1;
 
        k1 = _mod(k1 * c / d, abs(m[i] / d));
        //这里如果爆ll可以换龟速乘
        a[0] = k1 * m[0] + a[0];
        m[0] = abs(m[0] / d * m[i]);
        a[0] = _mod(a[0], m[0]);
    }
 
    return a[0];
}

```



#### 		莫比乌斯反演

![image-20220819182927022](C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20220819182927022.png)

莫比乌斯函数性质: 

​			Ⅰ: 

![img](https://bkimg.cdn.bcebos.com/formula/905e5d424a8756c5acef64d916b2cdb2.svg)

​			Ⅱ: 

![img](https://bkimg.cdn.bcebos.com/formula/18a1d7ba4d6c435888a5bd5547b95565.svg)

1. 莫比乌斯函数预处理法(同时求前缀和)

筛法, 时间复杂度$O(n)$

```C++
int primes[N], _cnt, mu[N], sum[N];
bool st[N];

void init(){
    mu[1] = 1;
    for (int i = 2; i < N; i ++ ){
        if (!st[i]) primes[_cnt ++ ] = i, mu[i] = -1;
        for (int j = 0; primes[j] * i < N; j ++ ){
            st[primes[j] * i] = true;
            if (i % primes[j] == 0) break;
            mu[primes[j] * i] = -mu[i];
        }
    }
    for (int i = 1; i < N; i ++ ) sum[i] = sum[i - 1] + mu[i];
}

```



#### 		BSGS算法

求$a^{x} ≡ b (mod\ p)$ 的最小非负整数$x$ 

**朴素版** $(a, p) = 1$ 

```C++
int qmi(int a, int k, int p){  // 求a^k mod p
    int res = 1 % p;
    while (k){
        if (k & 1) res = (LL)res * a % p;
        a = (LL)a * a % p;
        k >>= 1;
    }
    return res;
}

int bsgs(int a, int b, int p){	//求a^x = b (mod p)
    if (1 % p == b % p) return 0;	//注释掉这句话代表的是求最小正数解x
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

```



**扩展BSGS** (当$(a, p) != 1$时)

```C++
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

int bsgs(int a, int b, int p){
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

int exbsgs(int a, int b, int p){
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



#### 		Miller_Rabin素性检验

```c++
namespace Miller_Rabin{
//    typedef __int128 LL;		//如果范围较大可以使用__int128

    __int128 read() {
        __int128 x = 0, f = 1;
        char ch = getchar();
        while (!isdigit(ch) && ch != '-')ch = getchar();
        if (ch == '-')f = -1, ch = getchar();
        while (isdigit(ch))x = x * 10 + ch - '0', ch = getchar();
        return f * x;
    }

    void print(__int128 x) {
        if (x < 0) {
            putchar('-');
            x = -x;
        }
        if (x > 9)
            print(x / 10);
        putchar(x % 10 + '0');
    }
    LL gcd(LL a, LL b) {
        return b ? gcd(b, a % b) : a;
    }
    LL mul(LL a, LL b, LL mod) {		//龟速乘
        LL ret = 0;
        while (b) {
            if (b & 1) ret = (ret + a) % mod;
            a = (a + a) % mod;
            b >>= 1;
        }
        return ret;
    }
    LL pow(LL a, LL b, LL mod) {		//快速幂
        LL ret = 1;
        while (b) {
            if (b & 1) ret = mul(ret, a, mod);
            a = mul(a, a, mod);
            b >>= 1;
        }
        return ret;
    }
    bool check(LL a, LL n) {
        LL x = n - 1, y;
        int t = 0;
        while ((x & 1) == 0) {
            x >>= 1;
            t++;
        }
        x = pow(a, x, n);
        for (int i = 1; i <= t; i++) {
            y = mul(x, x, n);
            if (y == 1 && x != 1 && x != n - 1) return true;
            x = y;
        }
        if (y != 1) return true;
        return false;
    }
    bool isprime(LL n) {
        if (n == 2) return true;
        if (n == 1 || !(n & 1)) return false;
        const int arr[12] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
        for (int i = 0; i < 12; i++) {
            if (arr[i] >= n) break;
            if (check(arr[i], n)) return false;
        }
        return true;
    }
}

```



 #### 		Pollard_rho

```c++
mt19937 mt(time(0)); //随机化
inline LL PR(LL n) {
    LL x = uniform_int_distribution<LL>(0, n - 1)(mt), s, t, c = uniform_int_distribution<LL>(1, n - 1)(mt); //随机化
    for (int gol = 1; 1; gol <<= 1, s = t, x = 1) {
        for (int stp = 1; stp <= gol; ++stp) {
            t = (Miller_Rabin::mul(t, t, n) + c) % n;
            x = Miller_Rabin::mul(x, abs(s - t), n);
            if ((stp & 127) == 0) {
                LL d = gcd(x, n);
                if (d > 1) return d;
            }
        }
        LL d = gcd(x, n);
        if (d > 1) return d;
    }
}

```





### 	组合数学

#### 		组合数性质

![image-20220820112939142](C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20220820112939142.png)

#### 		求组合数四种情况

##### 			1.递推

**时间复杂度为$O(n^2)$.**

```c++
int c[N][N];
void init(){
    for(int i = 0; i < N; i ++)
        for(int j = 0; j <= i; j ++)
            if(!j) c[i][j] = 1;
            else c[i][j] = (c[i - 1][j] + c[i - 1][j - 1]) % mod;
}
```

##### 			2. 预处理法

**N在1e5附近，时间复杂度为$O(N)$。**

```c++
namespace Combin {
    const int N = 1e6 + 10, mod = 1e9 + 7;
    using ll = long long;
    ll infac[N], fac[N];
    int qpow(int a, int b) {
        ll res = 1;
        for (; b; b >>= 1) {
            if (b & 1) res = 1ll * res * a % mod;
            a = 1ll * a * a % mod;
        }
        return res;
    }
    
    void init(){
        fac[0] = 1, infac[0] = 1;
        for (int i = 1; i < N; i++) fac[i] = 1ll * fac[i - 1] * i % mod;
        infac[N - 1] = qpow(fac[N - 1], mod - 2);
        for(int i = N - 2; i; i -- ) infac[i] = 1ll * infac[i + 1] * (i + 1) % mod;
    }
    
    int C(int a, int b) {
        if(a < 0 || b < 0 || a < b) return 0;
        return 1ll * fac[a] * infac[b] % mod * infac[a - b] % mod;
    }
    
    int A(int a, int b) {
        if(a < 0 || b < 0 || a < b) return 0;
        return 1ll * fac[a] * infac[a - b] % mod;
    }
}
```

##### 			3. Lucas定理![image-20211225194056993](C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20211225194056993.png)

**这种情况下N很大，超出int，但是p很小（要mod的数很小）**

 这样就可以将N降到p的数量级，直接计算即可, 注意要求$p$ 为质数

```c++
int qpow(int a, int k){
    int res = 1;
    while (k){
        if (k & 1) res = (LL)res * a % mod;
        a = (LL)a * a % mod;
        k >>= 1;
    }
    return res;
}

int C(int a, int b, int p){
    if (a < b) return 0;

    LL down = 1, up = 1;
    for (int i = a, j = 1; j <= b; i --, j ++ ){
        up = (LL)up * i % p;
        down = (LL)down * j % p;
    }

    return (LL)up * qmi(down, p - 2) % p;
}

int Lucas(int a, int b, int p){
    if (a < p && b < p) return C(a, b);
    return (LL)Lucas(a / p, b / p) * C(a % p, b % p) % p;
}

```



#### 因式分解

![image-20230318203547320](C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20230318203547320.png)



#### 		隔板法

1. $x_1 + x_2 + ... + x_n = k\ \ \ (x_i > 0)$

   答案为$C_{k - 1}^{n - 1}$

2. $x_1 + x_2 + ... + x_n = k\ \ \ (x_i >= 0)$

   答案为$C_{n + k - 1}^{n - 1}$

3. $x_1 + x_2 + ... + x_n <= k\ \ \ (x_i >= 0)$

​		答案为$C_{n + k}^{n}$

(所有部分都加1， 转化为$n + k$ 个物品分为$n$ 个部分，最后一个部分大于等于0，这样的话将最后一个部分舍弃，就可以达到小于号的目的，这样的话有$n + k$ 个空, 插入$n$ 个板子)



#### 		康托展开/逆康托展开

​		Ⅰ: $a$ 为排列, $n$ 为排列长度

```c++
int cantor(int a[], int n) {	//求一个全排列在所有全排列中的次序
    int ans = 0, sum = 0;
    for(int i = 1; i < n; i ++ ) {
        for(int j = i + 1; j <= n; j ++ ) {
            if(a[j] < a[i]) sum ++ ;
        }
        ans += sum * fac[n - i];
        sum = 0;
    }
    return ans + 1;
}
```

​	Ⅱ: $x$ 为排列次序, $n$ 为排列长度

```c++
vector<int> decantor(int x, int n) {	// 求一个给定长度以及排名的排列
 	vector<int> v, a;
    for(int i = 1; i <= n; i ++ ) v.push_back(i);
    for(int i = n; i >= 1; i -- ) {
        int r = x % fac[i - 1], t = x / fac[i - 1];
        x = r;
        sort(v.begin(), v.end());
        a.push_back(v[t]);
        v.erase(v.begin() + t);
    }
    return a;
}
```



#### 		Catalan数

数列前几项: 1, 1, 2, 5, 14, 42, 132, 429, 1430, 4862...

典型应用: 火车进出站, 括号匹配, 满二叉树的数量(n+1个节点可以组成多少满二叉树)

通项: $C_n = \frac{c_{2n}^{n}}{n + 1}$ 

递推公式: $C_1 = 1, C_{n} = C_{n -1} * \frac{4*n - 2}{n + 1}$ 



#### 		Stirling 数

##### 			第一类斯特林数

 **将$n$个不同的元素构成$m$个圆排列的方案数(m个轮换)**

![image-20220907205323710](C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20220907205323710.png)



**快速求一行**

<img src="C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20221021215122459.png" alt="image-20221021215122459" style="zoom:67%;" />

根据生成函数, 我们考虑倍增: <img src="C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20220907123755917.png" alt="image-20220907123755917" style="zoom: 67%;" />

也就是$\prod_{i=0}^{2n-1}{(x + i)}$  = $\prod_{i=0}^{n - 1}{x}$  $\prod_{i=0}^{n-1}{(x + i + n)}$ 

 这相当于我们已知$f(x)$ , 求$f(x + c)$ 

可以推导出来, 如果$f(x) = f_0 + f_1x^1 + ...$ , 则$f(x + c) = \sum_{j=0}^{n}\frac{x^j}{j!}\sum_{i=j}^{n}{i!f_i\frac{c^{i-j}}{(i-j)!}}$ 

证: <img src="C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20221021223105867.png" alt="image-20221021223105867" style="zoom:67%;" />

这就可以转换为卷积的形式, $b(x) = \sum_{i}{\frac{c^i}{i!}x^{i}}$ , $a(x) = \sum_i{i!f_ix^i}$  

反转$a(x)$ 之后卷积, 再反转即可得到

```c++
int infac[N], fac[N];
Poly f;
void calc(int n){
    f.resize(n + 1);
    if(n == 1) {
        f[1] = 1;
        return ;
    }
    else if(n & 1){
        calc(n - 1);
        for(int i = n; i >= 1; i -- ) f[i] = (f[i - 1] + f[i] * (n - 1) % P) % P;
        f[0] = f[0] * (n - 1) % P;
    }
    else{
        int m = n / 2;
        int res = 1;
        calc(m);
        Poly a(n + 1), b(n + 1), c(n + 1);
        for(int i = 0; i <= m; i ++ ) {
            a[i] = f[i] * fac[i] % P;
            b[i] = res * infac[i] % P;
            res = res * m % P;
        }
        reverse(a.begin(), a.begin() + m + 1);
        a = a * b;
        for(int i = 0; i <= m; i ++ ) c[i] = a[m - i] * infac[i] % P;
        f = f * c;
    }
}

```

​	**快速求一列**

<img src="C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20221021223143791.png" alt="image-20221021223143791" style="zoom:67%;" />

当$k = 1$ 时, 因为我们知道$S_1(n, 1) = (n - 1)!$  所以

<img src="C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20220907171833374.png" alt="image-20220907171833374" style="zoom: 67%;" />

那么我们就可以得到对于$k$ 为任意数的指数生成函数为

<img src="C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20220907171913997.png" alt="image-20220907171913997" style="zoom:67%;" />

除以$k!$ 的原因为我们生成的环排列其实有顺序, 所以要除以一个排列, 然后对该式子化简为

<img src="C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20220907172122588.png" alt="image-20220907172122588" style="zoom:67%;" />

直接多项式快速幂计算即可

```c++
int fac[N], infac[N];

void init() {
    fac[0] = infac[0] = 1;
    for(int i = 1; i < N; i ++ ) fac[i] = fac[i - 1] * i % P;
    for(int i = 1; i < N; i ++ ) infac[i] = infac[i - 1] * inv[i] % P;
}

signed main() {
    ios::sync_with_stdio(0);
    cin.tie(0);

    init();
    int n, k; cin >> n >> k;
    Poly F;
    F.resize(n + 1);
    for(int i = 1; i <= n; i ++ ) F[i] = inv[i];
    F = Pow(F, k, k);	//加强版快速幂
    for(int i = 0; i <= n; i ++ ) {
        if(i < k) {
            cout << 0 << " ";
            continue;
        }
        int ans = fac[i] * F[i] % P * infac[k] % P;
        cout << ans << " ";
    }
    return 0;
}
```



##### 		第二类斯特林数: 

**将$n$ 个不同的元素划分为$k$个相同的集合, 保证集合非空的划分方案数**

![image-20220906101810702](C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20220906101810702.png)

**快速求一行:**

<img src="C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20221021230940076.png" alt="image-20221021230940076" style="zoom:67%;" />

 可以将通项公式转换为卷积的形式

原式 = $\frac{1}{k!}\sum_{i = 0}^{k}{(-1)^iC_{k}^{i}(k-i)^n}$  = $\sum_{i = 0}^{k}(-1)^i \frac{1}{i!} * \frac{1}{(k - i)!}(k - i)^n$ 我们转换成两个多项式的卷积的形式: $f(x) = \sum_{i>=0}{(-1)^i\frac{1}{i!}x^i}$  $g(x) = \sum_{j>=0}{\frac{1}{j!}j^nx^j}$ , 第$k$ 项就是答案

```c++
signed main() {
    ios::sync_with_stdio(0);
    cin.tie(0);

    int n; cin >> n;
    fac[0] = infac[0] = 1;
    for(int i = 1; i <= n; i ++ ) fac[i] = fac[i - 1] * i % P;
    for(int i = 1; i <= n; i ++ ) infac[i] = infac[i - 1] * inv[i] % P;

    Poly f(n + 1), g(n + 1);
    for(int i = 0; i <= n; i ++ ) {
        f[i] = infac[i] * (i & 1 ? P - 1 : 1) % P;
        g[i] = infac[i] * qpow(i, n, P) % P;
    }
    f = f * g;
    for(int i = 0; i <= n; i ++ ) cout << f[i] << " ";
    cout << endl;
    return 0;
}
```



**快速求一列:**

```c++
int fac[N], infac[N];

void init() {
    fac[0] = infac[0] = 1;
    for(int i = 1; i < N; i ++ ) fac[i] = fac[i - 1] * i % P;
    for(int i = 1; i < N; i ++ ) infac[i] = infac[i - 1] * inv[i] % P;
}

signed main() {
    ios::sync_with_stdio(0);
    cin.tie(0);

    init();
    int n, k; cin >> n >> k;
    Poly F;
    F.resize(n + 1);
    for(int i = 1; i <= n; i ++ ) F[i] = infac[i];
    F = Pow(F, k, k);
    for(int i = 0; i <= n; i ++ ) {
        if(i < k) {
            cout << 0 << " ";
            continue;
        }
        int ans = fac[i] * F[i] % P * infac[k] % P;
        cout << ans << " ";
    }
    return 0;
}
```



![image-20220907210154206](C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20220907210154206.png)

### 	线性代数

#### 		高斯消元

```c++
int gauss(){  // 高斯消元，答案存于a[i][n]中，0 <= i < n
    int c, r;
    for (c = 0, r = 0; c < n; c ++ ){
        int t = r;
        for (int i = r; i < n; i ++ )  // 找绝对值最大的行
            if (fabs(a[i][c]) > fabs(a[t][c]))
                t = i;
        
        if (fabs(a[t][c]) < eps) continue;
        
        for (int i = c; i <= n; i ++ ) swap(a[t][i], a[r][i]);  // 将绝对值最大的行换到最顶端
        for (int i = n; i >= c; i -- ) a[r][i] /= a[r][c];  // 将当前行的首位变成1
        for (int i = r + 1; i < n; i ++ )  // 用当前行将下面所有的列消成0
            if (fabs(a[i][c]) > eps)
                for (int j = n; j >= c; j -- )
                    a[i][j] -= a[r][j] * a[i][c];
        
        r ++ ;
    }
    
    if (r < n){
        for (int i = r; i < n; i ++ )
            if (fabs(a[i][n]) > eps)
                return 2; // 无解
        return 1; // 有无穷多组解
    }
    
    for (int i = n - 1; i >= 0; i -- )
        for (int j = i + 1; j < n; j ++ )
            a[i][n] -= a[i][j] * a[j][n];
    
    return 0; // 有唯一解
}

```

**01高斯消元**

```c++
int gauss(){  // 高斯消元，答案存于a[i][n]中，0 <= i < n
    int c, r;
    for (c = 0, r = 0; c < n; c ++ ){
        int t = r;
        for (int i = r; i < n; i ++ )  // 找绝对值最大的行
            if (a[i][c])
                t = i;

        if (!a[t][c]) continue;
        swap(a[t], a[r]);
        
        for (int i = r + 1; i < n; i ++ )  // 用当前行将下面所有的列消成0
            if (a[i][c])
                for (int j = n; j >= c; j -- )
                    a[i][j] ^= a[r][j];

        r ++ ;
    }

    if (r < n){
        for (int i = r; i < n; i ++ )
            if (a[i][n])
                return 2; // 无解
        return 1; // 有无穷多组解
    }
    for (int i = n - 1; i >= 0; i -- )
        for (int j = i + 1; j < n; j ++ )
            a[i][n] ^= a[i][j] & a[j][n];

    return 0; // 有唯一解
}
```



#### 		拉格朗日插值法

解决多项式求值问题: 我们有, $n+1$ 个坐标不同的点可以确定唯一的最高为$n$ 次的多项式,

时间复杂度: $O(n^2)$ 

假设多项式为$f(x)$ , 第$i$ 个点的坐标为$(x_{i}, y_{i})$ , 我们需要找到该多项式在$k$ 点的取值, 那么我们根据拉格朗日插值法, 就有

<img src="C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20220714202701315.png" alt="image-20220714202701315" style="zoom: 67%;" />

```C++
int qpow(int a, int k){  // 求a^k mod p
    int res = 1;
    while (k){
        if (k & 1) res = (LL)res * a % mod;
        a = (LL)a * a % mod;
        k >>= 1;
    }
    return res;
}

int lagrange(int n, int *x, int *y, int xi){ //求xi点处的f(xi)值
    int ans = 0;
    for(int i = 1; i <= n; i ++ ){
         int s1 = 1, s2 = 1;
         for(int j = 1; j <= n; j ++ ){
             if(i != j){
                 s1 = 1ll * s1 * (xi - x[j]) % mod;
                 s2 = 1ll * s2 * (x[i] - x[j]) % mod;
             }
         }
         ans = (1ll * ans + 1ll * y[i] * s1 % mod * qpow(s2, mod - 2) % mod) % mod;
    }
    return (ans + mod) % mod;
}
```



**拓展**

1. 在$x$取值连续的时候

这个时候我们可以将算法优化到$O(n)$ 复杂度, 直接将$x_{i}$ 换成$i$ 

我们可以得到
$$
f(k) = \sum_{i=0}^{n}y_{i}\Pi_{i!=j}\frac{k-j}{i-j}
$$
现在考虑如何计算后面的连乘部分, 先考虑分子, 我们维护出分子的前缀积$pre_{i}$和后缀积$suf_{i}$, 则

 

<img src="C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20220714214502061.png" alt="image-20220714214502061" style="zoom: 67%;" />

然后再考虑分母, 可以发现, 分母其实就是阶乘的形式, 对于这个我们直接用$fac[i]$ 表示$i!$ 

那么式子就可以变成

<img src="C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20220714214656346.png" alt="image-20220714214656346" style="zoom: 67%;" />

```C++
int pre[N], suf[N], inv[N];
int lagrange(int n, int *x, int *y, int xi) {//求xi点处的f(xi)值, n表示从0到n的n + 1个点确定的多项式
    if(xi < 0) return 0;
    if(xi <= n) return y[xi];
    int ans = 0;
    pre[0] = (xi - x[0]) % mod, suf[n + 1] = 1;
    for(int i = 1; i <= n; i ++ ) pre[i] = 1ll * pre[i - 1] * (xi - x[i] + mod) % mod;
    for(int i = n; i >= 0; i -- ) suf[i] = 1ll * suf[i + 1] * (xi - x[i]) % mod;

    inv[0] = inv[1] = 1;
    for(int i = 2; i <= n; i ++ ) inv[i] = -1ll * mod / i * inv[mod % i] % mod;
    for(int i = 2; i <= n; i ++ ) inv[i] = 1ll * inv[i] * inv[i - 1] % mod;

    for(int i = 0; i <= n; i ++ ) {
        ans = (ans + (1ll * y[i] * (i == 0 ? 1 : pre[i - 1])) % mod * suf[i + 1] % mod
                     * inv[i] % mod * (((n - i) & 1) ? -1 : 1) * inv[n - i] % mod + mod) % mod;
    }
    return (ans + mod) % mod;
}
```

#### 		矩阵快速幂

```C++
struct Matrix {
    int m, n;
    vector<vector<int>> data;

    Matrix() : m(0), n(0) { data.resize(0); }

    Matrix(int row, int col, int num = 0) : m(row), n(col) {
        data.resize(m, vector<int>(n, num));
    }

    Matrix(int n) : m(n), n(n) {
        data.resize(n, vector<int>(n, 0));
        for (int i = 0; i < n; i++)
            data[i][i] = 1;
    }

    Matrix operator+(const Matrix &b) const {
        if (m != b.m || n != b.n) return Matrix();
        Matrix ans(m, n);
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                ans.data[i][j] = (data[i][j] + b.data[i][j]) % mod;
        return ans;
    }

    Matrix operator-(const Matrix &b) const {
        if (m != b.m || n != b.n) return Matrix();
        Matrix ans(m, n);
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                ans.data[i][j] = (data[i][j] - b.data[i][j]) % mod;
        return ans;
    }

    Matrix operator*(const int &b) const {
        Matrix ans(m, n);
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                ans.data[i][j] = (data[i][j] * b) % mod;
        return ans;
    }

    Matrix operator*(const Matrix &b) const {
        if (n != b.m) return Matrix();
        Matrix ans(m, b.n);
        for (int i = 0; i < m; i++)
            for (int j = 0; j < b.n; j++)
                for (int k = 0; k < n; k++)
                    ans.data[i][j] = (ans.data[i][j] + data[i][k] * b.data[k][j] % mod) % mod;
        return ans;
    }

    Matrix operator%(const int &b) const {
        Matrix ans(m, n);
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                ans.data[i][j] = data[i][j] % b;
        return ans;
    }

    vector<int> &operator[](int rank) { return data[rank]; }

    int sum() const {
        int ans = 0;
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                ans = (ans + data[i][j]) % mod;
        return ans;
    }

    int mul() const {
        int ans = 1;
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                ans = ans * data[i][j] % mod;
        return ans;
    }


    Matrix inv() {
        Matrix res(n), a = *this;
        for (int i = 0, r; i < n; ++i) {
            r = i;
            for (int j = i + 1; j < n; ++j)
                if (a.data[j][i] > a.data[r][i]) r = j;
            swap(a.data[i], a.data[r]), swap(res.data[i], res.data[r]);
            if (!a.data[i][i]) return res.data[0][0] = -1, res;//矩阵的秩<n,无逆元

            long long invaii = qmi(a.data[i][i], mod - 2, mod);
            for (int k = 0; k < n; ++k) a.data[i][k] = a.data[i][k] * invaii % mod;
            for (int k = 0; k < n; ++k) res.data[i][k] = res.data[i][k] * invaii % mod;
            for (int j = 0; j < n; ++j)
                if (j != i) {
                    long long tmp = a.data[j][i];
                    for (int k = i; k < n; ++k)
                        a.data[j][k] = (a.data[j][k] - tmp * a.data[i][k] % mod + mod) % mod;
                    for (int k = 0; k < n; ++k)
                        res.data[j][k] = (res.data[j][k] - tmp * res.data[i][k] % mod + mod) % mod;
                }
        }
        return res;
    }

    friend Matrix qmi(Matrix a, int b) {
        if (a.m != a.n) return Matrix();
        Matrix ans(a.n);
        for (; b; a = a * a, b >>= 1)
            if (b & 1)
                ans = ans * a;
        return ans;
    }

    friend istream &operator>>(istream &in, Matrix &x) {
        in >> x.m >> x.n;
        x.data.resize(x.m, vector<int>(x.n));
        for (int i = 0; i < x.m; i++)
            for (int j = 0; j < x.n; j++)
                in >> x.data[i][j];
        return in;
    }

    friend ostream &operator<<(ostream &out, const Matrix &x) {
        for (int i = 0; i < x.m; i++)
            for (int j = 0; j < x.n; j++)
                out << x.data[i][j] << (j == x.n - 1 ? '\n' : ' ');
        return out;
    }

    bool operator==(const Matrix &b) {
        if (m != b.m || n != b.n) return false;
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                if (data[i][j] != b.data[i][j]) return false;
        return true;
    }

    bool operator!=(const Matrix &b) {
        if (m != b.m || n != b.n) return true;
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                if (data[i][j] != b.data[i][j]) return true;
        return false;
    }
};
```

#### 		线性基

应用: 在一些数中选择某些数, 使得异或和最大

板子:

1. 不显式地构造出集合$A$, 支持动态插入

```c++
struct LinearBasis{
    long long a[MAXL + 1];

    LinearBasis(){
        std::fill(a, a + MAXL + 1, 0);
    }

    LinearBasis(long long *x, int n){
        build(x, n);
    }

    void insert(long long t){
        for (int j = MAXL; j >= 0; j--){
            if (!t) return;
            if (!(t & (1ll << j))) continue;

            if (a[j]) t ^= a[j];
            else{
                for (int k = 0; k < j; k++) if (t & (1ll << k)) t ^= a[k];
                for (int k = j + 1; k <= MAXL; k++) if (a[k] & (1ll << j)) a[k] ^= t;
                a[j] = t;
                break;
            }
        }
    }

    // 数组 x 表示集合 S，下标范围 [1...n]
    void build(long long *x, int n){
        std::fill(a, a + MAXL + 1, 0);

        for (int i = 1; i <= n; i++){
            insert(x[i]);
        }
    }

    long long queryMax(){
        long long res = 0;
        for (int i = 0; i <= MAXL; i++) res ^= a[i];
        return res;
    }

    void mergeFrom(const LinearBasis &other){
        for (int i = 0; i <= MAXL; i++) insert(other.a[i]);
    }

    static LinearBasis merge(const LinearBasis &a, const LinearBasis &b){
        LinearBasis res = a;
        for (int i = 0; i <= MAXL; i++) res.insert(b.a[i]);
        return res;
    }
};
```

2. 直接给出集合$A$ , 不支持动态插入

```c++
struct LinearBasis{
    std::vector<long long> v;   //v[i]表示最高位为1的位为i的线性基(如果没有则为0)
    const int MAXL = 32;
    int n; // 原有集合 S 的大小

    // 数组 x 表示集合 S，下标范围 [1...n]
    void build(long long *x, int n){
        this->n = n;
        std::vector<long long> a(MAXL + 1);

        for (int i = 1; i <= n; i++){
            long long t = x[i];

            for (int j = MAXL; j >= 0; j--){
                if ((t & (1ll << j)) == 0) continue;

                if (a[j]) t ^= a[j];
                else{
                    for (int k = 0; k < j; k++) if (t & (1ll << k)) t ^= a[k];
                    for (int k = j + 1; k <= MAXL; k++) if (a[k] & (1ll << j)) a[k] ^= t;

                    a[j] = t;
                    break;
                }
            }
        }

        v.clear();
        for (int i = 0; i <= MAXL; i++) v.push_back(a[i]);
    }

    long long queryMax(){       //查询最大值
        long long x = 0;
        for (size_t i = 0; i < v.size(); i++) x ^= v[i];
        return x;
    }

    bool query(LL x){       //查询是否存在某个数可以异或得到
        for(int i = MAXL; i >= 0; i -- ){
            if(x >> i & 1) x ^= v[i];
        }
        return x == 0;
    }

};
```



### 	多项式

#### 		FFT

```c++
#define PI acos(-1.0)
struct Complex {
    double x, y;

    Complex(double xx = 0, double yy = 0) { x = xx, y = yy; }

    Complex friend operator+(Complex a, Complex b) {
        return Complex(a.x + b.x, a.y + b.y);
    }

    Complex friend operator-(Complex a, Complex b) {
        return Complex(a.x - b.x, a.y - b.y);
    }

    Complex friend operator*(Complex a, Complex b) {
        return Complex(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
    }
};
typedef Complex cp;

cp a[N], b[N];
int r[N]; // r[i]表示将i的二进制反转之后的结果

void fft(Complex *p, int op, int n) {
    for (int i = 1; i < n; i++) if (i > r[i]) swap(p[i], p[r[i]]);
    for (int i = 1; i < n; i <<= 1) { //表示操作区间集的每个区间的长度
        cp Wn(cos(PI / i), op * sin(PI / i)); // 单元根
        for (int r = i << 1, j = 0; j < n; j += r) {//表示每个区间集的最右边端位置
            cp w(1, 0); // 幂
            for (int k = 0; k < i; k++, w = w * Wn) { //只遍历左区间，右区间O(1)得到
                cp X = p[j + k], Y = w * p[j + k + i];
                p[j + k] = X + Y;
                p[j + k + i] = X - Y;
            }
        }
    }

    if(op == -1){
        for(int i = 0; i < n; i ++ ) a[i].x /= n;
    }
}

int mult(Complex *a, Complex *b, int m) {
    int n = 0, l = 0;
    for (n = 1; n <= m; n <<= 1) l++; 		//l代表的是多项式长度(n)变为二进制的长度
    for (int i = 0; i < n; i++) r[i] = (r[i >> 1] >> 1) | ((i & 1) << (l - 1));
    fft(a, 1, n); fft(b, 1, n);
    for (int i = 0; i < n; i++) a[i] = a[i] * b[i];
    fft(a, -1, n);
    return l;
}

int main() {
    int n, m;
    cin >> n >> m;
    for (int i = 0; i <= n; i++) cin >> a[i].x;
    for (int i = 0; i <= m; i++) cin >> b[i].x;
    //1 表示从系数变为点值
    //-1 表示从点值变为系数
    mult(a, b, n + m);
    for (int i = 0; i <= 2 * n; i++) cout << (int) (a[i].x + 0.5) << (i == m ? '\n' : ' ');
    return 0;
}

```

#### 		NTT

##### 			1.原根

**阶的定义**

​	当$(a, m) = 1$ 时, 满足同余式$a^n \equiv 1\ (mod\ \ m)$ 的最小的$n$, 我们记$n$ 为$a$ 模 $m$ 的阶, 表示为$\delta_{m}(a)$ 

**阶的性质**

<img src="C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20220903173919811.png" alt="image-20220903173919811" style="zoom:67%;" /><img src="C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20220903173756932.png" alt="image-20220903173756932" style="zoom:67%;" />         



**原根定义** 

​	$(a, m) = 1$ 且$\delta_{m}(a) = \phi(m)$ , 我们就称$a$ 为模$m$ 的原根



![image-20220825125522550](C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20220825125522550.png)



**获得n的所有原根**

获得$n$ 的最小原根$g$ 之后, 我们可以发现, 对于所有的 $x, (x, \phi(m)) = 1$ , $g^x$ 均为原根, 所以模$m$ 的原根有$\phi(\phi(m))$ 个

对于最小的原根$g$ 我们首先有个结论, $g$ 的数量级不超过$m^{0.25}$ , 所以暴力枚举用判定定理判定即可

```c++
LL qmi(LL a, LL k, LL p){
    LL res = 1 % p;
    while (k){
        if (k & 1) res = (LL)res * a % p;
        a = (LL)a * a % p;
        k >>= 1;
    }
    return res;
}

int phi[N];
int primes[N], cnt;
bool st[N], rt[N];	//rt[i] = true代表模i下有原根

void init(int n){
    phi[1] = 1;
    for(int i = 2; i <= n; i ++ ){
        if(!st[i]){
            primes[cnt ++ ] = i;
            phi[i] = i - 1;
        }
        for(int j = 0; primes[j] * i <= n; j ++ ){
            st[primes[j] * i] = 1;
            if(i % primes[j] == 0){
                phi[i * primes[j]] = phi[i] * primes[j];
                break;
            }
            phi[i * primes[j]] = phi[i] * (primes[j] - 1);
        }
    }
    rt[2] = rt[4] = 1;      //根据原根存在定理来筛选
    for(int i = 1; i < cnt; i ++ ){
        for(int j = 1; primes[i] * j <= n; j *= primes[i]) rt[j * primes[i]] = 1;
        for(int j = 2; primes[i] * j <= n; j *= primes[i]) rt[j * primes[i]] = 1;
    }
}

struct Get_Omg {
    LL prime[50000];
    int tot = 0;

    LL gcd(LL a, LL b){
        return b ? gcd(b, a % b) : a;
    }

    void get_prime_factor(LL x) {		//质因数分解
        tot = 0;
        for (LL i = 2; i <= x / i; i ++) {
            if (x % i == 0) {
                prime[++ tot] = i;
                while (x % i == 0) x /= i;
            }
        }
        if (x > 1) prime[++ tot] = x;
    }

    bool check(LL x, LL p, LL ph){				//原根判定定理
        if(qmi(x, ph, p) != 1) return false;
        for(int i = 1; i <= tot; i ++ ){
            if(qmi(x, ph / prime[i], p) == 1) return false;
        }
        return true;
    }

    LL getrt(LL p) {            //获得p的最小原根
        get_prime_factor(phi[p]);
        for (LL a = 1; a < p; a ++) {
            if(check(a, p, phi[p])) return a;
        }
        return 0;
    }

    vector<int> getall(LL p){       //获得p的所有原根
        vector<int> ans;
        if(rt[p]) {
            int g = getrt(p);
            LL k = 1;
            for (int i = 1; i <= phi[p]; i++) {
                k = k * g % p;
                if (gcd(i, phi[p]) == 1) ans.push_back(k);
            }
//        sort(all(ans));       //如果要求有序, 可以排一下序
        }
        return ans;
    }
};
```



##### 			2.指标

![image-20220825131210609](C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20220825131210609.png)

指标满足对数性质, 求的时候可以使用BSGS

例题: 

![image-20220904224843858](C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20220904224843858.png)







##### 			3. NTT模板![image-20220826124036243](C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20220826124036243.png)



![image-20220825140320303](C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20220825140320303.png)

```c++
const int N = 3e6 + 10, mod = 998244353;
int r[N]; // l -> 二进制的位数-1
LL a[N], b[N];

LL qmi(LL a, LL k, LL p) {
    int res = 1 % p;
    while (k) {
        if (k & 1) res = (LL) res * a % p;
        a = (LL) a * a % p;
        k >>= 1;
    }
    return res;
}

void ntt(LL *p, int n, int op) {
    for (int i = 1; i < n; i++) if (i > r[i]) swap(p[i], p[r[i]]);
    for (int i = 1; i < n; i <<= 1) { //表示操作区间集的每个区间的长度
        LL Wn = qmi(3, (mod - 1) / (i << 1), mod); //原根
        if (op == -1) Wn = qmi(Wn, mod - 2, mod); // 求逆元
        for (int r = i << 1, j = 0; j < n; j += r) {//表示每个区间集的最右边端位置
            LL w = 1; //幂
            for (int k = 0; k < i; k++, w = w * Wn % mod) { //只遍历左区间，右区间 O(1)得到
                LL X = p[j + k], Y = w * p[j + k + i] % mod;
                p[j + k] = (X + Y) % mod;
                p[j + k + i] = ((X - Y) % mod + mod) % mod;
            }
        }
    }
    if (op == -1)
        for (int i = 0, Inv = qmi(n, mod - 2, mod); i < n; ++i)
            p[i] = 1LL * p[i] * Inv % mod;
}

void multi(LL *a, LL *b, int m) {
    int n, l = 0;
    for (n = 1; n <= m; n <<= 1) l++;   //l代表长度
    for (int i = 0; i < n; i++)
        r[i] = (r[i >> 1] >> 1) | ((i & 1) << (l - 1));
    ntt(a, n, 1); ntt(b, n, 1);
    //1 表示从系数变为点值
    //-1 表示从点值变为系数
    for (int i = 0; i < n; i++) a[i] = a[i] * b[i] % mod;
    ntt(a, n, -1);
}

int main() {
    int n, m, l;
    cin >> n >> m;
    for (int i = 0; i <= n; i++) cin >> a[i];	//n + 1项
    for (int i = 0; i <= m; i++) cin >> b[i];

    multi(a, b, n + m);
    LL inv = qmi(n, mod - 2, mod);
    for (int i = 0; i <= m; i++) cout << a[i] * inv % mod << " ";
    return 0;
}
```



#### 		NTT大板子

```c++
#include <bits/stdc++.h>
using namespace std;

#define int long long
#define fp(i, a, b) for (int i = (a), i##_ = (b) + 1; i < i##_; ++i)
#define fd(i, a, b) for (int i = (a), i##_ = (b) - 1; i > i##_; --i)
const int N = 1e6 + 5, P = 167772161;
using ll = int64_t;
using Poly = vector<int>;
using MultiPoly = vector<Poly>;

//快读
template <typename T>void read(T& x){
    x = 0;
    int f = 1;
    char ch = getchar();
    while (ch < '0' || ch > '9') { if (ch == '-')f = -1; ch = getchar(); }
    while (ch >= '0' && ch <= '9') { x = x * 10 + ch - '0'; ch = getchar(); }
    x *= f;
}
//二次剩余
/*---------------------------------------------------------------------------*/
class Cipolla {
    int P, I2{};
    using pll = pair<ll, ll>;
#define X first
#define Y second
    ll mul(ll a, ll b) const { return a * b % P; }
    pll mul(pll a, pll b) const { return { (a.X * b.X + I2 * a.Y % P * b.Y) % P, (a.X * b.Y + a.Y * b.X) % P }; }
    template<class T> T POW(T a, int b, T x) { for (; b; b >>= 1, a = mul(a, a)) if (b & 1) x = mul(x, a); return x; }
public:
    Cipolla(int p = 0) : P(p) {}
    pair<int, int> sqrt(int n) {
        int a = rand(), x;
        if (!(n %= P)) return { 0, 0 };
        if (POW(n, (P - 1) >> 1, 1ll) == P - 1) return { -1, -1 };
        while (POW(I2 = ((ll)a * a - n + P) % P, (P - 1) >> 1, 1ll) == 1) a = rand();
        x = (int)POW(pll{ a, 1 }, (P + 1) >> 1, { 1, 0 }).X;
        if (2 * x > P) x = P - x;
        return { x, P - x };
    }
#undef X
#undef Y
};
/*---------------------------------------------------------------------------*/
#define MUL(a, b) (ll(a) * (b) % P)
#define ADD(a, b) (((a) += (b)) >= P ? (a) -= P : 0) // (a += b) %= P
#define SUB(a, b) (((a) -= (b)) < 0 ? (a) += P: 0)  // ((a -= b) += P) %= P

//预处理L以内的逆元(0 ~ L-1)
Poly getInv(int L) { Poly inv(L); inv[1] = 1; fp(i, 2, L - 1) inv[i] = MUL((P - P / i), inv[P % i]); return inv; }
auto inv = getInv(N);

//快速幂
int qpow(ll a, int b = P - 2, ll x = 1) { for (; b; b >>= 1, a = a * a % P) if (b & 1) x = x * a % P; return x; }
/*---------------------------------------------------------------------------*/
namespace NTT {
    const int g = 3;
    Poly Omega(int L) {
        int wn = qpow(g, P / L);
        Poly w(L); w[L >> 1] = 1;
        fp(i, L / 2 + 1, L - 1) w[i] = MUL(w[i - 1], wn);
        fd(i, L / 2 - 1, 1) w[i] = w[i << 1];
        return w;
    }
    auto W = Omega(1 << 23); // 注意这边的size，如果大于3e5，改成23；
    void DIF(int* a, int n) {
        for (int k = n >> 1; k; k >>= 1)
            for (int i = 0, y; i < n; i += k << 1)
                for (int j = 0; j < k; ++j)
                    y = a[i + j + k], a[i + j + k] = MUL(a[i + j] - y + P, W[k + j]), ADD(a[i + j], y);
    }
    void IDIT(int* a, int n) {
        for (int k = 1; k < n; k <<= 1)
            for (int i = 0, x, y; i < n; i += k << 1)
                for (int j = 0; j < k; ++j)
                    x = a[i + j], y = MUL(a[i + j + k], W[k + j]),
                    a[i + j + k] = x - y < 0 ? x - y + P : x - y, ADD(a[i + j], y);
        int Inv = P - (P - 1) / n;
        fp(i, 0, n - 1) a[i] = MUL(a[i], Inv);
        reverse(a + 1, a + n);
    }
}
/*-----------------------------------------------------------*/
namespace FWT {
    void FWTor(Poly& a, bool rev) {
        int n = a.size();
        for (int l = 2, m = 1; l <= n; l <<= 1, m <<= 1)
            for (int j = 0; j < n; j += l) fp(i, 0, m - 1) {
                if (!rev) a[i + j + m] = ADD(a[i + j + m], a[i + j]);
                else a[i + j + m] = SUB(a[i + j + m], a[i + j]);
            }
    }
    void FWTand(Poly& a, bool rev) {
        int n = a.size();
        for (int l = 2, m = 1; l <= n; l <<= 1, m <<= 1)
            for (int j = 0; j < n; j += l) fp(i, 0, m - 1) {
                if (!rev) a[i + j] = ADD(a[i + j], a[i + j + m]);
                else a[i + j] = SUB(a[i + j], a[i + j + m]);
            }
    }
    void FWTxor(Poly& a, bool rev)
    {
        int n = a.size(), inv2 = (P + 1) >> 1;
        for (int l = 2, m = 1; l <= n; l <<= 1, m <<= 1)
            for (int j = 0; j < n; j += l) fp(i, 0, m - 1) {
                int x = a[i + j], y = a[i + j + m];
                if (!rev) a[i + j] = ADD(x, y), a[i + j + m] = SUB(x, y);
                else a[i + j] = MUL(ADD(x, y), inv2), a[i + j + m] = MUL(SUB(x, y), inv2);
            }
    }
}
/*---------------------------------------------------------------------------*/
namespace Polynomial {
    // size确定以及NTT乘法
    int norm(int n) { return 1 << ((int)log2(n - 1) + 1); }
    void norm(Poly& a) { if (!a.empty()) a.resize(norm(a.size()), 0); else a = { 0 }; }
    void DFT(Poly& a) { NTT::DIF(a.data(), a.size()); }
    void IDFT(Poly& a) { NTT::IDIT(a.data(), a.size()); }
    Poly& dot(Poly& a, Poly& b) { fp(i, 0, a.size() - 1) a[i] = MUL(a[i], b[i]); return a; }
 
    // 和整数的乘除运算
    Poly& operator*=(Poly& a, int b) { for (auto& x : a) x = MUL(x, b); return a; }
    Poly operator*(Poly a, int b) { return a *= b; }
    Poly operator*(int a, Poly b) { return b * a; }
    Poly& operator/=(Poly& a, int b) { return a *= qpow(b); }
    Poly operator/(Poly a, int b) { return a /= b; }
 
    // 多项式之间的加减运算
    Poly& operator+=(Poly& a, Poly b) {
        a.resize(max(a.size(), b.size()));
        fp(i, 0, b.size() - 1) ADD(a[i], b[i]);
        return a;
    }
    Poly operator+(Poly a, Poly b) { return a += b; }
    Poly& operator-=(Poly& a, Poly b) {
        a.resize(max(a.size(), b.size()));
        fp(i, 0, b.size() - 1) SUB(a[i], b[i]);
        return a;
    }
    Poly operator-(Poly a, Poly b) { return a -= b; }
 
    // 多项式乘法
    Poly operator*(Poly a, Poly b) {
        int n = a.size() + b.size() - 1, L = norm(n);
        if (a.size() <= 30 || b.size() <= 30) {
            Poly c(n);
            fp(i, 0, a.size() - 1) fp(j, 0, b.size() - 1)
                c[i + j] = (c[i + j] + (ll)a[i] * b[j]) % P;
            return c;
        }
        a.resize(L), b.resize(L);
        DFT(a), DFT(b), dot(a, b), IDFT(a);
        return a.resize(n), a;
    }
 
    // 多项式逆元
    Poly Inv2k(Poly a) { // |a| = 2 ^ k
        int n = a.size(), m = n >> 1;
        if (n == 1) return { qpow(a[0]) };
        Poly b = Inv2k(Poly(a.begin(), a.begin() + m)), c = b;
        b.resize(n), DFT(a), DFT(b), dot(a, b), IDFT(a);
        fp(i, 0, n - 1) a[i] = i < m ? 0 : P - a[i];
        DFT(a), dot(a, b), IDFT(a);
        return move(c.begin(), c.end(), a.begin()), a;
    }
    Poly Inv(Poly a) {
        int n = a.size();
        norm(a), a = Inv2k(a);
        return a.resize(n), a;
    }
 
    // 多项式除法/取模
    Poly operator/(Poly a, Poly b) {
        int k = a.size() - b.size() + 1;
        if (k < 0) return { 0 };
        reverse(a.begin(), a.end());
        reverse(b.begin(), b.end());
        b.resize(k), a = a * Inv(b);
        a.resize(k), reverse(a.begin(), a.end());
        return a;
    }
    pair<Poly, Poly> operator%(Poly a, const Poly& b) {
        Poly c = a / b;
        a -= b * c, a.resize(b.size() - 1);
        return { c, a };
    }
 
    // 多项式求导积分
    Poly deriv(Poly a) {
        fp(i, 1, a.size() - 1) a[i - 1] = MUL(i, a[i]);
        return a.pop_back(), a;
    }
    Poly integ(Poly a) {
        a.push_back(0);
        fd(i, a.size() - 1, 1) a[i] = MUL(inv[i], a[i - 1]);
        return a[0] = 0, a;
    }
 
    // 取ln
    Poly Ln(Poly a) {
        int n = a.size();
        a = deriv(a) * Inv(a);
        return a.resize(n - 1), integ(a);
    }

    // 取exp
    Poly Exp(Poly a) {
        int n = a.size(), k = norm(n);
        Poly b = { 1 }, c, d; a.resize(k);
        for (int L = 2; L <= k; L <<= 1) {
            d = b, b.resize(L), c = Ln(b), c.resize(L);
            fp(i, 0, L - 1) c[i] = a[i] - c[i] + (a[i] < c[i] ? P : 0);
            ADD(c[0], 1), DFT(b), DFT(c), dot(b, c), IDFT(b);
            move(d.begin(), d.end(), b.begin());
        }
        return b.resize(n), b;
    }
 
    // 开根
    Poly Sqrt(Poly a) {
        int n = a.size(), k = norm(n); a.resize(k);
        Poly b = { (new Cipolla(P))->sqrt(a[0]).first, 0 }, c;
        for (int L = 2; L <= k; L <<= 1) {
            b.resize(L), c = Poly(a.begin(), a.begin() + L) * Inv2k(b);
            fp(i, L / 2, L - 1) b[i] = MUL(c[i], (P + 1) / 2);
        }
        return b.resize(n), b;
    }
 
   // 多项式快速幂
    Poly Pow1(Poly& a, int b) { return Exp(Ln(a) * b); } // a[0] = 1, 循环卷积
    Poly Pow2(Poly& a, int b) {
        int n = (a.size() - 1) * b + 1, L = norm(n);
        a.resize(L);
        DFT(a);
        fp(i, 0, L - 1) a[i] = qpow(a[i], b);
        IDFT(a);
        return a;
    }
    Poly Pow(Poly a, int b1, int b2) { // b1 = b % P, b2 = b % phi(P) and b >= n if a[0] > 0
        int n = a.size(), d = 0, k;
        while (d < n && !a[d]) ++d;
        if ((ll)d * b1 >= n) return Poly(n);
        a.erase(a.begin(), a.begin() + d);
        k = qpow(a[0]), norm(a *= k);
        a = Pow1(a, b1) * qpow(k, P - 1 - b2);
        a.resize(n), d *= b1;
        fd(i, n - 1, 0) a[i] = i >= d ? a[i - d] : 0;
        return a;
    }
    
    Poly Sin(Poly &a) {
        int i = qpow(3, (P - 1) / 4);
        Poly x(a * i);
        return (Exp(x) - Exp((P - 1) * x)) * qpow(2 * i % P);
    }

    Poly Cos(Poly &a) {
        int i = qpow(3, (P - 1) / 4);
        Poly x(a * i);
        return (Exp(x) + Exp((P - 1) * x)) * qpow(2);
    }

    Poly ASin(Poly &a) {
        int i = qpow(3, (P - 1) / 4);
        return (P - 1) * i % P * Ln(i * a + Sqrt(Poly{1} - a * a));
    }

    Poly ATan(Poly &a) {
        int i = qpow(3, (P - 1) / 4);
        return i * qpow(2) % P * (Ln(Poly{1} - i * a) - Ln(Poly{1} + i * a));
    }
 
    // Get [x ^ k](f / g)
    int divAt(Poly f, Poly g, ll k) {
        int n = max(f.size(), g.size()), m = norm(n);
        for (; k; k >>= 1) {
            f.resize(m * 2, 0), DFT(f);
            g.resize(m * 2, 0), DFT(g);
            fp(i, 0, 2 * m - 1) f[i] = MUL(f[i], g[i ^ 1]);
            fp(i, 0, m - 1) g[i] = MUL(g[2 * i], g[2 * i + 1]);
            g.resize(m), IDFT(f), IDFT(g);
            for (int i = 0, j = k & 1; i < n; i++, j += 2) f[i] = f[j];
            f.resize(n), g.resize(n);
        }
        return f[0];
    }
 
    // Get a[k] by a[n] = sum c[i] * a[n - i]
    int LinearRecur(Poly a, Poly c, ll k) {
        c[0] = P - 1, a = a * c, a.resize(c.size() - 1);
        return divAt(a, c, k);
    }
 
    //Binary convolution for &^|
    Poly operator|(Poly a, Poly b) {
        int n = std::max(a.size(), b.size()), N = norm(n);
        a.resize(N), FWT::FWTor(a, false);
        b.resize(N), FWT::FWTor(b, false);
        Poly A(N);
        fp(i, 0, N - 1) A[i] = MUL(a[i], b[i]);
        FWT::FWTor(A, true);
        return A;
    }
    Poly operator&(Poly a, Poly b) {
        int n = std::max(a.size(), b.size()), N = norm(n);
        a.resize(N), FWT::FWTand(a, false);
        b.resize(N), FWT::FWTand(b, false);
        Poly A(N);
        fp(i, 0, N - 1) A[i] = MUL(a[i], b[i]);
        FWT::FWTand(A, true);
        return A;
    }
    Poly operator^(Poly a, Poly b) {
        int n = std::max(a.size(), b.size()), N = norm(n);
        a.resize(N), FWT::FWTxor(a, false);
        b.resize(N), FWT::FWTxor(b, false);
        Poly A(N);
        fp(i, 0, N - 1) A[i] = MUL(a[i], b[i]);
        FWT::FWTxor(A, true);
        return A;
    }
}
using namespace Polynomial;
//Multi and Inv for f(x,y)
namespace MultiPolynomial {
    void DFT(MultiPoly& a) {
        int m = a.size(), n = a[0].size();
        for (auto& i : a) Polynomial::DFT(i);
        assert(a[0].size() == n);
        fp(i, 0, n - 1) {
            Poly v(m, 0);
            fp(j, 0, m - 1) v[j] = a[j][i];
            Polynomial::DFT(v);
            fp(j, 0, m - 1) a[j][i] = v[j];
        }
    }
 
    void IDFT(MultiPoly& a) {
        int m = a.size(), n = a[0].size();
        assert(a[0].size() == n);
        for (auto& i : a) Polynomial::IDFT(i);
        fp(i, 0, n - 1) {
            Poly v(m, 0);
            fp(j, 0, m - 1) v[j] = a[j][i];
            Polynomial::IDFT(v);
            fp(j, 0, m - 1) a[j][i] = v[j];
        }
    }
 
    MultiPoly& operator*=(MultiPoly& a, int b) {
        int m = a.size(), n = a[0].size();
        fp(i, 0, m - 1) fp(j, 0, n - 1) a[i][j] = MUL(a[i][j], b);
        return a;
    }
    MultiPoly& operator*=(int b, MultiPoly& a) {
        int m = a.size(), n = a[0].size();
        fp(i, 0, m - 1) fp(j, 0, n - 1) a[i][j] = MUL(a[i][j], b);
        return a;
    }
    MultiPoly operator*(MultiPoly a, int b) { a *= b; return a; }
    MultiPoly operator*(int b, MultiPoly a) { a *= b; return a; }
    MultiPoly& operator/=(MultiPoly& a, int b) {
        int m = a.size(), n = a[0].size(), c = qpow(b);
        fp(i, 0, m - 1) fp(j, 0, n - 1) a[i][j] = MUL(a[i][j], c);
        return a;
    }
    MultiPoly operator/(MultiPoly a, int b) { a /= b; return a; }
 
    MultiPoly operator*(MultiPoly a, MultiPoly b) {
        int m = a.size() + b.size(), n = a[0].size() + b[0].size();
        int M = norm(m), N = norm(n);
        a.resize(M), b.resize(M);
        fp(i, 0, M - 1) a[i].resize(N), b[i].resize(N);
        DFT(a), DFT(b);
        MultiPoly ans(M, Poly(N, 0));
        fp(i, 0, M - 1) fp(j, 0, N - 1) ans[i][j] = MUL(a[i][j], b[i][j]);
        IDFT(ans);
        ans.resize(m);
        fp(i, 0, m - 1) ans[i].resize(n);
        return ans;
    }
 
    MultiPoly Inv2k(MultiPoly a) {
        int n = a[0].size();
        int cur = 1, C = a.size();
        MultiPoly ans(1, Polynomial::Inv(a[0]));
        assert(ans[0].size() == n);
        while (cur < C) {
            MultiPoly a0(4 * cur, Poly(4 * n, 0));
            fp(i, 0, 2 * cur - 1) fp(j, 0, n - 1) a0[i][j] = a[i][j];
            DFT(a0);
            ans.resize(4 * cur);
            for (auto& v : ans) v.resize(4 * n);
            DFT(ans);
            fp(i, 0, 4 * cur - 1) fp(j, 0, 4 * n - 1) ans[i][j] = MUL(ans[i][j], 2 + MUL(P - a0[i][j], ans[i][j]));
            IDFT(ans);
            cur <<= 1;
            ans.resize(cur);
            for (auto& v : ans) v.resize(n);
        }
        return ans;
    }
    MultiPoly Inv(MultiPoly a) {
        int m = a.size(), n = a[0].size();
        int M = norm(m), N = norm(n);
        a.resize(M);
        assert(a.size() == M);
        fp(i, 0, M - 1) a[i].resize(N);
        auto b = Inv2k(a);
        b.resize(m);
        fp(i, 0, m - 1) b[i].resize(n);
        return b;
    }
    MultiPoly operator/(MultiPoly a, MultiPoly b) {
        auto c = Inv(b);
        return a * c;
    }
}
// using namespace MultiPolynomial;
/*---------------------------------------------------------------------------*/
// Poly BerlekampMassey(Poly &a) {
//     Poly c, las, cur;
//     int k, delta, d, w;
//     fp(i, 0, a.size() - 1) {
//         d = P - a[i];
//         fp(j, 0, c.size() - 1) d = (d + (ll) a[i - j - 1] * c[j]) % P;
//         if (!d) continue;
//         if (c.empty()) { k = i, delta = d, c.resize(i + 1); continue; }
//         cur = c, w = POW(delta, P - 2, P - d);
//         if (c.size() < las.size() + i - k) c.resize(las.size() + i - k);
//         SUB(c[i - k - 1], w);
//         fp(j, 0, las.size() - 1) c[i - k + j] = (c[i - k + j] + (ll) w * las[j]) % P;
//         if (cur.size() <= las.size() + i - k) las = cur, k = i, delta = d;
//     }
//     return c;
// }
/*---------------------------------------------------------------------------*/
```



#### 		三模数NTT

```c++
#include <bits/stdc++.h>
using namespace std;
typedef long long LL;

using u64 = unsigned long long;

/*
 * 使用注意: 1. 没有重载输入输出流, 要自己读入然后手动赋值
 *         2. 访问值的时候要使用f[i].val();
 *         3. 有的操作需要将其size变为2的幂次
 *         4. 爆long long的时候使用注释内容
 */ 

namespace P1 {
    typedef int DAT;
    constexpr DAT P = 469762049, 998244353, 1004535809;//三个模数, 复制粘贴即可
    template<typename T>
    constexpr T Pow(T a, DAT b = P - 2) {
        T res = 1;
        for (; b; b /= 2, a *= a)if (b % 2)res *= a;
        return res;
    }//必须是封装后的类型
    struct Z {
        DAT z;
        constexpr DAT norm(DAT x) {
            if (x < 0)x += P;
            if (x >= P)x -= P;
            return x;
        } // -P <= x < 2P
        constexpr Z(DAT x = 0) : z(norm(x)) {}
        constexpr DAT val() const { return z; }
        constexpr Z operator-() const { return Z(P - z); }
        constexpr Z inv() const {
            assert(z != 0);
            return Pow(*this);
        }
        constexpr Z &operator*=(const Z &r) {
            z = (LL) z * r.z % P;
            return *this;
        }// int
        //constexpr Z&operator*=(const Z&r){u64 res=(u64)z*r.z-(u64)((long double)z/P*r.z+0.5L)*P;z=(res<P?res:res+P);return *this;}// long long
        constexpr Z &operator+=(const Z &r) {
            z = norm(z + r.z);
            return *this;
        }
        constexpr Z &operator-=(const Z &r) {
            z = norm(z - r.z);
            return *this;
        }
        constexpr Z &operator/=(const Z &r) { return (*this) *= r.inv(); }
        constexpr friend Z operator*(const Z &l, const Z &r) {
            Z res = l;
            return res *= r;
        }
        constexpr friend Z operator+(const Z &l, const Z &r) {
            Z res = l;
            return res += r;
        }
        constexpr friend Z operator-(const Z &l, const Z &r) {
            Z res = l;
            return res -= r;
        }
        constexpr friend Z operator/(const Z &l, const Z &r) {
            Z res = l;
            return res /= r;
        }
    };

    constexpr Z G(3), INVG = Pow(G);// generator
    typedef vector <Z> Poly;
    int EX2(int n) { return 1 << (32 - __builtin_clz(n - 1)); }
    vector<int> Rev{0};
    Z Invlen(1);
    
    void fft(Poly &a, int type = 1) {// a.size == 2^k   type==1 -> dft   type==-1 -> idft
        int n = a.size();
        if (n != Rev.size()) {
            Rev.resize(n);
            int k = 1 << (__builtin_ctz(n) - 1);
            for (int i = 1; i < n; ++i)Rev[i] = (Rev[i / 2] / 2) | ((i & 1) * k);
            Invlen = Z(n).inv();
        }
        for (int i = 1; i < n; ++i)if (i < Rev[i])swap(a[i], a[Rev[i]]);
        for (int i = 1; i < n; i *= 2) {
            Z wn = Pow((type == 1 ? G : INVG), ((P - 1) / (2 * i))), w(1), t;
            for (int j = 0; j < n; j += 2 * i, w = 1)
                for (int k = j; k < i + j; ++k, w *= wn)
                    t = a[k + i] * w, a[k + i] = a[k] - t, a[k] += t;
        }
        if (type < 0)for (int i = 0; i < n; ++i)a[i] *= Invlen;
    }

    void dot(Poly &a, const Poly &b) {// a.size == b.size
        for (int i = 0; i < int(a.size()); ++i)a[i] *= b[i];
    }

    Poly operator*(Poly a, Poly b) {//for all size
        if (a.size() == 0 || b.size() == 0)return Poly();
        int n = a.size() + b.size() - 1;
        if (min(a.size(), b.size()) <= 8) {
            Poly c(n);
            for (int i = 0; i < a.size(); ++i)
                for (int j = 0; j < b.size(); ++j)
                    c[i + j] += a[i] * b[j];
            return c;
        }
        int m = EX2(n);
        a.resize(m); fft(a);
        b.resize(m); fft(b);
        dot(a, b);
        fft(a, -1); a.resize(n);
        return a;
    }
}


namespace MTT {
    typedef int DAT;
    constexpr DAT P = 19260817;// 注意赋值
    template<typename T>
    constexpr T Pow(T a, DAT b = P - 2) {
        T res = 1;
        for (; b; b /= 2, a *= a)if (b % 2)res *= a;
        return res;
    }//必须是封装后的类型
    struct Z {
        DAT z;
        constexpr DAT norm(DAT x) {
            if (x < 0)x += P;
            if (x >= P)x -= P;
            return x;
        } // -P <= x < 2P
        constexpr Z(DAT x = 0) : z(norm(x)) {}
        constexpr DAT val() const { return z; }
        constexpr Z operator-() const { return Z(P - z); }
        constexpr Z inv() const {/*assert(z!=0);*/return Pow(*this); }
        constexpr Z &operator*=(const Z &r) {
            z = (LL) z * r.z % P;
            return *this;
        }// int
        //Z&operator*=(const Z&r){u64 res=(u64)z*r.z-(u64)((long double)z/P*r.z+0.5L)*P;z=(res<P?res:res+P);return *this;}// long long
        constexpr Z &operator+=(const Z &r) {
            z = norm(z + r.z);
            return *this;
        }
        constexpr Z &operator-=(const Z &r) {
            z = norm(z - r.z);
            return *this;
        }
        constexpr Z &operator/=(const Z &r) { return (*this) *= r.inv(); }
        constexpr friend Z operator*(const Z &l, const Z &r) {
            Z res = l;
            return res *= r;
        }
        constexpr friend Z operator+(const Z &l, const Z &r) {
            Z res = l;
            return res += r;
        }
        constexpr friend Z operator-(const Z &l, const Z &r) {
            Z res = l;
            return res -= r;
        }
        constexpr friend Z operator/(const Z &l, const Z &r) {
            Z res = l;
            return res /= r;
        }
    };

    typedef vector <Z> Poly;

    int EX2(int n) { return 1 << (32 - __builtin_clz(n - 1)); }

    void dot(Poly &a, const Poly &b) {// a.size == b.size
        for (int i = 0; i < int(a.size()); ++i)a[i] *= b[i];
    }

    template<typename T>
    constexpr void exgcd(T a, T b, T &x, T &y) {
        if (b == 0) { x = 1, y = 0; }
        else {
            exgcd(b, a % b, y, x);
            y -= (a / b) * x;
        }
    }

    Poly operator*(Poly a, Poly b) {//for all size
        if (a.size() == 0 || b.size() == 0)return Poly();
        int n = a.size() + b.size() - 1;
        if (min(a.size(), b.size()) <= 12) {
            Poly c(n);
            for (int i = 0; i < a.size(); ++i)
                for (int j = 0; j < b.size(); ++j)
                    c[i + j] += a[i] * b[j];
            return c;
        }

        P1::Poly a1(n);
        for (int i = 0; i < a.size(); ++i)a1[i] = a[i].val();
        P1::Poly b1(n);
        for (int i = 0; i < b.size(); ++i)b1[i] = b[i].val();
        a1 = a1 * b1;

        P2::Poly a2(n);
        for (int i = 0; i < a.size(); ++i)a2[i] = a[i].val();
        P2::Poly b2(n);
        for (int i = 0; i < b.size(); ++i)b2[i] = b[i].val();
        a2 = a2 * b2;

        P3::Poly a3(n);
        for (int i = 0; i < a.size(); ++i)a3[i] = a[i].val();
        P3::Poly b3(n);
        for (int i = 0; i < b.size(); ++i)b3[i] = b[i].val();
        a3 = a3 * b3;

        Poly res(n);

        typedef __int128 LLL;
        static constexpr LLL p = (LLL) P1::P * (LLL) P2::P * (LLL) P3::P;
        static constexpr LLL t1 =
                LLL(P1::Pow(P1::Z(P2::P % P1::P)).val()) * LLL(P1::Pow(P1::Z(P3::P % P1::P)).val()) * P2::P * P3::P % p;
        static constexpr LLL t2 =
                LLL(P2::Pow(P2::Z(P1::P % P2::P)).val()) * LLL(P2::Pow(P2::Z(P3::P % P2::P)).val()) * P1::P * P3::P % p;
        static constexpr LLL t3 =
                LLL(P3::Pow(P3::Z(P1::P % P3::P)).val()) * LLL(P3::Pow(P3::Z(P2::P % P3::P)).val()) * P1::P * P2::P % p;
        for (int i = 0; i < n; ++i)res[i] = (t1 * a1[i].val() + t2 * a2[i].val() + t3 * a3[i].val()) % p % P;
        return res;
    }

    Poly inv2k(const Poly &f) {// f.size == 2^k
        int n = f.size();
        //assert(n==EX2(n));///
        //assert(f[0].val()!=0);///
        Poly g{f[0].inv()};
        g.reserve(n);
        for (int k = 2; k <= n; k *= 2) {
            Poly ff(f.begin(), f.begin() + k);
            ff = (g * g) * ff;
            g.resize(k);
            for (int i = 0; i < k / 2; ++i)g[i] = 2 * g[i] - ff[i];
            for (int i = k / 2; i < k; ++i)g[i] = -ff[i];
        }
        return g;
    }

    Poly inv(Poly f) {// for all size
        int n = f.size(), m = EX2(n);
        f.resize(m);
        f = inv2k(f);
        f.resize(n);
        return f;
    }

    Poly dx(const Poly &f) {// for all size 求导
        if (f.size() == 0)return Poly();
        Poly res(f.size() - 1);
        for (int i = 0; i < res.size(); ++i)res[i] = f[i + 1] * (i + 1);
        return res;
    }

    Poly Inv({Z(1)});

    Poly fdx(const Poly &f) {// for all size    积分
        Poly res(f.size() + 1);
        if (f.size() > Inv.size())
            for (int i = Inv.size() + 1; i <= f.size(); ++i)
                Inv.push_back((P - P / i) * Inv[P % i - 1]);
        for (int i = 0; i < f.size(); ++i)res[i + 1] = f[i] * Inv[i];
        return res;
    }

    Poly ln2k(const Poly &f) {// f.size == 2^k
        //assert(f[0].val()==1);///
        Poly res = dx(f) * inv2k(f);
        res.resize(f.size() - 1);
        return fdx(res);
    }

    Poly exp2k(const Poly &f) {// f.size == 2^k
        int n = f.size();
        Poly g({Z(1)});
        g.reserve(n);
        for (int k = 2; k <= n; k *= 2) {
            //assert(g[0].val()==1);///
            g.resize(k);
            Poly h = ln2k(g);
            h[0] = 1;
            for (int i = 1; i < k; ++i)h[i] = f[i] - h[i];
            g = g * h;
            //assert(g[0].val()==1);///
            g.resize(k);
            //assert(g[0].val()==1);///
        }
        return g;
    }

    Poly sqrt2k(const Poly &f, Z q = Z(1)) {// f.size == 2^k
        int n = f.size();
        Poly g{q};
        for (int k = 2; k <= n; k *= 2) {
            g.resize(k);
            Poly h = Poly(f.begin(), f.begin() + k) * inv2k(g);
            for (int i = 0; i < k; ++i)g[i] = (g[i] + h[i]) * ((P + 1) / 2);
        }
        return g;
    }

    pair <Poly, Poly> mod(const Poly &a, const Poly &b) {// for all size
        const int n = a.size() - 1, m = b.size() - 1, t = n - m + 1;
        assert(t > 0);
        Poly f(t);
        for (int i = 0; i < t; ++i)f[i] = a[n - i];
        Poly g(t);
        for (int i = 0; i < min(t, m + 1); ++i)g[i] = b[m - i];
        g = inv(g) * f;
        Poly q(t);
        for (int i = 0; i < t; ++i)q[i] = g[t - i - 1];
        f = b * q;
        Poly r(m);
        for (int i = 0; i < m; ++i)r[i] = a[i] - f[i];
        return {q, r};
    }
}

//可修改模数MTT
//namespace MTT {
//    typedef int DAT;
//    DAT P;
//    template<typename T>
//     T Pow(T a, DAT b = P - 2) {
//        T res = 1;
//        for (; b; b /= 2, a *= a)if (b % 2)res *= a;
//        return res;
//    }//必须是封装后的类型
//    struct Z {
//        DAT z;
//         DAT norm(DAT x) {
//            if (x < 0)x += P;
//            if (x >= P)x -= P;
//            return x;
//        } // -P <= x < 2P
//         Z(DAT x = 0) : z(norm(x)) {}
//         DAT val() const { return z; }
//         Z operator-() const { return Z(P - z); }
//         Z inv() const {/*assert(z!=0);*/return Pow(*this); }
//         Z &operator*=(const Z &r) {
//            z = (LL) z * r.z % P;
//            return *this;
//        }// int
//        //Z&operator*=(const Z&r){u64 res=(u64)z*r.z-(u64)((long double)z/P*r.z+0.5L)*P;z=(res<P?res:res+P);return *this;}// long long
//         Z &operator+=(const Z &r) {
//            z = norm(z + r.z);
//            return *this;
//        }
//         Z &operator-=(const Z &r) {
//            z = norm(z - r.z);
//            return *this;
//        }
//         Z &operator/=(const Z &r) { return (*this) *= r.inv(); }
//         friend Z operator*(const Z &l, const Z &r) {
//            Z res = l;
//            return res *= r;
//        }
//         friend Z operator+(const Z &l, const Z &r) {
//            Z res = l;
//            return res += r;
//        }
//         friend Z operator-(const Z &l, const Z &r) {
//            Z res = l;
//            return res -= r;
//        }
//         friend Z operator/(const Z &l, const Z &r) {
//            Z res = l;
//            return res /= r;
//        }
//    };
//
//    typedef vector <Z> Poly;
//
//    void read() {
//        cin >> P;
//    }
//    int EX2(int n) { return 1 << (32 - __builtin_clz(n - 1)); }
//
//    void dot(Poly &a, const Poly &b) {// a.size == b.size
//        for (int i = 0; i < int(a.size()); ++i)a[i] *= b[i];
//    }
//
//    template<typename T>
//    constexpr void exgcd(T a, T b, T &x, T &y) {
//        if (b == 0) { x = 1, y = 0; }
//        else {
//            exgcd(b, a % b, y, x);
//            y -= (a / b) * x;
//        }
//    }
//
//    Poly operator*(Poly a, Poly b) {//for all size
//        if (a.size() == 0 || b.size() == 0)return Poly();
//        int n = a.size() + b.size() - 1;
//        if (min(a.size(), b.size()) <= 12) {
//            Poly c(n);
//            for (int i = 0; i < a.size(); ++i)
//                for (int j = 0; j < b.size(); ++j)
//                    c[i + j] += a[i] * b[j];
//            return c;
//        }
//
//        P1::Poly a1(n);
//        for (int i = 0; i < a.size(); ++i)a1[i] = a[i].val();
//        P1::Poly b1(n);
//        for (int i = 0; i < b.size(); ++i)b1[i] = b[i].val();
//        a1 = a1 * b1;
//
//        P2::Poly a2(n);
//        for (int i = 0; i < a.size(); ++i)a2[i] = a[i].val();
//        P2::Poly b2(n);
//        for (int i = 0; i < b.size(); ++i)b2[i] = b[i].val();
//        a2 = a2 * b2;
//
//        P3::Poly a3(n);
//        for (int i = 0; i < a.size(); ++i)a3[i] = a[i].val();
//        P3::Poly b3(n);
//        for (int i = 0; i < b.size(); ++i)b3[i] = b[i].val();
//        a3 = a3 * b3;
//
//        Poly res(n);
//
//        typedef __int128 LLL;
//        static constexpr LLL p = (LLL) P1::P * (LLL) P2::P * (LLL) P3::P;
//        static constexpr LLL t1 =
//                LLL(P1::Pow(P1::Z(P2::P % P1::P)).val()) * LLL(P1::Pow(P1::Z(P3::P % P1::P)).val()) * P2::P * P3::P % p;
//        static constexpr LLL t2 =
//                LLL(P2::Pow(P2::Z(P1::P % P2::P)).val()) * LLL(P2::Pow(P2::Z(P3::P % P2::P)).val()) * P1::P * P3::P % p;
//        static constexpr LLL t3 =
//                LLL(P3::Pow(P3::Z(P1::P % P3::P)).val()) * LLL(P3::Pow(P3::Z(P2::P % P3::P)).val()) * P1::P * P2::P % p;
//        for (int i = 0; i < n; ++i)res[i] = (t1 * a1[i].val() + t2 * a2[i].val() + t3 * a3[i].val()) % p % P;
//        return res;
//    }

using namespace MTT;

```



#### 拆系数FFT（1e9 + 7)

需要`define int long long`

```c++
const int mod = 1e9 + 7, N = 1e6 + 10, M = 1 << 21;

#define poly vector<LL>

namespace conv {
#define ld long double
    const int _ = 1 << 19 | 1; const ld pi = acos(-1);
    struct comp{
        ld x , y; comp(ld _x = 0 , ld _y = 0) : x(_x) , y(_y) {}
        friend comp operator +(comp p , comp q) {
            return comp(p.x + q.x , p.y + q.y);
        }
        friend comp operator -(comp p , comp q) {
            return comp(p.x - q.x , p.y - q.y);
        }
        friend comp operator *(comp p , comp q) {
            return comp(p.x * q.x - p.y * q.y , p.x * q.y + p.y * q.x);
        }
        friend comp operator /(comp p , ld q) {
            return comp(p.x / q , p.y / q);
        }
        friend comp operator ~(comp p) {
            return comp(p.x , -p.y);
        }

    }A[_] , B[_] , C[_] , D[_];
    int N , M , P;

    comp w[_];
    int dir[_] , need;

    void init(int len){
        need = 1;
        while(need < len) need <<= 1;
        for(int i = 1 ; i < need ; ++i)
            dir[i] = (dir[i >> 1] >> 1) | (i & 1 ? need >> 1 : 0);

        for(int i = 1 ; i < need ; i <<= 1){
            w[i] = comp(1 , 0);
            comp wn(cos(pi / i) , sin(pi / i));
            for(int j = 1 ; j < i ; ++j)
                w[i + j] = w[i + j - 1] * wn;
        }
    }

    void DFT(comp *A , int t){
        for(int i = 1 ; i < need ; ++i)
            if(i < dir[i])
                swap(A[i] , A[dir[i]]);

        for(int i = 1 ; i < need ; i <<= 1)
            for(int j = 0 ; j < need ; j += i << 1)
                for(int k = 0 ; k < i ; ++k){
                    comp x = A[j + k];
                    comp y = A[i + j + k] * w[i + k];
                    A[j + k] = x + y;
                    A[i + j + k] = x - y;
                }
        if(t == -1) {
            reverse(A + 1 , A + need);
            for(int i = 0 ; i < need ; ++i)
                A[i] = A[i] / need;
        }
    }

    void DFT(comp *A , comp *B){
        static comp P[_] , Q[_];
        for(int i = 0 ; i < need ; ++i)
            P[i] = A[i] + B[i] * comp(0 , 1);

        DFT(P , 1);
        Q[0] = ~P[0];
        for(int i = 1 ; i < need ; ++i)
            Q[i] = ~P[need - i];
        for(int i = 0 ; i < need ; ++i) {
            A[i] = (P[i] + Q[i]) / 2;
            B[i] = (P[i] - Q[i]) * comp(0 , -0.5);
        }
    }

    void IDFT(comp *A , comp *B){
        static comp P[_];
        for(int i = 0 ; i < need ; ++i)
            P[i] = A[i] + B[i] * comp(0 , 1);

        DFT(P , -1);
        for(int i = 0 ; i < need ; ++i) {
            A[i] = comp(P[i].x , 0);
            B[i] = comp(P[i].y , 0);
        }
    }

    void main(int n, int m, int p, int *a, int *b, int *ans){
        N = n; M = m; P = p;
        for(int i = 0 ; i < need ; ++i)
            A[i] = B[i] = C[i] = D[i] = comp();

        for(int i = 0 ; i <= N ; ++i) {
            int x = a[i];
            A[i].x = x & 32767;
            B[i].x = x >> 15;
        }

        for(int i = 0 ; i <= M ; ++i) {
            int x = b[i];
            C[i].x = x & 32767;
            D[i].x = x >> 15;
        }

        init(N + M + 1);
        DFT(A , B); DFT(C , D);

        static comp A1[_] , B1[_] , C1[_] , D1[_];

        for(int i = 0 ; i < need ; ++i) {
            A1[i] = A[i] * C[i];
            B1[i] = A[i] * D[i];
            C1[i] = B[i] * C[i];
            D1[i] = B[i] * D[i];
        }

        IDFT(A1 , B1); IDFT(C1 , D1);

        for(int i = 0 ; i <= N + M ; ++i) {
            ans[i] = (long long)round(A1[i].x);
            ans[i] += (((long long)(round(B1[i].x) + round(C1[i].x)) % P) << 15) % P;
            ans[i] %= P;
            ans[i] += (((long long)round(D1[i].x) % P) << 30) % P;
            ans[i] %= P;
        }
    }
}

poly operator * (const poly &x, const poly &y) {
    static int f[M], g[M], h[M];
    memset(g, 0, sizeof(g));
    memset(h, 0, sizeof(h));
    for(int i = 0; i < x.size(); i ++) f[i] = x[i];
    for(int i = 0; i < y.size(); i ++) g[i] = y[i];

    conv::main(x.size() - 1, y.size() - 1, mod, f, g, h);

    poly z; z.resize(x.size()+y.size());
    for(int i = 0; i < x.size() + y.size() - 1; i ++) {
        z[i] = h[i];
    }
    return z;
}
```



#### 	多项式全家桶

##### 		1. 下降幂多项式乘法

下降幂多项式的点值与EGF有关, 我们先写出$F(x)$ 的点值EGF:

<img src="C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20221021204757436.png" alt="image-20221021204757436" style="zoom: 67%;" />

我们知道$x^{\underline{n}}$ 的点值EGF: $\sum_{i = 0}^{\infty}{\frac{i^{\underline{n}}}{i!}x^i}$经过推导可以化为$e^xx^n$ , 那么$F(n) = \sum_{i = 0}^{\infty}{F[i]n^{\underline{i}}}$ , 那么就相当于

<img src="C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20221021205351859.png" alt="image-20221021205351859" style="zoom:67%;" />

<img src="C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20221021205408137.png" alt="image-20221021205408137" style="zoom:67%;" />

<img src="C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20221021205421883.png" alt="image-20221021205421883" style="zoom:67%;" />

<img src="C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20221021205455344.png" alt="image-20221021205455344" style="zoom:67%;" />

结论: 将$F(x)$ $G(x)$ 当作普通多项式乘上$e^x$ 就可以变成点值表示, 然后直接相乘, 最后再乘上$e^{-x}$ , 就是答案

```c++
signed main() {
    ios::sync_with_stdio(0);
    cin.tie(0);

    int n, m; cin >> n >> m;

    fac[0] = infac[0] = 1;
    for(int i = 1; i <= m + n + 1; i ++ ) fac[i] = fac[i - 1] * i % P;
    for(int i = 1; i <= n + m + 1; i ++ ) infac[i] = infac[i - 1] * inv[i] % P;

    Poly f(n + 1), g(m + 1);
    for(int i = 0; i <= n; i ++ ) cin >> f[i];
    for(int i = 0; i <= m; i ++ ) cin >> g[i];
    Poly e(n + m + 1), ie(n + m + 1);
    for(int i = 0; i < n + m + 1; i ++ ) e[i] = infac[i], ie[i] = i & 1 ? P - infac[i] : infac[i];
    f = f * e; g = g * e;
    for(int i = 0; i < n + m + 1; i ++ ) f[i] = f[i] * g[i] % P * fac[i] % P;
    f = f * ie;
    for(int i = 0; i < n + m + 1; i ++ ) cout << f[i] << " ";
    cout << endl;

    return 0;
}
```



##### 		2.多项式多点求值

<img src="C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20221024183937531.png" alt="image-20221024183937531" style="zoom:67%;" />

```c++
Poly pre[N], F;
int a[N], len[N], n, m;
int bflim;
int ans[N];

void init(int l, int r, int u) {		//预处理(x - ai)的乘积
    if (l == r) {
        len[u] = 1;
        pre[u].resize(2);
        pre[u][0] = P - a[l], pre[u][1] = 1;
        return;
    }
    int mid = l + r >> 1;
    init(l, mid, u << 1), init(mid + 1, r, u << 1 | 1);
    len[u] = r - l + 1;
    Poly F = pre[u << 1], G = pre[u << 1 | 1];
    pre[u].resize(len[u] + 1);
    if (r - l > bflim) pre[u] = F * G;
    else {
        for (int i = 0; i <= len[u << 1]; ++ i)
            for (int j = 0; j <= len[u << 1 | 1]; ++ j)
                pre[u][i + j] = (pre[u][i + j] + F[i] * G[j] % P) % P;
    }
}

void solve(Poly F, int l, int r, int u, int n) {
    if (r - l <= bflim) {
        LL pw[17];
        int res, x, s1, s2, s3, s4;
        pw[0] = 1;
        for (int j = l; j <= r; j++) {
            res = F[n], x = a[j];
            int i = 1;
            for (; i <= 16; i++) pw[i] = pw[i - 1] * x % P;
            i = n - 1;
            while (i >= 15) {
                s1 = res * pw[16] + F[i] * pw[15] + F[i - 1] * pw[14] + F[i - 2] * pw[13];
                s2 = F[i - 3] * pw[12] + F[i - 4] * pw[11] + F[i - 5] * pw[10] + F[i - 6] * pw[9];
                s3 = F[i - 7] * pw[8] + F[i - 8] * pw[7] + F[i - 9] * pw[6] + F[i - 10] * pw[5];
                s4 = F[i - 11] * pw[4] + F[i - 12] * pw[3] + F[i - 13] * pw[2] + F[i - 14] * x;
                res = ((F[i - 15] + s1 + s2) % P + s3 + s4) % P;
                i -= 16;
            }
            i = (n & 15) - 1;
            for(; ~i; --i) res = (res * x + F[i]) % P;
            ans[j] = res;
        }
        return ;
    }
    int mid = l + r >> 1;
    auto res = F % pre[u << 1];
    Poly G = res.second;
    solve(G, l, mid, u << 1, len[u << 1] - 1);
    res = F % pre[u << 1 | 1];
    G = res.second;
    solve(G, mid + 1, r, u << 1 | 1, len[u << 1 | 1] - 1);
}

void Multi_calc() {
    bflim = log2(m);
    init(1, m, 1);
    solve(F, 1, m, 1, n);
}

signed main() {
    ios::sync_with_stdio(0);
    cin.tie(0);

    cin >> n >> m;
    F.resize(n + 1);
    for (int i = 0; i <= n; i++) cin >> F[i];	//多项式
    for (int i = 1; i <= m; i++) cin >> a[i];	//要求值的点
    
    Multi_calc();		//点值答案储存在ans[]中, 下标为1 - m
    
    for(int i = 1; i <= m; i ++ ) cout << ans[i] << endl;
    return 0;
}
```



##### 3.一阶线性微分方程

<img src="C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20221024222613596.png" alt="image-20221024222613596" style="zoom:67%;" />

```c++
Poly F;
void solve(Poly A, Poly B, int n) {
    static int stk[N], top = 0;
    while(n) {
        stk[++ top] = n;
        n >>= 1;
    }
    F[0] = 1;
    while(top -- ) {
        n = stk[top + 1];
        F.resize(n + 1);
        Poly DG = A * Exp(F - Poly{1});
        Poly G = DG + B;
        Poly _DG(n + 1);
        for(int i = 0; i <= n; i ++ ) _DG[i] = (DG[i] == 0 ? 0 : P - DG[i]);
        Poly P = Exp(integ(_DG));
        F = Inv(P) * (integ(P * (G - DG * F)) + Poly{1});
    }
}

signed main() {
    ios::sync_with_stdio(0);
    cin.tie(0);

    int n; cin >> n;
    Poly A(n + 1), B(n + 1);
    for(int i = 0; i <= n; i ++ ) cin >> A[i];
    for(int i = 0; i <= n; i ++ ) cin >> B[i];
    F.resize(n + 1);
    solve(A, B, n);
    for(int i = 0; i <= n; i ++ ) cout << F[i] << " ";
    cout << endl;
    return 0;
}
```



##### 		4.多项式快速插值

<img src="C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20221024224234643.png" alt="image-20221024224234643" style="zoom:67%;" />

```c++
Poly pre[N], H[N];
int a[N], len[N], n, m;
int bflim;
int ans[N], x[N], y[N];

void init(int l, int r, int u) {    //对于(x - ai)乘积的预处理
    if (l == r) {
        len[u] = 1;
        pre[u].resize(2);
        pre[u][0] = P - a[l], pre[u][1] = 1;
        return;
    }
    int mid = l + r >> 1;
    init(l, mid, u << 1), init(mid + 1, r, u << 1 | 1);
    len[u] = r - l + 1;
    Poly F = pre[u << 1], G = pre[u << 1 | 1];
    pre[u].resize(len[u] + 1);
    if (r - l > bflim) pre[u] = F * G;
    else {
        for (int i = 0; i <= len[u << 1]; ++ i)
            for (int j = 0; j <= len[u << 1 | 1]; ++ j)
                pre[u][i + j] = (pre[u][i + j] + F[i] * G[j] % P) % P;
    }
}

void solve(Poly F, int l, int r, int u, int n) {		//多点求值, 答案储存在ans中, 下标1 - n
    if (r - l <= bflim) {
        LL pw[17];
        int res, x, s1, s2, s3, s4;
        pw[0] = 1;
        for (int j = l; j <= r; j++) {
            res = F[n], x = a[j];
            int i = 1;
            for (; i <= 16; i++) pw[i] = pw[i - 1] * x % P;
            i = n - 1;
            while (i >= 15) {
                s1 = res * pw[16] + F[i] * pw[15] + F[i - 1] * pw[14] + F[i - 2] * pw[13];
                s2 = F[i - 3] * pw[12] + F[i - 4] * pw[11] + F[i - 5] * pw[10] + F[i - 6] * pw[9];
                s3 = F[i - 7] * pw[8] + F[i - 8] * pw[7] + F[i - 9] * pw[6] + F[i - 10] * pw[5];
                s4 = F[i - 11] * pw[4] + F[i - 12] * pw[3] + F[i - 13] * pw[2] + F[i - 14] * x;
                res = ((F[i - 15] + s1 + s2) % P + s3 + s4) % P;
                i -= 16;
            }
            i = (n & 15) - 1;
            for(; ~i; --i) res = (res * x + F[i]) % P;
            ans[j] = res;
        }
        return ;
    }
    int mid = l + r >> 1;
    auto res = F % pre[u << 1];
    Poly G = res.second;
    solve(G, l, mid, u << 1, len[u << 1] - 1);
    res = F % pre[u << 1 | 1];
    G = res.second;
    solve(G, mid + 1, r, u << 1 | 1, len[u << 1 | 1] - 1);
}

void Multi_calc() {
    bflim = log2(n);		//这里的长度是点的个数, 本题是n
    init(1, n, 1);
    solve(deriv(pre[1]), 1, n, 1, n);
}

void solve2(int l, int r, int u) {
    if (l == r) {
        len[u] = 1;
        H[u].resize(2);
        H[u][0] = y[l] * qpow(ans[l]) % P, H[u][1] = 0;
        return;
    }
    int mid = l + r >> 1;
    solve2(l, mid, u << 1), solve2(mid + 1, r, u << 1 | 1);
    len[u] = r - l + 1;
    Poly f = H[u << 1], g = H[u << 1 | 1];
    H[u].resize(len[u] + 1);
    H[u] = pre[u << 1 | 1] * f + pre[u << 1] * g;
}

signed main() {
    ios::sync_with_stdio(0);
    cin.tie(0);

    cin >> n;
    for(int i = 1; i <= n; i ++ ) {
        cin >> x[i] >> y[i];
        a[i] = x[i];		//ai代表的是要预处理的多项式, 这里就是x, 然后复制到a里
    }
    Multi_calc();			//多项式多点求值
    solve2(1, n, 1);		//快速插值
    for(int i = 0; i < n; i ++ ) cout << H[1][i] << " ";
    return 0;
}
```

##### 		5.普通多项式转下降幂多项式

![image-20221028155440155](C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20221028155440155.png)

```c++
Poly F, G;
Poly A[N];

void calc(int u, int l, int r) {    //处理下降幂, 运用的仍然是经典的分治做法
    if(l == r) {
        A[u].resize(2);
        A[u][0] = P - l, A[u][1] = 1;
        return ;
    }
    int mid = l + r >> 1;
    calc(u << 1, l, mid), calc(u << 1 | 1, mid + 1, r);
    A[u] = A[u << 1] * A[u << 1 | 1];
}

void solve(int u, int l, int r, Poly F) {
    if(l == r) {
        G[l] = F[0];
        return ;
    }
    int mid = l + r >> 1;
    pair<Poly, Poly> res = F % A[u << 1];
    solve(u << 1, l, mid, res.second);
    solve(u << 1 | 1, mid + 1, r, res.first);
}

signed main() {
    ios::sync_with_stdio(0);
    cin.tie(0);

    int n; cin >> n;
    F.resize(n); G.resize(n);
    for(int i = 0; i < n; i ++ ) cin >> F[i];
    calc(1, 0, n - 1);
    solve(1, 0, n - 1, F);
    for(int i = 0; i < n; i ++ ) cout << G[i] << " ";
    cout << endl;
    return 0;
}
```

##### 		6.下降幂多项式转普通多项式

![image-20221028180633605](C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20221028180633605.png)

```c++
Poly F, G;
Poly g[N], f[N];

void init(int u, int l, int r) {
    if(l == r) {
        g[u].resize(2);
        g[u][0] = P - l; g[u][1] = 1;
        return ;
    }

    int mid = l + r >> 1;
    init(u << 1, l, mid), init(u << 1 | 1, mid + 1, r);
    g[u] = g[u << 1] * g[u << 1 | 1];
}

void solve(int u, int l, int r) {
    if(l == r) {
        f[u].resize(1);
        f[u][0] = F[l];
        return ;
    }

    int mid = l + r >> 1;
    solve(u << 1, l, mid); solve(u << 1 | 1, mid + 1, r);
    f[u] = f[u << 1] + f[u << 1 | 1] * g[u << 1];
}

signed main() {
    ios::sync_with_stdio(0);
    cin.tie(0);

    int n; cin >> n;
    F.resize(n);
    for(int i = 0; i < n; i ++ ) cin >> F[i];
    init(1, 0, n - 1);
    solve(1, 0, n - 1);
    for(int i = 0; i < n; i ++ ) cout << f[1][i] << " ";
    cout << endl;
    return 0;
}
```

##### 		7.高阶前缀和

对于一个序列$A$ , 计算其高维前缀和$A^{(k)}$ 的第$i$ 项其实就相当于计算$A$ 中的每一项在第$i$ 项的贡献；我们发现$A^{(k)}_i$  一定是由$A_1...A_i$ 构成的，如果我们确定他们的系数，就可以确定$A^{(k)}_i$ 的值

首先引入dp式子, $dp_{i, k}$ 表示$A^{(k)}_i$ 的系数构成，这样的话可以发现转移式子为$dp_{i,j} = dp_{i-1,j} + dp_{i,j-1}$ 然后等价为$n * n$ 的网格中从$(1,1)$走到$(i,j)$ 的方案数, 这个时候发现可以用组合数来表示, 每一项的系数为$C_{k}^{0}, C_{k+1}^{1}, ....C_{k+i}^{i}$ , 然后$A^{(k)}_i$ 的系数：$A_{m}$ 的系数为$dp_{i-m,j}$,  所以我们对其进行卷积即可得到答案 

<img src="C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20221207211933394.png" alt="image-20221207211933394" style="zoom:50%;" />

```C++
signed main() {
    ios::sync_with_stdio(0);
    cin.tie(0);
	// 求一个长度为n的序列的k阶前缀和
    int n, k; cin >> n >> k;
    Poly a(n);
    for(int i = 0; i < n; i ++ ) cin >> a[i];
    Poly ki(n);
    k = k % P;
    ki[0] = 1;
    for(int i = 1; i < n; i ++ ) ki[i] = ki[i - 1] * inv[i] % P * ((k + i - 1) % P) % P;
    a = a * ki;
    for(int i = 0; i < n; i ++ ) cout << a[i] << " ";
    cout << endl;

    return 0;
}
```

##### 8.多项式平移

给定$f(x) = \sum_{i = 0}^{n}{f_ix^i}$ , 求$f(x + c)$

+ **利用二项式定理推导系数关系**

<img src="C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20230526141511282.png" alt="image-20230526141511282" style="zoom: 67%;" />

+ **分治法直接求解**

  <img src="C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20230526141836989.png" alt="image-20230526141836989" style="zoom:67%;" />

  **性质:**

  Ⅰ. 如果我们想求解$c$, 使得$A(x + c) = B(x)$ , 这个时候使用二项式定理推导出来的系数关系, 可以得到

  $A(x + c) = \sum_{i= 0}^{n}{a^i\sum_{j = 0}^{i}{C_{i}^{j}x^jc^{i - j}}}$

  $=\sum_{j = 0}^{n}{\sum_{i = j}^{n}{a^iC_{i}^{j}c^{i-j}x^j}}$

​		于是可以得到$b_i = \sum_{j = i}^{n}{C_{j}^{i}a^jc^{j - i}}$ 

​		这个时候我们如果令$i = n - 1$ 可以发现一个系数关系: $b_{n - 1}= a_{n - 1} + a_{n} * n * c$	

​		所以我们只要知道两个多项式的俩系数就可以可以知道$c$ 

#### 	多项式牛顿迭代

![image-20220826143020142](C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20220826143020142.png)

##### 		多项式求逆/开方

![image-20220826143148895](C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20220826143148895.png)



##### 		多项式对数/指数

<img src="C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20220903185601441.png" alt="image-20220903185601441" style="zoom:67%;" />

对于$lnf(x)$ , 可以对$f(x)$ 求导再积分获得

<img src="C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20220903185646569.png" alt="image-20220903185646569" style="zoom:67%;" />

而对于$exp\ f(x)$ 可以使用多项式牛顿迭代

##### 		多项式除法/取模

<img src="C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20221023174754479.png" alt="image-20221023174754479" style="zoom:67%;" />\

##### 		多项式快速幂(加强版)

```c++
signed main() {
    int n, k1 = 0, k2 = 0, k3 = 0;
    string K;
    cin >> n >> K;
    Poly f(n);
    for(int i = 0; i < n; i ++ ) cin >> f[i];
    int l = 0; while(f[l] == 0) l ++ ;
    for(int i = 0; i < K.size(); i ++ ) {
        k1 = (k1 * 10 + K[i] - '0') % P;
        k2 = (k2 * 10 + K[i] - '0') % (P - 1);
        if(k3 * 10 + K[i] - '0' <= P) k3 = k3 * 10 + K[i] - '0';
    }
    if(k3 * l >= n) {
        for(int i = 0; i < n; i ++ ) cout << 0 << " ";
        return 0;
    }
    Poly g = Pow(f, k1, k2);
    for(int i = 0; i < n; i ++ ) cout << g[i] << " ";
    return 0;
}
```



#### 	生成函数

斐波那契数列的生成函数为：$F(x) = \sum{f_ix^i} = \frac{x}{1 - x - x^2}$ 

##### 不同根的有理展开定理

<img src="C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20230414170209618.png" alt="image-20230414170209618" style="zoom:67%;" />

##### 						形式幂级数的其他运算

<img src="C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20220905220433940.png" alt="image-20220905220433940" style="zoom: 50%;" />



$\frac{1}{(1−x)^k}=∑_{i = 0}^{k}C_{k+i−1}^{i}x^i$



##### 		整数分拆(分配问题)

<img src="C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20220907180008369.png" alt="image-20220907180008369" style="zoom:67%;" />

**前置知识**

​	**k分拆数:**

**$n$ 个无标号的球分配到$k$ 个无标号的盒子, 且每个盒子非空的方案数**

这个时候不能简单的利用隔板法, 因为盒子无区分, 而且我们不知道具有相同球数量的盒子的数量, 无法除以排列

![image-20220906160952637](C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20220906160952637.png)

**递推关系**: 如果最后一行有一球, 那么前面的方案数为$p(n-1, k - 1)$ ; 如果最后一行有大于一个球, 那么我们可以选择每个盒子都去掉一个球, 就可以转换成$P(n-k,k)$ 

**生成函数**: 我们将每个盒子中的球拿出来排成一行, 每个盒子占一列, 遵循球的数量从上到下依次减少, 如图所示

<img src="C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20220907112406965.png" alt="image-20220907112406965" style="zoom:50%;" />

考虑有$i$ 个球的**列**的数量的生成函数$f_i(x) = 1 + x^i + x^{2i} + ...$  = $\frac{1}{1 - x^i}$ , 我们知道一列最多有$k$ 个球最少有1个球, 所以乘起来为

$\prod_{i=1}^k{\frac{1}{1 - x^i}}$ , 接着注意到我们要保证有一列必须为$k$ (有$k$ 个球) , 所以我们最后答案为: 

$\prod_{i=1}^k{\frac{1}{1 - x^i}} - \prod_{i=1}^{k-1}{\frac{1}{1 - x^i}}$  

= $(\frac{1}{1 - x^k} - 1)\prod_{i=1}^{k-1}{\frac{1}{1 - x^i}}$  

= $x^k\prod_{i=1}^k{\frac{1}{1 - x^i}}$



**分拆数**

通常情况下我们需要求:  $n$ 个无标号的球放入一些无标号的盒子, 盒子非空的方案数, 记为$p(n)$ , 易知这就是整数分拆本题答案

![image-20220907103633588](C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20220907103633588.png)

**生成函数**: 考虑每一**列**的球的数量的生成函数, 我们发现盒子数量是不确定的, 也就是行的数量不确定, 这个时候转化为**列**的意思就是每一列的球的数量不确定, 可以为无穷, 所以长度为$i$的列的数量的生成函数为$f_i(x) = 1 + x^i + x^{2i} + ...$ = $\frac{1}{1 - x^i}$

所以生成函数为$\prod_{i>=1}{\frac{1}{1 - x^i}}$ , 意思就是每一行至少有一个球, 至多有无限个球, 然后我们的第$n$ 项的意义就是有$n$ 个球的方案数, 行的数量我们并不关系

可以$O(nlogn + n\sqrt{n})$ 求解

```c++
signed main(){
//    DEBUG();
    int n, k; read(n);
    init();

    Poly f;
    f.resize(n + 1);
    for(int i = 1; i <= n; i ++ )
        for(int j = 1; j * i <= n; j ++ ) f[j * i] = f[j * i] + inv[j];

    f = f.exp(n + 1);
    for(int i = 1; i <= n; i ++ ) cout << f[i] << endl;

    return 0;
}
```



**递推关系** :由生成函数推出来, 证明略, 可以$O(n\sqrt{n})$ 求解

```c++
LL f[N];
vector<PII> V;
signed main() {
//    DEBUG();

    int n, k; read(n);

    for(int i = 1; i < N; i ++ ) {
        if(i * (3 * i - 1) / 2 > N) break;
        V.push_back({i * (3 * i - 1) / 2, 1});
    }

    for(int i = 1; i < N; i ++ ){
        if(i * (3 * i + 1) / 2 > N) break;
        V.push_back({i * (3 * i + 1) / 2, 1});
    }
    sort(all(V));

    for(int i = 0; i < V.size(); i ++ ){
        if(i % 4 <= 1) V[i].second = 1;
        else V[i].second = - 1;
    }
    f[0] = f[1] = 1;

    for(int i = 2; i <= n; i ++ ) {
        for(int j = 0; j < V.size(); j ++ ){
            if(V[j].first > i) break;
            f[i] = (f[i] + (f[i - V[j].first] * V[j].second % mod + mod) % mod) % mod;
        }
    }
    for(int i = 1; i <= n; i ++ ) cout << f[i] << endl;

    return 0;
}
```





**总结**:

把$n$ 个球放入$k$ 个盒子的方案:

![image-20220907110150702](C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20220907110150702.png)

第一行: 盒子不为空的话,第二类斯特林数就要求元素有标号, 盒子无标号且非空, 所以乘上盒子标号即可

第二行: 可以为空的话, 可以枚举第二类斯特林数, $i$ 代表有$i$ 个盒子非空

第三行: 转化为不定方程问题,

第四行: 分拆数, 如果盒子可以为空, 那么就在每个盒子都补上一个, 求分拆数即可







##### 		分配问题拓展

![image-20220907203224250](C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20220907203224250.png)

![image-20220907203527896](C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20220907203527896.png)











#### 	Pólya定理, Burnside引理

##### 		群对集合的作用

![image-20220908135143555](C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20220908135143555.png)



##### 		Burnside引理

<img src="C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20220908135046977.png" alt="image-20220908135046977" style="zoom:50%;" />

![image-20220908170429127](C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20220908170429127.png)



对于最后一步, 我们考虑有多少个$i\in [0, n)$,  满足$(n, i) = \frac{n}{d}$ , 两边同时除以$\frac{n}{d}$ 原问题可以转化为有多少个$j$ , 其中$j \in [0, d - 1)$ 满足$(d, j) = 1$ , 这和$\phi(d)$ 的定义一致

这和莫比乌斯反演得到的结果一致. 



接下来我们为了简化过程, 引入置换群的轮换指标这个概念

<img src="C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20220908160715726.png" alt="image-20220908160715726" style="zoom: 50%;" />

$x_i$ 代表$G$ 中的元素, 且长度为$i$ , 通俗一点就是长度为$i$ 的环有多少种情况;  一般要自己算, 如果放在上面问题上那么$x_i = m$ ,   

 $b_i$ 代表有$b_i$ 个长度为$i$ 的环

这个时候我们有常见的轮换指标

![image-20220908191158869](C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20220908191158869.png)



这里注意正$n$ 边形的二面体群的计算, 我们选择先给大括号乘$2n$ 最后再总体除去 

二面体轮换群 $<=>$ 旋转 + 翻转之后方案不变



**正方体的置换群**

<img src="C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20220908195641393.png" alt="image-20220908195641393" style="zoom:33%;" />

##### 		Pólya定理 

<img src="C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20220908201609308.png" alt="image-20220908201609308" style="zoom:50%;" />



<img src="C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20220908202244026.png" alt="image-20220908202244026" style="zoom: 50%;" />

<img src="C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20220908203438101.png" alt="image-20220908203438101" style="zoom:50%;" />



### 	数学结论

1. 从网格沿直线从$(1, 1)$ 到 $(n, m)$ 经过的格子数量为$n + m - gcd(n, m)$ 

   

2. **小数循环节定理:** 

   Ⅰ:  $1 \leq b < a$ , 且$a$ 没有2 或者 5 的质因子, 并且$a$ 与$b$ 互质, 那么$\frac{b}{a}$ 的循环节位数$d$ 恰好满足$10^d ≡ 1\ (mod\ a)$ 

   Ⅱ: $1 \leq b < a$ , 且$a$ 没有2 或者 5 的质因子, 并且$a$ 与$b$ 互质, 那么$\frac{b}{a}$ 的循环节位数$d$ |$\phi(a)$ 

   Ⅲ: $n, m > 2$ , 2, 5都不整除$nm$ , 并且$n$ 与$m$ 是互质的正整数, 则$\frac{1}{nm}$ 的循环节位数为$\frac{1}{n}$ 与$\frac{1}{m}$ 的最小公倍数

   Ⅳ: 若$p ≥ 7$ , 且$p$ 为质数, $n, m$ 是任意正整数并且$p$ 不整除$m$ , 那么$\frac{m}{p^n}$ 的循环节有偶数位, 将此循环节分为前后两段, 则两段对应位置的和均为9

   

3. **反素数:**

​		Ⅰ: 对于一个正整数$x$ , 其约数的个数记为$g(x)$ , 例如$g(1) = 1, g(6) = 4$ , 如果某个正整数$x$ 满足对于任意$i \in (0, x)$ , 均有$g(i) < g(x)$ , 则称$x$ 为反素数, 其实也就是约数个数最多的数

​		性质:

​			Ⅰ: 一个反素数的因子必然是从2开始的连续的质数

​			Ⅱ: $p = 2^{t_1}3^{t_2}5^{t_3}...$ 序列$\{t\}$ 一定单调递减, 即$t_1 >= t_2>= t_3..$ 



​	**4.四柱汉诺塔:** 

![image-20221106202349009](C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20221106202349009.png)

其中$f(n)$ 为四柱情况下移动$n$ 层的次数, $g(n)$为三柱下移动$n$ 层的次数

求解方法: 1. 决策单调性. 设最优转移点为$k$ , 那么对于$i \in [1, k)$ , 有$2f(k) + g(n-k) < 2f(i) + g(n-i)$ , 当$n$ 增大时, $f$ 不变

$g(n-i) > g(n-k)$ , 所以前面$k - 1$ 个点不会成为最优转移点, 可以发现最优转移点具有决策单调性:

​				2.结论: 最优转移点为$n - \lfloor{\sqrt{2n + 1}}\rfloor + 1$ 



​	**5.裴蜀定理**

若$a, b$ 为整数, 且$(a, b) = d$ , 那么对于任意整数$x, y$ , $ax + by$ 都一定是$d$ 的倍数, 特别的, 一定存在$x, y$ 使得$ax + by = d$ 成立

推论: 

​		Ⅰ: $ax ≡ b\ mod(m)$ , $x$ 有解的充要条件为$(a, m) | b$



​	**6.**

设素数$P_n$ , 不超过$\sqrt{P_n}$ 的最大素数为$P_m$ , 则$P_{n + 1} - P_n \leq 2P_m$ 

​	

​    **7** .

给定数字$d$ , 在范围$[0, n - 1]$ 内所有$d$ 的倍数mod $n$ 的不同数字个数为$n / gcd(n, d)$



​    **8.**

**经典问题: 对于一个数组$a$, 我们每一次可以选择$k$ 个互不相同的数减去, 问最大操作次数:** 

方法: 反向思考, 如果假设其函数关系为$g(x):=$ 选择$x$ 个互不相同的数的最大操作次数

我们设$f(x):=$ 操作$x$ 次, 最大可以选择的区间长度;  可以发现$f(x) = \frac{\sum_{i = 1}^{n}{min(C_i, x)}}{x}$ , 这个时候反向赋值$g(x)$ 即可

 



# **图论**

## **最短路**

#### **SPFA算法判断负环**

```c++
bool spfa(){  // 如果存在负环，则返回true，否则返回false。
    // 不需要初始化dist数组
    // 原理：如果某条最短路径上有n个点（除了自己），那么加上自己之后一共有n+1个点，
    // 由抽屉原理一定有两个点相同，所以存在环。
    //count为一个经验值，如果循环了n次，可能就会有一个点多次入队，那么就有可能存在负环
    //这个时候我们从父节点遍历一遍看看有无负环即可
    memset(st, 0, sizeof st);
    memset(pre, -1, sizeof pre);

    int hh = 0, tt = 0;

    for (int i = 1; i <= n; i++) q[tt++] = i, st[i] = true;

    auto detectCycle = [&]() {
        vector<int> vec;
        vector<bool> inStack(676, false);
        vector<bool> vis(676, false);
        for (int i = 1; i <= n; i++)
            if (!vis[i]) {
                for (int j = i; j != -1; j = pre[j]) {
                    if (!vis[j]) {
                        vis[j] = true;
                        vec.push_back(j);
                        inStack[j] = true;
                    } else {
                        if (inStack[j]) return true;
                        break;
                    }
                }
                for (int j: vec) inStack[j] = false;
                vec.clear();
            }
        return false;
    };
    int count = 0;
    while (hh != tt) {
        int t = q[hh++];
        if (hh == N) hh = 0;
        st[t] = false;

        for (int i = h[t]; ~i; i = ne[i]) {
            int j = e[i];
            if (dist[j] > dist[t] + w[i]) {
                dist[j] = dist[t] + w[i];
                pre[j] = t;
                if (count >= n) {
                    count = 0;
                    if (detectCycle()) return true;
                }
                if (!st[j]) {
                    st[j] = true;
                    q[tt++] = j;
                    if (tt == N) tt = 0;
                }
            }
        }
    }

    return false;
}
```

```c++
bool spfa(){  // 如果存在负环，则返回true，否则返回false。
    // 不需要初始化dist数组
    // 原理：如果某条最短路径上有n个点（除了自己），那么加上自己之后一共有n+1个点，
    // 由抽屉原理一定有两个点相同，所以存在环。
    int hh = 0, tt = 0;
    for (int i = 1; i <= n; i++) q[tt++] = i, st[i] = true;

    while (hh != tt) {
        int t = q[--tt];
        st[t] = false;

        for (int i = h[t]; ~i; i = ne[i]) {
            int j = e[i];
            if (dist[j] > dist[t] + w[i]) {
                dist[j] = dist[t] + w[i];
                cnt[j] = cnt[t] + 1;
                if (cnt[j] >= n) return true;
                if (!st[j]) {
                    st[j] = true;
                    q[tt++] = j;
                }
            }
        }
    }
    return false;
}

```

#### **Bellman-Ford算法  **

**时间复杂度$O(nm)$**    **注：可求最多经过k条边的最短距离**

```c++
int n, m, k;            //n个点，m条边，最多经过k条边
int dist[N];
int last[N];            //拷贝数组

void bellman_ford() {
    memset(dist, 0x3f, sizeof dist);

    dist[1] = 0;
    for (int i = 0; i < k; i++) {
        memcpy(last, dist, sizeof dist);
        for (int j = 0; j < m; j++) {
            auto e = edges[j];
            dist[e.b] = min(dist[e.b], last[e.a] + e.c);
        }
    }
}

```



### **最小生成树**

**注：n为点数， m为边数**

#### **Prim算法**

**朴素版 $O(n^2)$    邻接矩阵（稠密图）**

```c++
int prim(){      //res返回最小生成树的边的总权重, dist[i]表示从i到目前最小生成树的最短距离
    int res = 0;
    memset(dist, 0x3f, sizeof dist);
    dist[1] = 0;
    for (int i = 0; i < n; i++){           //迭代n次，每次加入一个点到最小生成树中
        int t = -1;
        for (int j = 1; j <= n; j++)
            if (!st[j] && (t == -1 || dist[t] > dist[j]))
                t = j;

        res += dist[t];                    //选出这个点t，加入最小生成树
        st[t] = 1;

        for (int j = 1; j <= n; j++) dist[j] = min(dist[j], w[t][j]);
    }

    return res;
}
```

#### **Kruskal算法      **

**$O(mlogm)$     邻接表（稀疏图）**

**特点 ：**

1.算法无论进行到什么时候都是正确的，所以可以用来求“最小生成树林”；

2.选出来的最小生成树具有在所有生成树中“最大的边的权值最小”这个性质

**应用 ：** 

1.拓展完全图，可以每一次合并令两个连通块内所有的点都连一条边权为w + 1的点（w为算法现在枚举的这条最小生成树的边），这样拓展出来的完全图边权最小

##### **1.次小生成树**

朴素做法：先求出最小生成树，然后预处理出（dfs）最小生成树的任意两个点之间的最大边权值和严格次大边权值，我们有结论，次小生成树与最小生成树只有一条边不一样，所以我们只需要枚举剩余的非树边，然后如果比最大值大，就更新最大值，要不然就更新严格次大值，每一次取一个min即可

```c++
struct Edge {
    int a, b, w;
    bool f;

    bool operator<(const Edge &t) const {
        return w < t.w;
    }
} edge[M];

int p[N];
int e[N * 2], ne[N * 2], w[N * 2], h[N], idx;
int dist1[N][N], dist2[N][N];

void dfs(int u, int fa, int maxd1, int maxd2, int d1[],
         int d2[]) {                                       //目的是预处理最小生成树上两个点之间的最大值和严格次大值
    d1[u] = maxd1, d2[u] = maxd2;
    for (int i = h[u]; ~i; i = ne[i]) {
        int j = e[i];
        if (j == fa) continue;

        int td1 = maxd1, td2 = maxd2;
        if (w[i] > td1) td2 = td1, td1 = w[i];
        else if (w[i] < td1 && w[i] > td2) td2 = w[i];
        dfs(j, u, td1, td2, d1, d2);
    }
}

int main() {
    cin >> n >> m;
    memset(h, -1, sizeof h);
    for (int i = 0; i < m; i++) {
        int a, b, c;
        cin >> a >> b >> c;
        edge[i] = {a, b, c};
    }
    sort(edge, edge + m);
    for (int i = 1; i <= n; i++) p[i] = i;

    LL sum = 0;
    for (int i = 0; i < m; i++) {     //做最小生成树，并且把最小生成树用邻接表存起来，标记这些边
        int a = edge[i].a, b = edge[i].b, w = edge[i].w;
        int pa = find(a), pb = find(b);
        if (pa != pb) {
            p[pa] = pb;
            sum += w;
            add(a, b, w), add(b, a, w);
            edge[i].f = 1;
        }
    }

    for (int i = 1; i <= n; i++) dfs(i, -1, -1e9, -1e9, dist1[i], dist2[i]);

    LL res = 1e18;
    for (int i = 0; i < m; i++)
        if (!edge[i].f) {
            int a = edge[i].a, b = edge[i].b, w = edge[i].w;
            LL t;
            if (w > dist1[a][b]) t = sum + w - dist1[a][b];
            else if (w > dist2[a][b]) t = sum + w - dist2[a][b];
            res = min(res, t);
        }
    cout << res << endl;
    return 0;
}
```

**LCA优化做法**

```c++
#include <bits/stdc++.h>
using namespace std;
typedef long long LL;

const int N = 1e6 + 10, M = N * 3, INF = 0x3f3f3f3f;
int n, m;

struct Edge {
    int a, b, w;
    bool used;

    bool operator<(const Edge &t) const &{
        return w < t.w;
    }
} edges[M];

int p[N], h[N], e[M], ne[M], w[M], idx;
int depth[N], fa[N][17], d1[N][17], d2[N][17]; //log2(1e5) = 16;

void add(int a, int b, int c){
    e[idx] = b, w[idx] = c, ne[idx] = h[a], h[a] = idx++;
}

int find(int x){      //并查集
    return x == p[x] ? p[x] : p[x] = find(p[x]);
}

LL kruscal(){    //kruscal算法求最小生成树
    for (int i = 1; i <= n; i++) p[i] = i;
    sort(edges, edges + m);
    LL res = 0;

    for (int i = 0; i < m; i++) {
        int a = find(edges[i].a), b = find(edges[i].b), w = edges[i].w;
        if (a != b) {
            p[a] = b;
            res += w;
            edges[i].used = 1;
        }
    }

    return res;
}

void build(){    //建立最小生成树这个图
    memset(h, -1, sizeof h);
    for (int i = 0; i < m; i++) {
        if (edges[i].used) {
            int a = edges[i].a, b = edges[i].b, w = edges[i].w;
            add(a, b, w), add(b, a, w);
        }
    }
}

void bfs(){      //预处理出depth数组和d2,d1数组
    memset(depth, 0x3f, sizeof depth);
    depth[0] = 0, depth[1] = 1;
    queue<int> q;
    q.push(1);

    while (q.size()) {
        int t = q.front();
        q.pop();
        for (int i = h[t]; ~i; i = ne[i]) {
            int j = e[i];
            if (depth[j] > depth[t] + 1) {
                depth[j] = depth[t] + 1;
                q.push(j);
                fa[j][0] = t;
                d1[j][0] = w[i], d2[j][0] = -INF;

                for (int k = 1; k <= 16; k++) {
                    int anc = fa[j][k - 1];
                    fa[j][k] = fa[anc][k - 1];
                    int distance[4] = {d1[j][k - 1], d2[j][k - 1], d1[anc][k - 1], d2[anc][k - 1]};
                    //最大值和次大值肯定是这四个值之中的一个
                    d1[j][k] = d2[j][k] = -INF;
                    for (int u = 0; u < 4; u++){     //然后遍历求最大值和次大值即可
                        int d = distance[u];
                        if (d > d1[j][k]) d2[j][k] = d1[j][k], d1[j][k] = d;
                            //严格次大值
                        else if (d != d1[j][k] && d > d2[j][k]) d2[j][k] = d;
                    }
                }
            }
        }
    }
}

int lca(int a, int b, int w) {
    static int distance[N * 2];
    int cnt = 0;

    if (depth[a] < depth[b]) swap(a, b);
    for (int k = 16; k >= 0; k--) {
        if (depth[fa[a][k]] >= depth[b]) {
            distance[cnt++] = d1[a][k];
            distance[cnt++] = d2[a][k];
            a = fa[a][k];
        }
    }

    if (a != b) {
        for (int k = 16; k >= 0; k--) {
            if (fa[a][k] != fa[b][k]) {
                distance[cnt++] = d1[a][k];
                distance[cnt++] = d2[a][k];
                distance[cnt++] = d1[b][k];
                distance[cnt++] = d2[b][k];

                a = fa[a][k], b = fa[b][k];
            }
        }

        distance[cnt++] = d1[a][0];
        distance[cnt++] = d1[b][0];
    }

    int dist1 = -INF, dist2 = -INF;
    for (int i = 0; i < cnt; i++) {
        int d = distance[i];
        if (d > dist1) dist2 = dist1, dist1 = d;
        else if (d != dist1 && d > dist2) dist2 = d;
    }

    if (w > dist1) return w - dist1;
    if (w > dist2) return w - dist2;

    return INF;
}

int main() {
    cin >> n >> m;
    for (int i = 0; i < m; i++) {
        int a, b, c;
        cin >> a >> b >> c;
        edges[i] = {a, b, c};
    }

    LL sum = kruscal();
    build();

    bfs();

    LL res = 1e18;

    for (int i = 0; i < m; i++)
        if (!edges[i].used) {
            int a = edges[i].a, b = edges[i].b, w = edges[i].w;
            res = min(res, sum + lca(a, b, w));
        }
    printf("%lld\n", res);

    return 0;
}
```



### 拓扑排序

```c++
void topsort() {
    int hh = 0, tt = -1;
    for (int i = 1; i <= n; i ++ )
        if (!d[i]) q[ ++ tt] = i;

    while (hh <= tt) {
        int t = q[hh ++ ];
        for (int i = h[t]; ~i; i = ne[i]) {
            int j = e[i];
            if ( -- d[j] == 0)
                q[ ++ tt] = j;	//q中0 -> n - 1 为拓扑排序结果
        }
    }
}
```

检查是否是DAG：可以直接拓扑排序之后遍历点，检查度数是否有不为0的。

### LCA(最近公共祖先)

#### 1.倍增求LCA

​		时间复杂度$O(nlogn + logm)$

```C++
vector<int> dep(n, 0x3f3f3f3f);
vector<vector<int>> fa(n, vector<int>(21));
auto bfs = [&](int root) {
    queue<int> q;
    q.push(root);
    dep[root] = 1;
    while (q.size()) {
        int u = q.front(); q.pop();
        for (auto v : g[u]) {
            if (dep[v] > dep[u] + 1) {
                dep[v] = dep[u] + 1;
                q.push(v);
                fa[v][0] = u;
                for (int k = 1; k <= 20; k++) {
                    fa[v][k] = fa[fa[v][k - 1]][k - 1];
                }
            }
        }
    }
};

auto lca = [&](int a, int b) {
    if (dep[a] < dep[b]) swap(a, b);
    for (int k = 20; k >= 0; k--) {
        if (dep[fa[a][k]] >= dep[b]) a = fa[a][k];
    }
    if (a == b) return a;
    for (int k = 20; k >= 0; k--) {
        if (fa[a][k] != fa[b][k]) {
            a = fa[a][k];
            b = fa[b][k];
        }
    }
    return fa[a][0];
};

auto dis = [&](int a, int b) {
    return dep[a] + dep[b] - 2 * dep[lca(a, b)];
};
```



#### Tarjan算法离线求lca

时间复杂度$O(n + m)$

```c++
//两点之间的距离，设d[i]表示i到根节点的距离，那么i到j的距离就是d[i] + d[j] - 2 * d[lca(i ,j)]
int n, m;
int e[M], ne[M], w[M], h[N], idx;
int res[M], p[N], dist[N];
int st[N];
vector <PII> query[N];

void add(int a, int b, int c){
    e[idx] = b, w[idx] = c, ne[idx] = h[a], h[a] = idx++;
}

int find(int x){      //并查集
    return x == p[x] ? p[x] : p[x] = find(p[x]);
}

void dfs(int u, int fa) {
    for (int i = h[u]; ~i; i = ne[i]) {
        int j = e[i];
        if (j == fa) continue;
        dist[j] = dist[u] + w[i];
        dfs(j, u);
    }
}

void tarjan(int u) {
    st[u] = 1;          //1表示正在搜索的点
    for (int i = h[u]; ~i; i = ne[i]) {
        int j = e[i];
        if (!st[j]){     //0表示还没有被搜到的点
            tarjan(j);
            p[j] = u;
        }
    }
    for (auto &[k, id]: query[u]) {
        if (st[k] == 2) {
            int anc = find(k);
            res[id] = dist[k] + dist[u] - 2 * dist[anc];
        }
    }
    st[u] = 2;          //2表示已经搜过的点
}

int main() {
    cin >> n >> m;
    memset(h, -1, sizeof h);

    for (int i = 0; i < n - 1; i++) {
        int a, b, c;
        cin >> a >> b >> c;
        add(a, b, c), add(b, a, c);
    }

    for (int i = 0; i < m; i++) {    //记录每一组的询问
        int a, b;
        cin >> a >> b;
        if (a != b) {
            query[a].push_back({b, i});
            query[b].push_back({a, i});
        }
    }

    for (int i = 1; i <= n; i++) p[i] = i;

    dfs(1, -1);                         //计算每个节点到根节点的距离
    tarjan(1);

    for (int i = 0; i < m; i++) cout << res[i] << endl;
    return 0;
}
```

### 树上问题

#### 1.树的性质与“心”

1. **树的直径**

```c++
int res;
int dfs(int u, int father) {//*dfs的返回值是以u为节点的下面所有路径的最大值
    int d1 = 0, d2 = 0;     //维护最大值和次大值
    for (auto [v, w] : g[u]) {
        if(v == father) continue;
        
        int d = dfs(v, u) + w;
        if (d >= d1) d2 = d1, d1 = d;
        else if (d > d2) d2 = d;
    }
    res = max(res, d1 + d2);

    return d1;
}
```

2. **树的中心**

定义: 树的中心距离其他节点的最远距离最近

```c++
int n;
vector<PII> g[N];
int d1[N], d2[N]; //*d1[i]储存的是节点i往下走距离最远的距离是多少  d2[i]储存的就是第二大的路径
int p1[N];        //*p1[i]储存的是节点i往下走走最远的路径是从哪个节点下去的
int up[N];        //*up[i]储存的是节点i往上走最远是多少


int dfs_d(int u, int father) { //*这个dfs求的是从点u往下走的最大距离
    d1[u] = d2[u] = -INF;
    for (auto [v, w] : g[u]) {
        if (v == father) continue;
        
        int d = dfs_d(v, u) + w;
        if (d >= d1[u]) { //*这边是同时存储从点u往下走所能获得的最大距离和第二大的距离
            d2[u] = d1[u];
            d1[u] = d;
            p1[u] = v; //*注意同时存一下u走最大距离往下走的时候经历的第一个点是什么
        } 
        else if (d > d2[u]) d2[u] = d;
    }

    if (d1[u] == -INF) d1[u] = d2[u] = 0; //*如果两个距离都没有被更新，说明这个点是叶子节点，直接置为0
    return d1[u];
}

//这个dfs求的是从点u往上走所能获得的最大距离，这个问题可以转化为father点 
//往下走所能获得的最大距离，然后这个最大距离如果经过j就用次大距离，如果不经过就用最大距离

void dfs_u(int u, int father) {
    for (auto [v, w] : g[u]) {
        if (v == father) continue;

        if (p1[u] == v) up[v] = max(up[u], d2[u]) + w;
        else up[v] = max(up[u], d1[u]) + w;
        
        dfs_u(v, u);
    }
}

int main() {
    cin >> n;
    for (int i = 0; i < n - 1; i++) {
        int a, b, c;
        cin >> a >> b >> c;
        g[a].push_back({b, c});
        g[b].push_back({a, c});
    }
    dfs_d(1, -1);
    dfs_u(1, -1);

    int res = d1[1];
    for (int i = 2; i <= n; i++)
        res = min(res, max(d1[i], up[i])); //*如果一个点
    cout << res << endl;
    return 0;
}
```

4. **树的重心**

重心定义：重心是指树中的一个结点，如果将这个点删除后，剩余各个连通块中点数的最大值最小，那么这个节点被称为树的重心。

```c++
int n, ans = N;
vector<int> g[N];
bool st[N];

int dfs(int u) {             //* u返回的是以u为根的子树中的所有节点
    st[u] = true;
    int size = 0, sum = 0;
    for (auto v : g[u]) {
        if (st[v]) continue;
        int s = dfs(v);
        size = max(size, s);   //* 除去u节点外，各个连通块（子树）中的点数的最大值
        sum += s;             //* 统计以u为根的子树的数量
    }
    size = max(size, n - sum - 1);
    ans = min(ans, size);
    return sum + 1;           //* 返回以u为根的子树的数量

}

int main() {
    cin >> n;
    for (int i = 0; i < n - 1; i++) {
        int a, b;
        cin >> a >> b;
        g[a].push_back(b);
        g[b].push_back(a);
    }
    dfs(1);
    cout << ans << endl;    //删除重心之后各个连通块中点数的最大值

    return 0;
}
```

5.**树的质心**

对于一棵树，节点$u$ 为质心的条件为：对于$u$ 的所有的子节点$v$ ， 所有以$v$ 为根节点的子树的节点个数最多是以$u$ 为根的子树的一半

性质：每一棵树都有1或者2个质心

```c++
//质心分解：将一棵树T重构为R，使得满足以下两个条件：
// 1. 对于R中每一个节点，对于其子树都是质心（有根树）
// 2. 对于R中每一对节点的lca，都在T的x -> y的路径上
signed main() {
    ios::sync_with_stdio(0);
    cin.tie(0);

    int n; cin >> n;
    vector<vector<int>> g(n);
    for(int i = 1; i < n; i ++ ) {
        int a, b; cin >> a >> b;
        a --, b --;
        g[a].push_back(b);
        g[b].push_back(a);
    }
    vector<int> p(n, -1);
    vector<bool> used(n, false);
    vector<int> sz(n, 0);

    function<void(int, int)> dfs_sz = [&](int u, int fa) {
        sz[u] = 1;
        for(auto v : g[u]) {
            if(v == fa || used[v]) continue;
            dfs_sz(v, u);
            sz[u] += sz[v];
        }
    };

    function<int(int, int, int)> find = [&](int u, int fa, int s) {
        for(auto v : g[u]) {
            if(v == fa || used[v] || sz[v] * 2 <= s) continue;
            return find(v, u, s);
        }
        return u;
    };

    function<int(int)> dfs = [&](int u) {
        dfs_sz(u, -1);
        int x = find(u, -1, sz[u]);
        used[x] = true;
        for(auto v : g[x]) {
            if(!used[v]) {
                int y = dfs(v);
                p[y] = x;
            }
        }
        return x;
    };

    dfs(0);
    for(int i = 0; i < n; i ++ ) {
        if(p[i] == -1) cout << -1 << " ";
        else cout << p[i] + 1 << " ";
    }
    cout << endl;

    return 0;
}
```





#### 2. 点分治

时间复杂度$O(nlogn)$ , 每次选择子树的重心作为根节点，重新选择根节点之后一定要重新计算子树大小

```c++
// 点分治板子，q次询问，查询树上是否存在距离为k的两个点

int sz[N], mxp[N], q[N];
vector<pair<int, int>> g[N];
int n, m, root, dist[N];
bool vis[N], res[N];

void calcSz(int u, int fa, int sum) {
    sz[u] = 1;
    mxp[u] = 0;
    for(auto [v, w] : g[u]) {
        if(v == fa || vis[v]) continue;
        calcSz(v, u, sum);
        mxp[u] = max(mxp[u], sz[v]);
        sz[u] += sz[v];
    }
    mxp[u] = max(mxp[u], sum - sz[u]);
    if(!root || mxp[u] < mxp[root]) root = u;
}

int A[N], cnt, B[N];
// 记当前分治的点为root， A数组记录从root能到的点，dist[i]对应的是A[i]到root的距离
// B数组记录A[i]属于root的哪一颗子树（即当B[A[i]] = B[A[j]] 时，说明A[i]，A[j]
// 属于root的同一棵子树
void calcDist(int u, int fa, int dis, int from) {
    A[++ cnt] = u;
    dist[u] = dis;
    B[u] = from;
    for(auto [v, w] : g[u]) {
        if(v == fa || vis[v]) continue;
        calcDist(v, u, dis + w, from);
    }
}

void dfs(int u) {
    cnt = 0; vis[u] = true;
    A[++ cnt] = u;
    dist[u] = 0;
    B[u] = u;
    for(auto [v, w] : g[u]) {
        if(vis[v]) continue;
        dist[v] = w;
        calcDist(v, u, w, v);
    }
    sort(A + 1, A + cnt + 1, [&](int x, int y) {
        return dist[x] < dist[y];
    });
    for(int i = 0; i < m; i ++ ) {
        int l = 1, r = cnt;
        if(res[i]) continue;
        while(l < r) {
            if(dist[A[l]] + dist[A[r]] > q[i]) r --;
            else if(dist[A[l]] + dist[A[r]] < q[i]) l ++ ;
            else if(B[A[l]] == B[A[r]]) {
                if(dist[A[r]] == dist[A[r - 1]]) r --;
                else l ++ ;
            } else {
                res[i] = true;
                break;
            }
        }
    }

    for(auto [v, w] : g[u]) {
        if(vis[v]) continue;
        root = 0;
        calcSz(v, 0, sz[v]);
        dfs(root);
    }
}


signed main() {
    ios::sync_with_stdio(0);
    cin.tie(0);

    cin >> n >> m;
    for(int i = 0; i < n - 1; i ++ ) {
        int u, v, w; cin >> u >> v >> w;
        g[u].push_back({v, w});
        g[v].push_back({u, w});
    }
    for(int i = 0; i < m; i ++ ) {
        cin >> q[i];
        if(!q[i]) res[i] = true;
    }

    root = 0;
    mxp[root] = INF;
    calcSz(1, -1, n);
    dfs(root);

    for(int i = 0; i < m; i ++ ) {
        if(res[i]) cout << "AYE" << endl;
        else cout << "NAY" << endl;
    }

    return 0;
}
```



#### 3.前序 and 中序 确定二叉树

```c++
int n;
int P[N], I[N], J[N];
int L[N], R[N];

void dfs(bool &ok, int l, int r, int p) {
    if(J[P[p]] < l || J[P[p]] >= r) {
        ok = false;
        return ;
    }
    if(J[P[p]] > l) {
        if(p == n - 1) {
            ok = false;
            return ;
        }
        L[P[p]] = P[p + 1];
        dfs(ok, l, J[P[p]], p + 1);
        if(!ok) return ;
    }

    if(J[P[p]] < r - 1) {
        int q = p + (J[P[p]] - l + 1);
        if(q >= n) {
            ok = false;
            return ;
        }
        R[P[p]] = P[q];
        dfs(ok, J[P[p]] + 1, r, q);
        if(!ok) return ;
    }
}

signed main() {
    ios::sync_with_stdio(0);
    cin.tie(0);

    cin >> n;
    for(int i = 0; i < n; i ++ ) {  // 前序
        cin >> P[i];
        P[i] --;
    }
    for(int i = 0; i < n; i ++ ) {  //中序
        cin >> I[i];
        I[i] --;
        J[I[i]] = i;
    }
    if(P[0] != 0) cout << -1 << endl;   //如果不以0为根，输出no
    else {
        fill(L, L + n, -1);
        fill(R, R + n, -1);
        bool ok = true;
        dfs(ok, 0, n, 0);
        if(!ok) {
            cout << -1 << endl;
        }
        else {
            for(int i = 0; i < n; i ++ )
                cout << L[i] + 1 << " " << R[i] + 1 << endl;
        }
    }
    return 0;
}
```



### 网络流

#### 1. 总纲

```c++
1. 基本概念
    1.1 流网络，不考虑反向边
    1.2 可行流，不考虑反向边
        1.2.1 两个条件：容量限制、流量守恒
        1.2.2 可行流的流量指从源点流出的流量 - 流入源点的流量
        1.2.3 最大流是指最大可行流
    1.3 残留网络，考虑反向边，残留网络的可行流f' + 原图的可行流f = 原题的另一个可行流
        (1) |f' + f| = |f'| + |f|
        (2) |f'| 可能是负数
    1.4 增广路径
    1.5 割
        1.5.1 割的定义
        1.5.2 割的容量，不考虑反向边，“最小割”是指容量最小的割。
        1.5.3 割的流量，考虑反向边，f(S, T) <= c(S, T)
        1.5.4 对于任意可行流f，任意割[S, T]，|f| = f(S, T)
        1.5.5 对于任意可行流f，任意割[S, T]，|f| <= c(S, T)
        1.5.6 最大流最小割定理
            (1) 可以流f是最大流
            (2) 可行流f的残留网络中不存在增广路
            (3) 存在某个割[S, T]，|f| = c(S, T)
2. 算法
    2.1 EK O(nm^2)
    2.2 Dinic O(n^2m)
3 应用
    3.1 二分图
        (1) 二分图匹配
        (2) 二分图多重匹配
    3.2 上下界网络流
        (1) 无源汇上下界可行流
        (2) 有源汇上下界最大流
        (3) 有源汇上下界最小流
    3.3 多源汇最大流
```



#### 2. 最大流

##### EK算法

时间复杂度$O(nm^2)$ 

```c++
int h[N], e[M], f[M], ne[M], idx;
int n, m, S, T;
int d[N], pre[N];
bool st[N];

void add(int a, int b, int c) {
    e[idx] = b, f[idx] = c, ne[idx] = h[a], h[a] = idx ++ ;
    e[idx] = a, f[idx] = 0, ne[idx] = h[b], h[b] = idx ++ ;
}

bool bfs() {
    queue<int> q;
    memset(st, 0, sizeof st);
    q.push(S), st[S] = true;
    d[S] = INF;
    while(!q.empty()) {
        int u = q.front(); q.pop();
        for(int i = h[u]; ~i; i = ne[i]) {
            int v = e[i];
            if(!st[v] && f[i]) {
                st[v] = true;
                d[v] = min(d[u], f[i]);
                pre[v] = i;
                if(v == T) return true;
                q.push(v);
            }
        }
    }
    return false;
}

int EK() {
    int flow = 0;
    while(bfs()) {
        flow += d[T];
        for(int i = T; i != S; i = e[pre[i] ^ 1])
            f[pre[i]] -= d[T], f[pre[i] ^ 1] += d[T];
    }
    return flow;
}
```



##### Dinic算法

时间复杂度$O(n^2m)$ 

```c++
const int N = 10010, M = 200010, INF = 1e9;

int h[N], e[M], f[M], ne[M], idx;
int n, m, S, T;
int d[N], cur[N];

void add(int a, int b, int c) {
    e[idx] = b, f[idx] = c, ne[idx] = h[a], h[a] = idx ++ ;
    e[idx] = a, f[idx] = 0, ne[idx] = h[b], h[b] = idx ++ ;
}

bool bfs() {
    queue<int> q;
    memset(d, -1, sizeof d);
    q.push(S), d[S] = 0, cur[S] = h[S];
    while(!q.empty()) {
        int u = q.front(); q.pop();
        for(int i = h[u]; ~i; i = ne[i]) {
            int v = e[i];
            if(d[v] == -1 && f[i]) {
                d[v] = d[u] + 1;
                cur[v] = h[v];
                if(v == T) return true;
                q.push(v);
            }
        }
    }
    return false;
}

int find(int u, int maxFlow) {
    if(u == T) return maxFlow;
    int Flow = 0;
    for(int i = cur[u]; ~i && Flow < maxFlow; i = ne[i]) {
        cur[u] = i;
        int v = e[i];
        if(d[v] == d[u] + 1 && f[i]) {
            int t = find(v, min(f[i], maxFlow - Flow));
            if(!t) d[v] = -1;
            f[i] -= t, f[i ^ 1] += t, Flow += t;
        }
    }
    return Flow;
}

int Dinic() {
    int maxFlow = 0, Flow = 0;
    while(bfs())
        while(Flow = find(S, INF)) maxFlow += Flow;
    return maxFlow;
}
```



##### 无源汇上下界可行流

<img src="C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20230304153950064.png" alt="image-20230304153950064" style="zoom: 67%;" />

```c++
int h[N], e[M], f[M], l[M], ne[M], idx;
int n, m, S, T;
int d[N], cur[N], D[N];

void add(int a, int b, int c, int d) {
    e[idx] = b, f[idx] = d - c, l[idx] = c, ne[idx] = h[a], h[a] = idx ++ ;
    e[idx] = a, f[idx] = 0, ne[idx] = h[b], h[b] = idx ++ ;
}

bool bfs() {
    queue<int> q;
    memset(d, -1, sizeof d);
    q.push(S), d[S] = 0, cur[S] = h[S];
    while(!q.empty()) {
        int u = q.front(); q.pop();
        for(int i = h[u]; ~i; i = ne[i]) {
            int v = e[i];
            if(d[v] == -1 && f[i]) {
                d[v] = d[u] + 1;
                cur[v] = h[v];
                if(v == T) return true;
                q.push(v);
            }
        }
    }
    return false;
}

int find(int u, int maxFlow) {
    if(u == T) return maxFlow;
    int Flow = 0;
    for(int i = cur[u]; ~i && Flow < maxFlow; i = ne[i]) {
        cur[u] = i;
        int v = e[i];
        if(d[v] == d[u] + 1 && f[i]) {
            int t = find(v, min(f[i], maxFlow - Flow));
            if(!t) d[v] = -1;
            f[i] -= t, f[i ^ 1] += t, Flow += t;
        }
    }
    return Flow;
}

int Dinic() {
    int maxFlow = 0, Flow = 0;
    while(bfs())
        while(Flow = find(S, INF)) maxFlow += Flow;
    return maxFlow;
}

signed main() {
    ios::sync_with_stdio(0);
    cin.tie(0);

    cin >> n >> m;
    memset(h, -1, sizeof h);
    S = 0, T = n + 1;
    for(int i = 0; i < m; i ++ ) {
        int a, b, c, d; cin >> a >> b >> c >> d;
        add(a, b, c, d);
        D[a] -= c, D[b] += c;
    }
    int tot = 0;
    for(int i = 1; i <= n; i ++ ) {
        if(D[i] > 0) add(S, i, 0, D[i]), tot += D[i];
        else if(D[i] < 0) add(i, T, 0, -D[i]);
    }
    int ans = Dinic();
    if(ans != tot) cout << "NO" << endl;
    else {
        cout << "YES" << endl;
        for(int i = 0; i < m * 2; i += 2) {
            cout << f[i ^ 1] + l[i] << endl;
        }
    }

    return 0;
}
```



# 动态规划

## 线性dp

1. **最长上升子序列**

```c++
vector<int>stk;//模拟堆栈
stk.push_back(w[1]);

for (int i = 2; i <= n; ++i) {
    if (w[i] > stk.back())//如果该元素大于栈顶元素,将该元素入栈
        stk.push_back(w[i]);
    else//替换掉第一个大于或者等于这个数字的那个数
        *lower_bound(stk.begin(), stk.end(), w[i]) = w[i];
}
cout << stk.size() << endl;
```

**2. 最长公共子序列LCS**

$f_{i,j}$ 表示$a$ 中前$i$ 个数字, 以及$b$ 中前$j$ 个数字中, 最长的公共序列, 那么状态转移方程：

$f_{i,j} = max(f_{i,i-1}, f_{i-1,j})\ \ [a_i != b_j]$      $f_{i,j} = max(f_{i,i}, f_{i-1,j-1} + 1)\ \ [a_i == b_j]$ 



3. **最长公共上升子序列LCIS**

$f_{i,j}$ 表示$a$ 中前$i$ 个数字, 以及$b$ 中前$j$ 个数字, 且以$b_j$ 结尾的公共上升子序列的集合, 状态转移方程为

如果不包含$a_i$ 可以写出:   $f_{i,j} = f_{i-1,j}$ 

如果包含$a_i$ (也就是$a_i = b_j$), 我们可以考虑倒数第二个数是$b$ 中的第几个数, 然后取最大值

$f_{i,j} = max(f_{i-1, k} + 1), k \in[1, j - 1]$ 

 暴力$O(n^3)$

```c++
for(int i = 1; i <= n; i ++ ) {
    for(int j = 1; j <= n; j ++ ) {
        f[i][j] = f[i - 1][j];
        if(a[i] == b[j]) {
         	int maxv = 1;
        	for(int k = 1; k < j; k ++ ) 
                if(a[i] > b[k]) maxv = max(maxv, f[i - 1][k] + 1);
        }
        f[i][j] = max(f[i][j], maxv);
    }
}
```

优化$O(n^2)$ 

```c++
for(int i = 1; i <= n; i ++ ) {
    for(int j = 1; j <= n; j ++ ) {
        f[i][j] = f[i - 1][j];
        if(a[i] == b[j]) {
            f[i][j] = max(f[i][j], maxv);
        }
        if(a[i] > b[j]) maxv = max(maxv, f[i - 1][j] + 1);
    }
}
```





# 数据结构

### 双链表

```c++
const int N = 100010;
int m;
int e[N], l[N], r[N], idx;
// 在节点a的右边插入一个数x
void insert(int a, int x) {
    e[idx] = x;
    l[idx] = a, r[idx] = r[a];
    l[r[a]] = idx, r[a] = idx++;
}
// 删除节点a
void remove(int a) {
    l[r[a]] = l[a];
    r[l[a]] = r[a];
}

int main() {
    cin >> m;
    // 0是左端点，1是右端点
    r[0] = 1, l[1] = 0;
    idx = 2;

    while (m--) {
        string op;
        cin >> op;
        int k, x;
        if (op == "L") {
            cin >> x;
            insert(0, x);
        } else if (op == "R") {
            cin >> x;
            insert(l[1], x);
        } else if (op == "D") {
            cin >> k;
            remove(k + 1);
        } else if (op == "IL") {
            cin >> k >> x;
            insert(l[k + 1], x);	//idx = k + 1 代表是第k个插入的数, l[k + 1] 代表是第k个插入的数左边的数
        } else {
            cin >> k >> x;
            insert(k + 1, x);
        }
    }

    for (int i = r[0]; i != 1; i = r[i]) cout << e[i] << ' ';
    cout << endl;

    return 0;
}

```



### 滑动窗口

```c++
int a[N], q[N];
int n, k;	//滑动窗口长度k
int main() {
    scanf("%d%d", &n, &k);
    for (int i = 0; i < n; i++) scanf("%d", &a[i]);

    int hh = 0, tt = -1;	//滑动窗口最小值
    for (int i = 0; i < n; i++) {
        if (hh <= tt && i - k + 1 > q[hh]) hh++;

        while (hh <= tt && a[q[tt]] >= a[i]) tt--;
        q[++tt] = i;

        if (i >= k - 1) printf("%d ", a[q[hh]]);
    }
    return 0;
}
```



### FenwickTree

```C++
template <typename T>
struct Fenwick {
    int n;
    std::vector<T> a;   //下标从0开始

    Fenwick(int n = 0) {
        init(n);
    }

    Fenwick(vector<T>& b) {
        init(b);
    }

    void init(int n) {
        this->n = n;
        a.assign(n, T());
    }

    void init(vector<T> &b) {   //根据b O(n) 建树
        this->n = b.size();
        a.assign(n, T());
        for (int i = 1; i <= n; i++) {
            a[i - 1] += b[i - 1];
            int j = i + (i & -i);
            if (j <= n) a[j - 1] += a[i - 1];
        }
    }

    void add(int x, T v) {
        for (int i = x + 1; i <= n; i += i & -i) {
            a[i - 1] += v;
        }
    }

    T sum(int x) {      //求[0, x)的前缀和
        auto ans = T();
        for (int i = x; i > 0; i -= i & -i) {
            ans += a[i - 1];
        }
        return ans;
    }

    T rangeSum(int l, int r) {  //求[l, r)的区间和
        return sum(r) - sum(l);
    }

    int kth(T k) {
        int x = 0;
        for (int i = 1 << std::__lg(n); i; i /= 2) {
            if (x + i <= n && k >= a[x + i - 1]) {
                x += i;
                k -= a[x - 1];
            }
        }
        return x;
    }
};
```



### 动态开点线段树

处理大值域, 或者带负数的值域, 可以动态开点, 

```c++
#define ls(x) tree[x].ls
#define rs(x) tree[x].rs
#define val(x) tree[x].val
#define mark(x) tree[x].mark
int cnt = 1;
struct node {
    int val, mark;
    int ls, rs;
} tree[N << 6];

void pushup(int p) {
    val(p) = val(ls(p)) + val(rs(p));
}

void pushdown(int p, int len) {
    if (len <= 1) return;
    if (!ls(p)) ls(p) = ++cnt;	//没有左右儿子就开一个点.
    if (!rs(p)) rs(p) = ++cnt;
    if(mark(p)) {
        val(ls(p)) += mark(p) * (len / 2);
        mark(ls(p)) += mark(p);
        val(rs(p)) += mark(p) * (len - len / 2);
        mark(rs(p)) += mark(p);
        mark(p) = 0;
    }
}

LL query(int l, int r, int p = 1, int cl = 1, int cr = 1e9) {
    if (cl >= l && cr <= r) return val(p);
    pushdown(p, cr - cl + 1);
    LL mid = (cl + cr - 1) / 2, ans = 0;
    if (mid >= l) ans += query(l, r, ls(p), cl, mid);
    if (mid < r) ans += query(l, r, rs(p), mid + 1, cr);
    return ans;
}

void update(int l, int r, int d, int p = 1, int cl = 1, int cr = 1e9) {
    if (cl >= l && cr <= r) return val(p) += d * (cr - cl + 1), mark(p) += d, void();
    pushdown(p, cr - cl + 1);
    int mid = (cl + cr - 1) / 2;
    if (mid >= l) update(l, r, d, ls(p), cl, mid);
    if (mid < r) update(l, r, d, rs(p), mid + 1, cr);
    pushup(p);
}
```



PS: 清空代码:

```c++
memset(tree, 0, sizeof (tree[0]) * (cnt + 2));
```



### 平衡树(Treap)

```C++
const int N = 100010, INF = 1e8;

int n;
struct Node {
    int l, r;
    int key, val;
    int cnt, size;
} tr[N];

int root, idx;

void pushup(int p) {
    tr[p].size = tr[tr[p].l].size + tr[tr[p].r].size + tr[p].cnt;
}

int get_node(int key) {
    tr[++idx].key = key;
    tr[idx].val = rand();
    tr[idx].cnt = tr[idx].size = 1;
    return idx;
}

void zig(int &p){ // 右旋
    int q = tr[p].l;
    tr[p].l = tr[q].r, tr[q].r = p, p = q;
    pushup(tr[p].r), pushup(p);
}

void zag(int &p){ // 左旋
    int q = tr[p].r;
    tr[p].r = tr[q].l, tr[q].l = p, p = q;
    pushup(tr[p].l), pushup(p);
}

void build() {
    get_node(-INF), get_node(INF);
    root = 1, tr[1].r = 2;
    pushup(root);

    if (tr[1].val < tr[2].val) zag(root);
}


void insert(int &p, int key) {
    if (!p) p = get_node(key);
    else if (tr[p].key == key) tr[p].cnt++;
    else if (tr[p].key > key) {
        insert(tr[p].l, key);
        if (tr[tr[p].l].val > tr[p].val) zig(p);
    } else {
        insert(tr[p].r, key);
        if (tr[tr[p].r].val > tr[p].val) zag(p);
    }
    pushup(p);
}

void remove(int &p, int key) {
    if (!p) return;
    if (tr[p].key == key) {
        if (tr[p].cnt > 1) tr[p].cnt--;
        else if (tr[p].l || tr[p].r) {
            if (!tr[p].r || tr[tr[p].l].val > tr[tr[p].r].val) {
                zig(p);
                remove(tr[p].r, key);
            } else {
                zag(p);
                remove(tr[p].l, key);
            }
        } else p = 0;
    } else if (tr[p].key > key) remove(tr[p].l, key);
    else remove(tr[p].r, key);

    pushup(p);
}

int get_rank_by_key(int p, int key){ // 通过数值找排名
    if (!p) return 0;	//如果这个数值不存在就返回0
    if (tr[p].key == key) return tr[tr[p].l].size + 1;
    if (tr[p].key > key) return get_rank_by_key(tr[p].l, key);
    return tr[tr[p].l].size + tr[p].cnt + get_rank_by_key(tr[p].r, key);
}

int get_key_by_rank(int p, int rank){   // 通过排名找数值
    if (!p) return INF;     //如果这个排名不存在就返回INF
    if (tr[tr[p].l].size >= rank) return get_key_by_rank(tr[p].l, rank);
    if (tr[tr[p].l].size + tr[p].cnt >= rank) return tr[p].key;
    return get_key_by_rank(tr[p].r, rank - tr[tr[p].l].size - tr[p].cnt);
}

int get_prev(int p, int key){ // 找到严格小于key的最大数
    if (!p) return -INF;
    if (tr[p].key >= key) return get_prev(tr[p].l, key);
    return max(tr[p].key, get_prev(tr[p].r, key));
}

int get_next(int p, int key){ // 找到严格大于key的最小数
    if (!p) return INF;
    if (tr[p].key <= key) return get_next(tr[p].r, key);
    return min(tr[p].key, get_next(tr[p].l, key));
}

```

### 笛卡尔树

笛卡尔树是一种二叉树，每一个结点由一个键值二元组$(k, w)$ 构成。要求$k$ 满足二叉搜索树的性质，而$w$ 满足堆的性质。如果笛卡尔树的$(k, w)$键值确定，且$k$ 互不相同，$w$ 互不相同，那么这个笛卡尔树的结构是唯一的。

以两个数组下标为$k$ 构建的两个笛卡尔树同构  <=>  构建笛卡尔树的两个数组的$w$ 值的大小关系相同 （$w$ 不一定相同）

<img src="C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20221118223643419.png" alt="image-20221118223643419" style="zoom:50%;" />

构建过程： 使用栈构建，

我们考虑将元素按照键值$k$ 排序。然后一个一个插入到当前的笛卡尔树中。那么每次我们插入的元素必然在这个树的右链（右链：即从根结点一直往右子树走，经过的结点形成的链）的末端。于是我们执行这样一个过程，从下往上比较右链结点与当前结点$u$的$w$，如果找到了一个右链上的结点$x$满足$x_w < u_w$，就把$u$接到$x$ 的右儿子上，而$x$ 原本的右子树就变成$u$ 的左子树。

所以我们使用一个栈来维护右链即可

<img src="C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20221118223318414.png" alt="image-20221118223318414" style="zoom: 80%;" />

```c++
int stk[N], a[N], ls[N], rs[N];
int n, m;

int build(int n) {   //返回根节点, 参数是数组的长度，下标[1, n]
    int top = 0;
    for (int i = 1; i <= n; i++) {
        int k = top;
        while (k && a[stk[k]] < a[i]) k --; // < 为构建大根堆， > 为构建小根堆
        if (k) rs[stk[k]] = i;  // rs代表笛卡尔树每个节点的右儿子
        if (k < top) ls[i] = stk[k + 1];  // ls代表笛卡尔树每个节点的左儿子
        stk[top = ++k] = i;
    }
    return stk[1];
}
```



### ZKW线段树

李超线段树是一种用于维护平面直角坐标系内线段关系的数据结构。它常被用来处理这样一种形式的问题：给定一个平面直角坐标系，支持动态插入一条线段，询问从某一个位置 $(x,+∞)$ 向下看能看到的最高的一条线段（也就是给一条竖线，问这条竖线与所有线段的最高的交点）

<img src="https://pic4.zhimg.com/v2-87357e8a2b67b615616f2ead6b63f01b_r.jpg" alt="img" style="zoom:67%;" />

### 莫队算法

1. **基础莫队**

离线询问区间，左端点分块处理，右端点在块内递增，时间复杂度为$O(n\sqrt n)$ 

```c++
struct PPP{
    int id, l, r;
};

signed main() {
    int n, q; cin >> n >> q;
    vector<int> a(n + 1), cnt(N + 1), ans(q);
    vector<PPP> query;
    for (int i = 1; i <= n; i ++ ) cin >> a[i];

    for (int i = 0; i < q; i ++ ) {
        int l, r; cin >> l >> r;
        query.push_back({i, l, r});
    }
    
    int B = sqrt(n);
    sort(query.begin(), query.end(), [&](PPP &a, PPP &b){
        if ((a.l / B) == (b.l / B)) return a.r < b.r;
        return (a.l / B) < (b.l / B);
    });

    auto add = [&](int x, int &res) {
        res += cnt[x] * (cnt[x] - 1) / 2;
        cnt[x] ++ ;
    };

    auto del = [&](int x, int &res) {
        cnt[x] --;
        res -= cnt[x] * (cnt[x] - 1) / 2;
    };

    for (int i = 0, L = 1, R = 0, res = 0; i < q; i ++ ) {
        auto [id, l, r] = query[i];
        while (R < r) add(a[++ R], res);
        while (R > r) del(a[R --], res);
        while (L > l) add(a[-- L], res);
        while (L < l) del(a[L ++ ], res);
        ans[id] = res;
    }

    for (int i = 0; i < q; i ++ ) cout << ans[i] << endl;

    return 0;
}
```

2. **带修莫队**



3. **回滚莫队**

   

4. **树上莫队**



5. **二次离线莫队**

### 树链剖分

#### 		轻重链剖分

```c++
LL lazy[N << 2], sum[N << 2];
vector<int> g[N];
int n, w[N];
int fa[N], son[N], dep[N], top[N], sz[N];//top重链顶点，son重儿子
int id[N], dfn, nw[N];	// id[i] 是节点i的dfn序编号. nw[id[i]] = w[i];

void pushup(int u) {
    sum[u] = sum[u << 1] + sum[u << 1 | 1];
}

void pushdown(int u, int L, int R) {
    if(lazy[u]) {
        int mid = L + R >> 1;
        lazy[u << 1] += lazy[u], lazy[u << 1 | 1] += lazy[u];
        sum[u << 1] += lazy[u] * (mid - L + 1);
        sum[u << 1 | 1] += lazy[u] * (R - mid);
        lazy[u] = 0;
    }
}

void build(int u, int l, int r) {
    if(l == r) {
        sum[u] = nw[l];
        return ;
    }
    int mid = l + r >> 1;
    build(u << 1, l, mid), build(u << 1 | 1, mid + 1, r);
    pushup(u);
}

void update(int u, int L, int R, int l, int r, int k) {
    if(l <= L && r >= R) {
        lazy[u] += k;
        sum[u] += 1ll * k * (R - L + 1);
        return ;
    }
    pushdown(u, L, R);
    int mid = L + R >> 1;
    if(l <= mid) update(u << 1, L, mid, l, r, k);
    if(r > mid) update(u << 1 | 1, mid + 1, R, l, r, k);
    pushup(u);
}

LL query(int u, int L, int R, int l, int r) {
    if(L >= l && R <= r) {
        return sum[u];
    }
    pushdown(u, L, R);
    int mid = L + R >> 1;
    LL ans = 0;
    if(l <= mid) ans += query(u << 1, L, mid, l, r);
    if(r > mid) ans += query(u << 1 | 1, mid + 1, R, l, r);
    return ans;
}

void dfs1(int u, int father, int depth) {
    sz[u] = 1, fa[u] = father, dep[u] = depth;
    for (auto v : g[u]) {
        if(v == father) continue;
        dfs1(v, u, depth + 1);
        sz[u] += sz[v];
        if(sz[son[u]] < sz[v]) son[u] = v;
    }
}

void dfs2(int u, int t) {
    id[u] = ++ dfn, top[u] = t, nw[dfn] = w[u];
    if(!son[u]) return ;
    dfs2(son[u], t);
    for (auto v : g[u]) {
        if(v == fa[u] || v == son[u]) continue;
        dfs2(v, v);
    }
}

signed main() {
    ios::sync_with_stdio(0);
    cin.tie(0);

    cin >> n;
    for (int i = 1; i <= n; i ++ ) cin >> w[i];
    for (int i = 1; i < n; i ++ ) {
        int a, b; cin >> a >> b;
        g[a].push_back(b);
        g[b].push_back(a);
    }
    dfs1(1, -1, 1);
    dfs2(1, 1);
    build(1, 1, n);

    int m; cin >> m;
    while (m -- ) {
        int t, u, v, k; cin >> t >> u;
        if(t == 1) {    //更新路径
            cin >> v >> k;
            while (top[u] != top[v]) {
                if(dep[top[u]] < dep[top[v]]) swap(u, v);
                update(1, 1, n, id[top[u]], id[u], k);
                u = fa[top[u]];
            }
            if(dep[u] < dep[v]) swap(u, v);
            update(1, 1, n, id[v], id[u], k);
            
        } else if (t == 2) {   //更新子树
            cin >> k;
            update(1, 1, n, id[u], id[u] + sz[u] - 1, k);
            
        } else if (t == 3) {    // 查询路径
            cin >> v;
            LL res = 0;
            while (top[u] != top[v]) {
                if(dep[top[u]] < dep[top[v]]) swap(u, v);
                res += query(1, 1, n, id[top[u]], id[u]);
                u = fa[top[u]];
            }
            if(dep[u] < dep[v]) swap(u, v);
            res += query(1, 1, n, id[v], id[u]);
            cout << res << endl;
            
        } else {    // 查询子树
            cout << query(1, 1, n, id[u], id[u] + sz[u] - 1) << endl;
        }
    }

    return 0;
}
```



树剖求LCA：

```c++
int lca(int x, int y) {
    while (top[x] != top[y]) {
        if (dep[top[x]] > dep[top[y]])
            x = fa[top[x]];
        else y = fa[top[y]];
    }
    return dep[x] < dep[y] ? x : y;
}
```







# 字符串

### 字符串hash

```c++
const int P = 131, N = 1e5 + 10;
ULL h[N], p[N];     //h[i]存储的是从开头到i的字符串hash值, p[i]存储的是P^i的值
                    //用ULL的好处:可以自动对2^64取模                   
void Hash(){
    p[0] = 1;
    for(int i = 1; i <= n; i ++ ) {
        h[i] = h[i - 1] * P + str[i];
        p[i] = p[i - 1] * P;
    }
}

ULL get(int l, int r) {
    return h[r] - h[l - 1] * p[r - l + 1];
}
```

二维hash

```C++
char c[N][N];
int n, m;
ull pow1[N], pow2[N];
ull hash1[N][N];
void init() {
    pow1[0] = pow2[0] = 1;
    for(int i = 1; i < N; i ++ ) {
        pow1[i] = pow1[i - 1] * P1;
        pow2[i] = pow2[i - 1] * P2;
    }

    for(int i = 1; i <= n; i ++ )
        for(int j = 1; j <= m; j ++ ) {
            hash1[i][j] = hash1[i - 1][j] * P1 + hash1[i][j - 1] * P2 - hash1[i - 1][j - 1] * P1 * P2 + c[i][j];
        }
}

ull get(int lx, int ly, int rx, int ry) {
    return hash1[rx][ry] - hash1[lx - 1][ry] * pow1[rx - lx + 1] - hash1[rx][ly - 1] * pow2[ry - ly + 1] + hash1[lx - 1][ly - 1] * pow1[rx - lx + 1] * pow2[ry - ly + 1];
}

```



### Trie树

```C++
int son[N][26], cnt[N], idx;

void insert(char *str){  // 插入字符串
    int p = 0;
    for (int i = 0; str[i]; i ++ ){
        int u = str[i] - 'a';
        if (!son[p][u]) son[p][u] = ++ idx;
        p = son[p][u];
    }
    cnt[p] ++ ;	//记录哪个节点处有一个字符串结束
}

int query(char *str){  // 查询字符串出现次数
    int p = 0;
    for (int i = 0; str[i]; i ++ ){
        int u = str[i] - 'a';
        if (!son[p][u]) return 0;
        p = son[p][u];
    }
    return cnt[p];
}

```



#### 01Trie树

```c++
int idx = 1;
void insert(int x){  // 插入
    int p = 1;
    for (int i = 31; i >= 0; i -- ){
        int u = x >> i & 1;
        if (!son[p][u]) son[p][u] = ++ idx;
        p = son[p][u];
        cnt[p] ++ ;
    }
}

int query(int x){  // 查询出现次数
    int p = 1, res = 0;
    for (int i = 31; i >= 0; i -- ){
        int u = x >> i & 1;
        if (son[p][!u] && cnt[son[p][!u]]) {
            res = res * 2 + 1;
            p = son[p][!u];
        }
        else {
            p = son[p][u];
            res *= 2;
        }
    }
    return res;
}

void del(int x) {
    int p = 1;
    for(int i = 31; i >= 0; i -- ) {
        int u = x >> i & 1;
        p = son[p][u];
        cnt[p] --;
    }
}
```



### KMP

时间复杂度$O(n)$

```C++
int ne[N]; //ne[i]代表的是长度为i的前缀字符串中, 最长公共前后缀的长度 一定要注意, 下标从1开始
void Get(char *p) {
    ne[0] = ne[1] = 0;
    int n = strlen(p + 1);
    for (int i = 2, j = 0; i <= n; i++) {
        while (j && p[i] != p[j + 1]) j = ne[j];
        if (p[i] == p[j + 1]) j++;
        ne[i] = j;
    }
}

void kmp(char *s, char *p) {
    int n = strlen(p + 1), m = strlen(s + 1);
    Get(p);
    for (int i = 1, j = 0; i <= m; i++) {
        while (j && s[i] != p[j + 1]) j = ne[j];
        if (s[i] == p[j + 1]) j++;
        if (j == n) {
            //return 需要的信息

            j = ne[j];
        }
    }
}

```

#### kmp + dp

例: 求出对于某个字符串的前缀字符串在该串的出现次数

方法: dp

`dp[i]` 表示以`i` 结尾的前缀的出现次数, 那么

`dp[i] = dp[ne[i]] + 1`

#### 最小循环元

如果一个字符串$S$ 是由一个字符串$T$ 重复$K$ 次形成的，则$T$ 是$S$ 的循环元。使K最大的字符串$T$ 称为$S$ 的最小循环元，此时的$K$ 称为最大循环次数。现在给一个长度为$N$ 的字符串$ S$，对$S$ 的每一个前缀，如果它的最大循环次数大于1，则输出该前缀的最小循环元长度和最大循环次数。

**结论** : 当对于每个$i$, `i - ne[i] % i == 0`时, `S[1 ~ i - ne[i]]`就是$S$的循环元, 最大循环次数为$i / (i - ne[i])$

​			$S$ 的最小循环元为$n - ne[n]$

 ```C++
for(int i = 2; i <= n; i ++ ){
    int k = i - ne[i];
    if(i % k == 0 && i / k > 1) cout << i << " " << i / k << endl;
}
 ```



### AC自动机

AC自动机建图过程

```C++
void insert(char str[])  // 将str插入Trie中
{
    int p = 0;
    for (int i = 0; str[i]; i++) {
        int u = str[i] - 'a';
        if (!tr[p][u]) tr[p][u] = ++idx;
        p = tr[p][u];
    }
    cnt[p]++;  // 记录单词出现次数
}

void build()  // 创建AC自动机
{
    int hh = 0, tt = -1;
    for (int i = 0; i < 26; i++)
        if (tr[0][i])
            q[++tt] = tr[0][i];
    while (hh <= tt) {
        int t = q[hh++];
        for (int i = 0; i < 26; i++) {
            int p = tr[t][i];
            if (!p) tr[t][i] = tr[ne[t]][i];
            else {
                ne[p] = tr[ne[t]][i];
                cnt[p] += cnt[ne[p]];
                q[++tt] = p;
            }
        }
    }
}

```



例题: 

![image-20220712201414516](C:\Users\86150\AppData\Roaming\Typora\typora-user-images\image-20220712201414516.png)

```C++
#include <bits/stdc++.h>

using namespace std;
const int N = 1e4 + 10, S = 55, M = 1e6 + 10;
int tr[N * S][26], idx, cnt[N * S];
int ne[N * S];          //next数组
int n, m;
string str;

void insert()           //正常的trie树的插入
{
    int p = 0;
    for (int i = 0; str[i]; i++) {
        int t = str[i] - 'a';
        if (!tr[p][t]) tr[p][t] = ++idx;
        p = tr[p][t];
    }
    cnt[p]++;
}

void build() {
    queue<int> q;       //用bfs来实现建立trie图的过程
    for (int i = 0; i < 26; i++)
        if (tr[0][i]) q.push(tr[0][i]);

    while (q.size()) {
        int t = q.front();
        q.pop();

        for (int i = 0; i < 26; i++) {
            int p = tr[t][i];
            if (!p) tr[t][i] = tr[ne[t]][i];
            else {
                ne[p] = tr[ne[t]][i];
                q.push(p);
            }
        }
    }
}

int main() {
    int T;
    cin >> T;
    while (T--) {
        cin >> n;
        memset(tr, 0, sizeof tr);
        memset(cnt, 0, sizeof cnt);
        memset(ne, 0, sizeof ne);
        idx = 0;

        for (int i = 0; i < n; i++) {
            cin >> str;
            insert();
        }

        build();
        cin >> str;
        int res = 0;

        for (int i = 0, j = 0; str[i]; i++) {
            int t = str[i] - 'a';
            j = tr[j][t];

            int p = j;
            while (p) {
                res += cnt[p];
                cnt[p] = 0;
                p = ne[p];
            }
        }
        printf("%d\n", res);
    }
    return 0;
}
```



# 其他

### $[0, x]$中数字$d$出现的次数

```c++
LL Count(int x, int d)  //求 0 - x 中, 数字d出现的次数
{
    LL ret = d == 0;
    for(int k, i = 1; k = x / i; i *= 10)
    {
        int high = k / 10, cur = k % 10;
        if(d == 0) { if(high) --high; else break; }
        ret += high * i;
        if(cur > d) ret += i;
        else if(cur == d) ret += x - k * i + 1;
    }
    return ret;
}
```





### MEX操作

```C++
struct MEX{
    set<int> st;
    int cot[N];             //*N是所有插入的数的范围,不是数的个数
    multiset<int> mset;

    void init(){            //*第一次初始化为nlogn,之后清空使用clear()
        memset(cot, 0, sizeof cot);
        for(int i = 0; i < N; i ++ ) st.insert(i);
    }

    void insert(int x){        //*插入操作
        if(cot[x] == 0) st.erase(x);
        cot[x] ++ ;
        mset.insert(x);
    }

    void del(int x){
        if(cot[x] == 1) st.insert(x);
        cot[x] -- ;
        mset.erase(mset.find(x));
    }

    int mex(){return *st.begin();}      //*查询操作

    int size(){return mset.size();}

    void clear(){           //*初始化
        while(mset.size()) del(*mset.begin());
    }
};

```





### 小技巧

```c++
__builtin_popcount(x)	//求x中有多少个1
__builtin_popcountll(x)	//求x中有多少个1（long long）
```



+ 给定长度为64的一个数列$a$ ，然后定义一个序列$s$ ， $s[i] = $$i$ 的二进制对应的$a$ 的异或和，求$[l, r]$的$s[i]$ 

```C++
for(int j = 0; j < 64; j ++ ) {
    if(i >> j & 1) s[i] ^= a[j];
}
```

简单求法：可以对$a$ 求前缀异或和，然后直接求

```c++
LL s = 0;
for(int i = 0; i < 64; i ++ ) {
    if(l >> i & 1) s ^= a[i];
}
for(int i = 1; i < 64; i ++ ) a[i] ^= a[i - 1];
cout << s << " ";
for(LL x = l, x < r; x ++ ) {
    s ^= a[__builtin_popcountll(x ^ (x + 1)) - 1];
    cout << s << " ";
}
```

这揭示了一个性质，其实每两个相邻的数之间都可以通过异或一个形如$2^x - 1$ 的数得到， 其中$x$ 的大小正是`i^(i + 1)` 中1的个数



+ **区间加等差数列**

```c++
signed main() {
	int n, m; cin >> n >> m
    vector<int> a(n), b(n);
    // m次操作，在[0, n - 1] 上对区间[l, r] 加0 ~ r - l;
    while(m -- ) {
        int l, r; cin >> l >> r;
        //首项
        b[l] -= l;
        b[r + 1] += l;
        //公差
        a[l] ++;
        a[r + 1] --;
    }

    for(int i = 0; i < n; i ++ ) {
        a[i + 1] += a[i];
        b[i + 1] += b[i];
    }
    for(int i = 0; i < n; i ++ ) {
        cout << b[i] + i * a[i] << " \n"[i == n - 1];
    }
    
    return 0;
}
```



+ 一种交互题:

给定一个排列长度为$n$ , 一开始从编号1开始, 每一次询问可以向前走或者向后走$k$ , 然后获得当前的排列的值, 问排列的长度(排列成环)

解法: 有一种$2 * \sqrt{n}$ 次的解法, 本质的解法都是寻找相同的数, 因为是一个排列, 所以我们找到两个相同的数, 两者距离相减就是长度

所以我们先获得前$\sqrt{n}$ 的值, 如果有相同就解决了, 否则就继续每次跳$\sqrt{n}$ , 直到获得答案, 可以发现这样最多进行$2 * \sqrt{n}$ 次就可以获得





### Java 手册

1. **读入输出问题 与 BigInteger/BigDecimal**

```java
import java.io.*;
import java.util.*;

public class Main {

    static Scanner in;
    static PrintWriter out;

    public static void main(String[] args) {
        in = new Scanner(new BufferedInputStream(System.in));
        out = new PrintWriter(System.out);
        new Main().solve();
        out.flush();
    }
    
    public static void solve() {
        // ---------------读入输出问题------------------
        int n = in.nextInt(); 	// 读入整数
        String s = in.nextLine();	//读入一整行
        
		int[] a = {1, 2, 4};		//一行带空格输出
        for(int i = 0; i < 3; i ++ ) {
            System.out.print(new Integer(a[i]).toString() + " ");
        }
        
        // 输出浮点数保留位数问题，#意思为除了0以外的数显示，0的话不显示
        DecimalFormat f = new DecimalFormat("#.00#");
        DecimalFormat g = new DecimalFormat("0.000");
        double a = 132.123213, b = 0.2342312;
        System.out.println(f.format(a));
        System.out.println(g.format(a));
        System.out.println(f.format(b));
        System.out.println(g.format(b));
        
        // 输出：132.123， 132.123， .234， 0.234 
        
        
        
        // ---------------BigInteger的应用-------------
        // 1. 转换为大整数
        int a = 3;
        BigInteger A = BigInteger.valueOf(a);
        String s = "12345"
        BigInteger S = BigInteger.valueOf(s);
        BigInteger a = new BigInteger("101", 2); // 以二进制转化为大整数（十进制）
        
        // 2.基本运算
        BigInteger a = new BigInteger("23");
        BigInteger b = new BigInteger("45");
        int c = 4;
        a = a.add(b); // a += b
        a = a.subtract(b); // a -= b
        a = a.divide(b) 	// a = a / b(取整)
        a = a.multiply(b) 	// a = a * b;
        a = a.pow(c)		// a = a ^ c;
        a = a.abs()			// abs(a)
        a = a.negate();		// a = -a
        a = a.mod(b)		// a %= b;
        a = a.max(b)		// a = max(a, b)
        a.equals(b / c)	 	// a == b or c
        
        // 一些常量
        A = BigInteger.ONE    // 1  
		B = BigInteger.TEN    // 10  
		C = BigInteger.ZERO   // 0 
        
        // 读入问题
        Scanner in = new Scanner(System.in);
        int n;
        BigInteger m;
        n = in.nextInt();
        m = in.BigInteger();
        
        // BigDecimal大浮点数，类似           
    }
}
```

2. String 类

```java
public class Main {
    public static void main(String[] args) {
        //----------------String------------------
        String a = "Hello";
        a.substring(0, 4); 	//截取子串

        // java中String不可修改，只能转化为字符数组
        char[] ch = a.toCharArray();

    }
}

```





