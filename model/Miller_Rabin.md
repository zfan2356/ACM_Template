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