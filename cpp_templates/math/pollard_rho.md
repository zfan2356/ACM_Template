```c++
mt19937 mt(time(0)); //随机化
i64 PR(i64 n) {
    i64 x = uniform_int_distribution<LL>(0, n - 1)(mt);
    i64 s, t, c = uniform_int_distribution<LL>(1, n - 1)(mt); //随机化
    for (int gol = 1; 1; gol <<= 1, s = t, x = 1) {
        for (int stp = 1; stp <= gol; ++stp) {
            t = (Miller_Rabin::mul(t, t, n) + c) % n;
            x = Miller_Rabin::mul(x, abs(s - t), n);
            if ((stp & 127) == 0) {
                i64 d = std::gcd(x, n);
                if (d > 1) {
                    return d;
                }
            }
        }
        i64 d = std::gcd(x, n);
        if (d > 1) {
            return d;
        }
    }
}
```
