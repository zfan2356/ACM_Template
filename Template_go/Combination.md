```go
// MInt ————————————————————————————
type MInt int

func getMod() MInt { return MInt(Mod) }

func (m MInt) norm() {
    if m < 0 {
        m += getMod()
    }
    if m >= getMod() {
        m -= getMod()
    }
}

func (m MInt) neg() MInt {
    return getMod() - m
}

func (m MInt) inv() MInt {
    if m == 0 {
        panic("division by zero")
    }
    return power(m, int64(getMod()-2))
}

func (m MInt) mul(rhs MInt) MInt {
    return (m * rhs) % getMod()
}

func (m MInt) add(rhs MInt) MInt {
    return (m + rhs) % getMod()
}

func (m MInt) sub(rhs MInt) MInt {
    return (m - rhs + getMod()) % getMod()
}

func (m MInt) div(rhs MInt) MInt {
    return m.mul(rhs.inv())
}

func power(a MInt, b int64) MInt {
    res := MInt(1)
    for ; b > 0; b /= 2 {
        if b%2 == 1 {
            res = res.mul(a)
        }
        a = a.mul(a)
    }
    return res
}

type Z = MInt

// Comb ————————————————————————————
type Comb struct {
    n       int
    _fac    []Z
    _invfac []Z
    _inv    []Z
}

func NewComb() *Comb {
    return &Comb{
        n:       0,
        _fac:    []Z{Z(1)},
        _invfac: []Z{Z(1)},
        _inv:    []Z{Z(0)},
    }
}

func (c *Comb) resize(m int) {
    if m <= c.n {
        return
    }
    c._fac = append(c._fac, make([]Z, m-c.n)...)
    c._invfac = append(c._invfac, make([]Z, m-c.n)...)
    c._inv = append(c._inv, make([]Z, m-c.n)...)

    for i := c.n + 1; i <= m; i++ {
        c._fac[i] = c._fac[i-1].mul(Z(i))
    }

    c._invfac[m] = c._fac[m].inv()
    for i := m; i > c.n; i-- {
        c._invfac[i-1] = c._invfac[i].mul(Z(i))
        c._inv[i] = c._invfac[i].mul(c._fac[i-1])
    }
    c.n = m
}

func (c *Comb) fac(m int) Z {
    if m > c.n {
        c.resize(2 * m)
    }
    return c._fac[m]
}

func (c *Comb) invfac(m int) Z {
    if m > c.n {
        c.resize(2 * m)
    }
    return c._invfac[m]
}

func (c *Comb) inv(m int) Z {
    if m > c.n {
        c.resize(2 * m)
    }
    return c._inv[m]
}

func (c *Comb) binom(n int, m int) Z {
    if n < m || m < 0 {
        return Z(0)
    }
    return c.fac(n).mul(c.invfac(m)).mul(c.invfac(n - m))
}

func (c *Comb) perm(n int, m int) Z {
    if n < m || m < 0 {
        return Z(0)
    }
    return c.fac(n).mul(c.invfac(n - m))
}

var Mod int = 1000000007
var comb *Comb = NewComb()


```

#### 说明
本段代码实现了自动取模类Z, 并且可以计算组合数, 适用于不断地对一个`int`类型取模

###### 设计点
1. 这里的`MInt`采用的是对`int`重命名的方式, 来实现取模类, `getMod()`返回`Mod`值, 不采取结构体方式, 这样输入输出更加自然
2. `Comb`类实现的是组合数类, 这里采用了自动扩容的方式, 如果调用的范围超过当前范围, 就会自动扩容两倍, 同时采用结构体指针的方式, 缩小内存压力, 增加运行速度


###### 尚未解决的难点
1. `MInt`不能与`int`变量运算, go语言少了方法重载, 所以只能将`int`强转为`MInt`, 不但麻烦, 而且比较担心运行效率
2. `MInt`如果利用泛型会更好一些, 也就是模数可变