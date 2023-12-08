```go
type MInt struct {
	x int
}

func NewMInt(x int) *MInt {
	return &MInt{x: x % getMod()}
}

func (m *MInt) norm() {
	if m.x < 0 {
		m.x += getMod()
	}
	if m.x >= getMod() {
		m.x -= getMod()
	}
}

func (m *MInt) val() int {
	return m.x
}

func (m *MInt) Neg() *MInt {
	res := &MInt{
		x: getMod() - m.x,
	}
	return res
}

func (m *MInt) Inv() *MInt {
	if m.x == 0 {
		panic("division by zero")
	}
	a := &MInt{m.x}
	return power(a, int64(getMod()-2))
}

func (m *MInt) Mul(rhs *MInt) *MInt {
	return NewMInt((m.x * rhs.x) % getMod())
}

func (m *MInt) Add(rhs *MInt) *MInt {
	return NewMInt((m.x + rhs.x) % getMod())
}

func (m *MInt) Sub(rhs *MInt) *MInt {
	return NewMInt((m.x - rhs.x + getMod()) % getMod())
}

func (m *MInt) Div(rhs *MInt) *MInt {
	return m.Mul(rhs.Inv())
}

func power(a *MInt, b int64) *MInt {
	res := NewMInt(1)
	for ; b > 0; b /= 2 {
		if b%2 == 1 {
			res = res.Mul(a)
		}
		a = a.Mul(a)
	}
	return res
}

func getMod() int {
	return Mod
}

type Z = MInt

type Comb struct {
	n       int
	_fac    []Z
	_invfac []Z
	_inv    []Z
}

func NewComb() *Comb {
	return &Comb{
		n:       0,
		_fac:    []Z{Z{1}},
		_invfac: []Z{Z{1}},
		_inv:    []Z{Z{0}},
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
		c._fac[i] = *(c._fac[i-1].Mul(&Z{i}))
	}

	c._invfac[m] = *(c._fac[m].Inv())
	for i := m; i > c.n; i-- {
		c._invfac[i-1] = *(c._invfac[i].Mul(&Z{i}))
		c._inv[i] = *(c._invfac[i].Mul(&c._fac[i-1]))
	}
	c.n = m
}

func (c *Comb) fac(m int) *Z {
	if m > c.n {
		c.resize(2 * m)
	}
	return &c._fac[m]
}

func (c *Comb) invfac(m int) *Z {
	if m > c.n {
		c.resize(2 * m)
	}
	return &c._invfac[m]
}

func (c *Comb) inv(m int) *Z {
	if m > c.n {
		c.resize(2 * m)
	}
	return &c._inv[m]
}

func (c *Comb) binom(n int, m int) *Z {
	if n < m || m < 0 {
		return NewMInt(0)
	}
	return c.fac(n).Mul(c.invfac(m)).Mul(c.invfac(n - m))
}

func (c *Comb) perm(n int, m int) *Z {
	if n < m || m < 0 {
		return NewMInt(0)
	}
	return c.fac(n).Mul(c.invfac(n - m))
}

var Mod int = 1000000007

```