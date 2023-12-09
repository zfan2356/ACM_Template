### 线性基

```go
// Basis ————————————————————————————
type Basis struct {
    a [20]int
    t [20]int
}

func (b *Basis) Init()  {
    for i := 0; i < 20; i++ {
        b.t[i] = -1
    }
}

func (b *Basis) insert(x int, y int)  {
    for i := 0; i < 20; i++ {
        if (x >> i & 1) > 0 {
            if y > b.t[i] {
                b.a[i], x = x, b.a[i]
                b.t[i], y = y, b.t[i]
            }
            x ^= b.a[i]
        }
    }
}

func (b *Basis) isExist(x int, y int) bool {
    for i := 0; i < 20; i++ {
        if (x >> i & 1) > 0 && b.t[i] >= y {
            x ^= b.a[i]
        }
    }
    if x == 0 {
        return true
    }
    return false
}
```

### 注意
实际应用中也较为灵活, 没有必要存
这里的板子额外维护了一个`t`数组, 代表对应基底的优先级, 也就是我们有多个地位等同的基底
的时候, 要按照优先级选择
`isExist(x, y)` 的作用是询问`x` 是否可以被表示, 然后`y`代表其优先级下限, 也就是我们
只能选取优先级大于等于`y`的 基底