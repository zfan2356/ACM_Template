```go
// LazySegmentTree 懒标记线段树————————————————————————————
type LazySegmentTree struct {
    n int
    info []*Info
    tag  []*Tag
}

// Init 初始化
func (lst *LazySegmentTree) Init(_init *[]*Info) {
    lst.n = len(*(_init))

    M := bits.Len(uint(lst.n)) - 1
    lst.info = make([]*Info, 4 << M)
    lst.tag = make([]*Tag, 4 << M)

    for i := 0; i < (4 << M); i++ {
        lst.info[i] = newInfo()
        lst.tag[i] = newTag()
    }

    var build func(int, int, int)
    build = func(p int, l int, r int) {
        if r - l == 1 {
            lst.info[p] = (*_init)[l]
            return
        }
        m := (l + r) / 2
        build(2 * p, l, m)
        build(2 * p + 1, m, r)
    }
    build(1, 0, lst.n)
}

func (lst *LazySegmentTree) pull(p int)  {
    lst.info[p] = lst.info[2 * p].merge(lst.info[2 * p + 1])
}

func (lst *LazySegmentTree) apply(p int, t *Tag)  {
    lst.info[p].apply(t)
    lst.tag[p].apply(t)
}

func (lst *LazySegmentTree) push(p int)  {
    lst.apply(2 * p, lst.tag[p])
    lst.apply(2 * p + 1, lst.tag[p])
    lst.tag[p] = newTag()
}

// 单点修改
func (lst *LazySegmentTree) modify(p int, l int, r int, x int, v *Info)  {
    if r - l == 1 {
        lst.info[p] = v
        return
    }
    m := (l + r) / 2
    lst.push(p)
    if x < m {
        lst.modify(2 * p, l, m, x, v)
    } else {
        lst.modify(2 * p + 1, m, r, x, v)
    }
    lst.pull(p)
}
func (lst *LazySegmentTree) Modify(x int, v *Info)  {
    lst.modify(1, 0, lst.n, x, v)
}

// 区间查询
func (lst *LazySegmentTree) rangeQuery(p int, l int, r int, x int, y int) *Info {
    if l >= y || r <= x {
        return newInfo()
    }
    if l >= x && r <= y {
        return lst.info[p]
    }
    m := (l + r) / 2
    lst.push(p)
    return lst.rangeQuery(2 * p, l, m, x, y).merge(lst.rangeQuery(2 * p + 1, m, r, x, y))
}
func (lst *LazySegmentTree) RangeQuery(l int, r int) *Info {
    return lst.rangeQuery(1, 0, lst.n, l, r)
}

// 区间修改
func (lst *LazySegmentTree) rangeApply(p int, l int, r int, x int, y int, v *Tag)  {
    if l >= y || r <= x {
        return
    }
    if l >= x && r <= y {
        lst.apply(p, v)
        return
    }
    m := (l + r) / 2
    lst.push(p)
    lst.rangeApply(2 * p, l, m, x, y, v)
    lst.rangeApply(2 * p + 1, m, r, x, y, v)
    lst.pull(p)
}
func (lst *LazySegmentTree) RangeApply(l int, r int, v *Tag)  {
    lst.rangeApply(1, 0, lst.n, l, r, v)
}

// 线段树上二分
type F func(*Info)bool
func (lst *LazySegmentTree) findFirst(p int, l int, r int, x int, y int, pred F) int {
    if l >= y || r <= x || !pred(lst.info[p]) {
        return -1
    }
    if r - l == 1 {
        return l
    }
    m := (l + r) / 2
    lst.push(p)
    res := lst.findFirst(2 * p, l, m, x, y, pred)
    if res == -1 {
        res = lst.findFirst(2 * p + 1, m, r, x, y, pred)
    }
    return res
}
func (lst *LazySegmentTree) FindFirst(l int, r int, pred F) int {
    return lst.findFirst(1, 0, lst.n, l, r, pred)
}
func (lst *LazySegmentTree) findLast(p int, l int, r int, x int, y int, pred F) int {
    if l >= y || r <= x || !pred(lst.info[p]) {
        return -1
    }
    if r - l == 1 {
        return l
    }
    m := (l + r) / 2
    lst.push(p)
    res := lst.findLast(2 * p + 1, m, r, x, y, pred)
    if res == -1 {
        res = lst.findLast(2 * p, l, m, x, y, pred)
    }
    return res
}
func (lst *LazySegmentTree) FindLast(l int, r int, pred F) int {
    return lst.findLast(1, 0, lst.n, l, r, pred)
}

type Info struct {

}
func newInfo() *Info {
    return &Info {
    }
}

type Tag struct {
}
func newTag() *Tag {
    return &Tag {
    }
}

func (a *Info) merge(b *Info) *Info {
    return &Info{
    }
}

func (a *Info) apply(v *Tag)  {
}

func (t *Tag) apply(v *Tag)  {
}
```

### 注意
1. 因为根据地址索引较为耗费时间, 所以封装起来之后的效率并不理想, 写题时尽量手写, 也可能是因为封装细节出现问题
