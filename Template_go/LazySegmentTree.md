```go
// LazySegmentTree 懒标记线段树————————————————————————————
type LazySegmentTree struct {
    n int
    info []*Info
    tag  []*Tag
}

// Init 初始化
func (this *LazySegmentTree) Init(_init *[]*Info) {
    this.n = len(*(_init))

    M := bits.Len(uint(this.n)) - 1
    this.info = make([]*Info, 4 << M)
    this.tag = make([]*Tag, 4 << M)

    for i := 0; i < (4 << M); i++ {
        this.info[i] = newInfo()
        this.tag[i] = newTag()
    }

    var build func(int, int, int)
    build = func(p int, l int, r int) {
        if r - l == 1 {
            this.info[p] = (*_init)[l]
            return
        }
        m := (l + r) / 2
        build(2 * p, l, m)
        build(2 * p + 1, m, r)
    }
    build(1, 0, this.n)
}

func (this *LazySegmentTree) pull(p int)  {
    this.info[p] = this.info[2 * p].merge(this.info[2 * p + 1])
}

func (this *LazySegmentTree) apply(p int, t *Tag)  {
    this.info[p].apply(t)
    this.tag[p].apply(t)
}

func (this *LazySegmentTree) push(p int)  {
    this.apply(2 * p, this.tag[p])
    this.apply(2 * p + 1, this.tag[p])
    this.tag[p] = newTag()
}

// 单点修改
func (this *LazySegmentTree) modify(p int, l int, r int, x int, v *Info)  {
    if r - l == 1 {
        this.info[p] = v
        return
    }
    m := (l + r) / 2
    this.push(p)
    if x < m {
        this.modify(2 * p, l, m, x, v)
    } else {
        this.modify(2 * p + 1, m, r, x, v)
    }
    this.pull(p)
}
func (this *LazySegmentTree) Modify(x int, v *Info)  {
    this.modify(1, 0, this.n, x, v)
}

// 区间查询
func (this *LazySegmentTree) rangeQuery(p int, l int, r int, x int, y int) *Info {
    if l >= y || r <= x {
        return newInfo()
    }
    if l >= x && r <= y {
        return this.info[p]
    }
    m := (l + r) / 2
    this.push(p)
    return this.rangeQuery(2 * p, l, m, x, y).merge(this.rangeQuery(2 * p + 1, m, r, x, y))
}
func (this *LazySegmentTree) RangeQuery(l int, r int) *Info {
    return this.rangeQuery(1, 0, this.n, l, r)
}

// 区间修改
func (this *LazySegmentTree) rangeApply(p int, l int, r int, x int, y int, v *Tag)  {
    if l >= y || r <= x {
        return
    }
    if l >= x && r <= y {
        this.apply(p, v)
        return
    }
    m := (l + r) / 2
    this.push(p)
    this.rangeApply(2 * p, l, m, x, y, v)
    this.rangeApply(2 * p + 1, m, r, x, y, v)
    this.pull(p)
}
func (this *LazySegmentTree) RangeApply(l int, r int, v *Tag)  {
    this.rangeApply(1, 0, this.n, l, r, v)
}

// 线段树上二分
type F func(*Info)bool
func (this *LazySegmentTree) findFirst(p int, l int, r int, x int, y int, pred F) int {
    if l >= y || r <= x || !pred(this.info[p]) {
        return -1
    }
    if r - l == 1 {
        return l
    }
    m := (l + r) / 2
    this.push(p)
    res := this.findFirst(2 * p, l, m, x, y, pred)
    if res == -1 {
        res = this.findFirst(2 * p + 1, m, r, x, y, pred)
    }
    return res
}
func (this *LazySegmentTree) FindFirst(l int, r int, pred F) int {
    return this.findFirst(1, 0, this.n, l, r, pred)
}
func (this *LazySegmentTree) findLast(p int, l int, r int, x int, y int, pred F) int {
    if l >= y || r <= x || !pred(this.info[p]) {
        return -1
    }
    if r - l == 1 {
        return l
    }
    m := (l + r) / 2
    this.push(p)
    res := this.findLast(2 * p + 1, m, r, x, y, pred)
    if res == -1 {
        res = this.findLast(2 * p, l, m, x, y, pred)
    }
    return res
}
func (this *LazySegmentTree) FindLast(l int, r int, pred F) int {
    return this.findLast(1, 0, this.n, l, r, pred)
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
