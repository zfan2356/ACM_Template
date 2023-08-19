```c++
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

01trie

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

int query(int x){  // 查询当前的数与x异或最大值是多少
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