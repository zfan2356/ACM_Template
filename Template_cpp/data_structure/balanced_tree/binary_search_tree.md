## 二叉搜索树

struct Node {
int key;
Node* l;
Node* r;
// 维护其他信息
int size, cnt;
Node(int value) : key(value), size(1), cnt(1), l(nullptr), r(nullptr) {}
};

// 利用子节点更新父亲节点的信息
void pull(Node* p) {
p->size = p->cnt + (p->l ? p->l->size : 0) + (p->r ? p->r->size : 0);
}

// BST的中序遍历为非降序列
void inorderTraversal(Node* p) {
if (p == nullptr) {
return ;
}
inorderTraversal(p->l);
std::cout << p->key << ' ';
inorderTraversal(p->r);
}

// 查找最小最大值, 时间复杂度为O(h)
int findMin(Node* p) {
if (p == nullptr) {
return -1;
}
while (p->l != nullptr) {
p = p->l;
}
return p->key;
}
Node* findMinNode(Node* p) {
if (p == nullptr) {
return p;
}
while (p->l != nullptr) {
p = p->l;
}
return p;
}
int findMax(Node* p) {
if (p == nullptr) {
return -1;
}
while (p->r != nullptr) {
p = p->r;
}
return p->key;
}
Node* findMaxNode(Node* p) {
if (p == nullptr) {
return p;
}
while (p->r!= nullptr) {
p = p->r;
}
return p;
}

// 查找key值为val的元素
bool search(Node* p, int val) {
if (p == nullptr) {
return false;
}
if (p->key == val) {
return true;
} else if (p->key < val) {
return search(p->r, val);
} else {
return search(p->l, val);
}
}

// 插入一个元素
Node* insert(Node* p, int val) {
if (p == nullptr) {
return new Node(val);
}
if (val < p->key) {
p->l = insert(p->l, val);
} else if (val > p->key) {
p->r = insert(p->r, val);
} else {
p->cnt++;
}
pull(p);
return p;
}

// 删除一个元素
Node* remove(Node* p, int val) {
if (p == nullptr) {
return p;
}
if (val < p->key) {
p->l = remove(p->l, val);
} else if (val > p->key) {
p->r = remove(p->r, val);
} else {
if (p->cnt > 1) {
p->cnt--;
} else {
if (p->l == nullptr) {
auto tmp = p->r;
delete p;
return tmp;
} else if (p->r == nullptr) {
auto tmp = p->l;
delete p;
return tmp;
} else {
// 后继
auto suc = findMinNode(p->r);
p->key = suc->key;
p->cnt = suc->cnt;
p->r = remove(p->r, suc->key);
}
}
}
return p;
}

// 求元素的排名
int rank(Node* p, int val) {
if (p == nullptr) {
return 0;
}
if (p->key == val) {
return (p->l ? p->l->size : 0) + 1;
} else if (p->key > val) {
return rank(p->l, val);
} else {
return rank(p->r, val) + (p->l ? p->l->size : 0) + p->cnt;
}
}

// 查找排名为k的元素
int kth(Node* p, int k) {
if (p == nullptr) {
return -1;
}
int size_l = p->l ? p->l->size : 0;
if (size_l >= k) {
return kth(p->l, k);
} else if (size_l < k && k <= size_l + p->cnt) {
return p->key;
} else {
return kth(p->r, k - size_l - p->cnt);
}
}
