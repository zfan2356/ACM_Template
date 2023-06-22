```c++
const double EPS = 1e-6;
struct point{
    double x, y;
    point(){}
    point(double x, double y): x(x), y(y){}

    point operator +(point P){
        return point(x + P.x, y + P.y);
    }
    point operator -(point P){
        return point(x - P.x, y - P.y);
    }
    point operator *(double k){
        return point(x * k, y * k);
    }
    point operator /(double k){
        return point(x / k, y / k);
    }

    double abs() {
        return sqrt(x * x + y * y);
    }
};

point rotate90(point P){
    return point(P.y, -P.x);
}

point midpoint(point P, point Q){
    return (P + Q) / 2;
}

double dist(point P, point Q){
    return (Q - P).abs();
}

bool cmp(double x, double y) {
    return fabs(x - y) < EPS;
}

double dot(point P, point Q){
    return P.x * Q.x + P.y * Q.y;
}

double cross(point P, point Q){
    return P.x * Q.y - P.y * Q.x;
}

double is_collinear(point P, point Q, point R){
    return abs(cross(Q - P, R - P)) < EPS;
}

struct line{
    point A, B;
    line(point A, point B): A(A), B(B){}
};

point vec(line L){ 
    return L.B - L.A;
}

line perpendicular_bisector(point P, point Q){  //垂直平分线
    point A = midpoint(P, Q);
    point B = A + rotate90(Q - P);
    return line(A, B);
}

point line_intersection(line L1, line L2){  //两条线段交点
    return L1.A + vec(L1) * cross(L2.A - L1.A, vec(L2)) / cross(vec(L1), vec(L2));
}

point circumcenter(point A, point B, point C){  //外接圆
    return line_intersection(perpendicular_bisector(A, B), perpendicular_bisector(A, C));
}
```

### 随机增量法求最小圆覆盖
线性时间复杂度
```c++
void solve() {
    int n;
    cin >> n;
    vector<point> p(n);
    for (int i = 0; i < n; i++) {
        cin >> p[i].x >> p[i].y;
    }

    for (int i = 0; i < n; i++) {
        swap(p[rand() % n], p[rand() % n]);
    }

    point o = p[0];
    double r = 0;
    for (int i = 0; i < n; i++) {
        if (dist(o, p[i]) < r || cmp(dist(o, p[i]), r)) {
            continue;
        }
        o = midpoint(p[i], p[0]);
        r = dist(p[i], p[0]) / 2;
        for (int j = 1; j < i; j++) {
            if (dist(o, p[j]) < r || cmp(dist(o, p[j]), r)) {
                continue;
            }
            o = midpoint(p[i], p[j]);
            r = dist(p[i], p[j]) / 2;
            for (int k = 0; k < j; k++) {
                if (dist(o, p[k]) < r || cmp(dist(o, p[k]), r)) {
                    continue;
                }
                o = circumcenter(p[i], p[j], p[k]);
                r = dist(o, p[i]);
            }
        }
    }
    cout << fixed << setprecision(10) << r << '\n';
}
```