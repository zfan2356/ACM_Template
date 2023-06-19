```c++
const double EPS = 0.0001;
const double INF = 1000000;
struct point{
  double x, y;
  point(){
  }
  point(double x, double y): x(x), y(y){
  }
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
};
point rotate90(point P){
  return point(P.y, -P.x);
}
point midpoint(point P, point Q){
  return (P + Q) / 2;
}
double abs(point P){
  return sqrt(P.x * P.x + P.y * P.y);
}
double dist(point P, point Q){
  return abs(Q - P);
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
  line(point A, point B): A(A), B(B){
  }
};
point vec(line L){
  return L.B - L.A;
}
line perpendicular_bisector(point P, point Q){
  point A = midpoint(P, Q);
  point B = A + rotate90(Q - P);
  return line(A, B);
}
point line_intersection(line L1, line L2){
  return L1.A + vec(L1) * cross(L2.A - L1.A, vec(L2)) / cross(vec(L1), vec(L2));
}
point circumcenter(point A, point B, point C){
  return line_intersection(perpendicular_bisector(A, B), perpendicular_bisector(A, C));
}
```