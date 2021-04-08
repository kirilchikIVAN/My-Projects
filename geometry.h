#include <iostream>
#include <vector>
#include <cmath>


const double eps = 1e-1;
const double PI = std::acos(-1.0);


double abs(const double& d) {
    if (d > 0) return d;
    else return -d;
}


class Point {
public:
    double x = 0, y = 0;
public:
    Point() {}
    Point(const double& x_, const double& y_): x(x_), y(y_) {}
    bool operator==(const Point& p) const;
    bool operator!=(const Point& p) const {return !(*this == p);}
    double dist(const Point& p) const;
    Point operator+(const Point& p) const {return {x + p.x, y + p.y};}
    Point operator-(const Point& p) const {return {x - p.x, y - p.y};}
    void rotate(const double& angle);
    Point operator*(const double& d) const {return {x * d, y * d};}
    Point operator/(const double& d) const {return {x / d, y / d};}
    friend std::ostream& operator<<(std::ostream& out, const Point& p);
    double operator*(const Point& p) const {return x * p.x + y * p.y;}
    double operator^(const Point& p) const {return x * p.y - p.x * y;}
    void rotate(const Point& center, double angle);
    Point reflex(const Point& center);
    Point scale(Point center, double coefficient);
    double size() const { return std::sqrt(x * x + y * y); }
};


Point Point::scale(Point center, double coefficient) {
    Point tmp = *this - center;
    tmp = tmp * coefficient;
    *this = center + tmp;
    return *this;
}


Point Point::reflex(const Point& center) {
    Point tmp = center - *this;
    *this = *this + tmp * 2;
    return *this;
}


void Point::rotate(const Point& center, double angle) {
    Point tmp = *this - center;
    tmp.rotate(angle);
    *this = center + tmp;
}


std::ostream& operator<<(std::ostream& out, const Point& p) {
    out << p.x << " " << p.y;
    return out;
}


void Point::rotate(const double& angle) {
    double x_ = x;
    x = x * std::cos(angle) - y * std::sin(angle);
    y = x_ * std::sin(angle) + y * std::cos(angle);
}


double Point::dist(const Point& p) const {
    return std::sqrt((x - p.x) * (x - p.x) + (y - p.y) * (y - p.y));
}


bool Point::operator==(const Point& p) const {
    return (std::max(x - p.x, p.x - x) < eps) && (std::max(y - p.y, p.y - y) < eps);
}


class Line {
public:
    double a = 0, b = 0, c = 0;
public:
    Line() {}
    Line(const Point& f, const Point& s);
    Line(const double& k, const double& bb): a(k), b(-1), c(bb) {}
    Line(const Point& p, const double& k): a(k), b(-1), c(-p.x + p.y) {}
    Line(const Point& p, const Point& n, int): a(n.y), b(-n.x), c(n ^ p) {}
    bool operator==(const Line& l) const;
    bool operator!=(const Line& l) const {return !(*this == l);}
    Point intersec(const Line& l) const;
    double dist(const Point& p) const {return abs(a * p.x + b * p.y + c) / std::sqrt(a * a + b * b);}
    bool contains(const Point& p) const;
    Point reflex(const Point& p) const;
};


Point Line::reflex(const Point& p) const {
    if (contains(p)) return p;
    Line perp(p, Point(a, b), 1);
    Point H = intersec(perp);
    Point pH = H - p;
    return p + pH * 2;
}


bool Line::contains(const Point& p) const {
    double tmp = a * p.x + b * p.y + c;
    if (abs(tmp) < eps) return 1;
    else return false;
}


Point Line::intersec(const Line& l) const {
    double x_ = (b * l.c - l.b * c) / (a * l.b - l.a * b);
    double y_ = (a * l.c - l.a * c) / (a * l.b - l.a * b);
    return {x_, y_};
}


bool Line::operator==(const Line& l) const {
    if ((a == 0) ^ (l.a == 0)) return false;
    if ((b == 0) ^ (l.b == 0)) return false;
    if ((c == 0) ^ (l.c == 0)) return false;
    if (a == 0 && l.a == 0) {
        if (c == 0 && l.c == 0) return true;
        if (b / l.b - c / l.c < eps) return true;
        return false;
    }
    if (b == 0 && l.b == 0) {
        if (c == 0 && l.c == 0) return true;
        if (a / l.a - c / l.c < eps) return true;
        return false;
    }
    if (b / l.b - c / l.c < eps && a / l.a - c / l.c < eps) return true;
    return false;
}


Line::Line(const Point& f, const Point& s) {
    if (f.y == s.y) {
        a = 0;
        b = 1;
        c = -f.y;
        return;
    }
    a = 1;
    b = -(s.x - f.x) / (s.y - f.y);
    c = -b * f.y - f.x;
}


class Shape {
public:
    virtual double perimeter() const = 0;
    virtual double area() const = 0;
    virtual bool operator!=(const Shape& another) const = 0;
    virtual bool operator==(const Shape& another) const = 0;
    virtual bool isCongruentTo(const Shape& another) const = 0;
    virtual bool isSimilarTo(const Shape& another) const = 0;
    virtual bool containsPoint(Point point) const = 0;
    virtual void rotate(const Point& center, const double& angle) = 0;
    virtual void reflex(const Point& center) = 0;
    virtual void reflex(const Line& axis) = 0;
    virtual void scale(const Point& center, const double& coefficient) = 0;
    virtual ~Shape() = default;
};


class Polygon: public Shape {
protected:
    std::vector<Point> vertices = {};
public:
    Polygon() {}
    Polygon (const std::vector<Point>& vertices_): vertices(vertices_) {}
    Polygon (const std::initializer_list<Point>& vertices_): vertices(vertices_) {}
    const std::vector<Point> getVertices() const {return vertices;}
    int verticesCount() const {return vertices.size();}
    bool isConvex() const;
    double perimeter() const override;
    double area() const override;
    void rotate(const Point& center, const double& angle) override;
    void reflex(const Point& center) override;
    void reflex(const Line& axis) override;
    void scale(const Point& center, const double& coefficient) override;
    bool operator==(const Shape& another) const override;
    bool operator!=(const Shape& another) const override { return !(*this == another); }
    bool isCongruentTo(const Shape& another) const override;
    bool isSimilarTo(const Shape& another) const override;
    bool containsPoint(Point point) const override;
    void print() const { for (auto u : vertices) std::cerr << u.x << " " << u.y << '\n'; }
    virtual ~Polygon() override = default;
};


bool Polygon::containsPoint(Point point) const {
    double angle = 0;
    Polygon cpy = *this;
    cpy.vertices.push_back(vertices[0]);
    for (size_t i = 0; i < cpy.vertices.size(); ++i) {
        if (Line(cpy.vertices[i], cpy.vertices[i + 1]).contains(point)) return true;
    }
    for (size_t i = 0; i < cpy.vertices.size() - 1; ++i) {
        Point vec1 = cpy.vertices[i] - point;
        Point vec2 = cpy.vertices[i + 1] - point;
        double tmp = std::acos((vec1 * vec2) / (vec1.size() * vec2.size()));
        if ((vec1 ^ vec2) >= 0) tmp = std::max(tmp, -tmp);
        else tmp = std::min(tmp, -tmp);
        angle += tmp;
    }
    if (abs(angle) < eps) return false;
    else return true;
}


bool Polygon::isSimilarTo(const Shape& another) const {
    std::cerr << "PolygonSim ";
    const Polygon* anoter_ = dynamic_cast<const Polygon*>(&another);
    if (anoter_ == nullptr) return false;
    const Polygon& other = dynamic_cast<const Polygon&>(another);
    if (verticesCount() != other.verticesCount()) return false;

    Polygon other_ = other;
    for (int k = 0; k < 2; ++k, other_.reflex(Line({0, 0}, {0, 1}))){
        for (size_t i = 0; i < vertices.size(); ++i) {
            Polygon cpy = other_;
            Point move_ = vertices[i] - cpy.vertices[0];
            for (auto &v : cpy.vertices) {
                v = v + move_;
            }

            for (size_t w = 0, j; w < 2; ++w) {
                if (w) j = (i == vertices.size() - 1) ? 0 : i + 1;
                else j = (i == 0) ? vertices.size() - 1 : i - 1;

                Point vec1 = vertices[j] - vertices[i];
                Point vec2 = cpy.vertices[1] - cpy.vertices[0];

                double angle = std::acos((vec1 * vec2) / (vec1.size() * vec2.size()));

                for (int k = 0; k < 2; ++k, angle *= -1) {
                    Polygon cpy_ = cpy;
                    cpy_.rotate(cpy_.vertices[0], angle);
                    cpy_.scale(cpy_.vertices[0], vec1.size() / vec2.size());

                    if (*this == cpy_) {
                        return true;
                    }
                }
            }
        }
    }
    return false;
}


bool Polygon::isCongruentTo(const Shape& another) const {
    const Polygon* anoter_ = dynamic_cast<const Polygon*>(&another);
    if (anoter_ == nullptr) return false;
    const Polygon& other = dynamic_cast<const Polygon&>(another);
    if (verticesCount() != other.verticesCount()) return false;

    if (isSimilarTo(other) && abs(perimeter() - other.perimeter()) < eps) {
        return true;
    } else {
        return true;
    }
}


bool Polygon::operator==(const Shape& another) const {
    const Polygon* anoter_ = dynamic_cast<const Polygon*>(&another);
    if (anoter_ == nullptr) return false;
    const Polygon& other = dynamic_cast<const Polygon&>(another);
    if (verticesCount() != other.verticesCount()) return false;

    for (auto u : vertices) {
        bool ret = false;
        for (auto v : other.vertices) {
            if (u == v) {
                ret = true;
                break;
            }
        }
        if (!ret) return false;
    }
    return true;
}


void Polygon::scale(const Point& center, const double& coefficient) {
    for (size_t i = 0; i < vertices.size(); ++i) {
        vertices[i].scale(center, coefficient);
    }
}


void Polygon::reflex(const Line& axis) {
    for (size_t i = 0; i < vertices.size(); ++i) {
        vertices[i] = axis.reflex(vertices[i]);
    }
}


void Polygon::reflex(const Point& center) {
    for (size_t i = 0; i < vertices.size(); ++i) {
        vertices[i].reflex(center);
    }
}


void Polygon::rotate(const Point& center, const double& angle) {
    for (size_t i = 0; i < vertices.size(); ++i) {
        vertices[i].rotate(center, angle);
    }
}


double Polygon::area() const {
    double ret = 0;
    for (size_t i = 0; i < vertices.size(); ++i) {
        size_t j = i + 1;
        if (i == vertices.size() - 1) j = 0;
        Point a = vertices[i], b = vertices[j];
        ret += a.x * b.y - a.y * b.x;
    }
    ret = std::max(ret, -ret);
    ret /= 2;
    return ret;
}


double Polygon::perimeter() const {
    double ret = 0;
    Point a, b;
    for (size_t i = 0; i < vertices.size() - 1; ++i) {
        a = vertices[i];
        b = vertices[i + 1];
        ret += std::sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
    }
    a = vertices[0];
    b = vertices.back();
    ret += std::sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
    return ret;
}


bool Polygon::isConvex() const {
    for (size_t i = 0, ii; i < vertices.size(); ++i) {
        if (i == vertices.size() - 1) {
            ii = 0;
        } else {
            ii = i + 1;
        }
        Line tmp(vertices[i], vertices[ii]);
        int bal = 0;
        for (size_t j = 0; j < vertices.size(); ++j) {
            double check = tmp.a * vertices[j].x + tmp.b * vertices[j].y + tmp.c;
            if (j == i || j == ii) continue;
            if (check > 0) {
                if (bal == -1) return false;
                if (bal == 0) {
                    bal = 1;
                    continue;
                }
            }
            if (check < 0) {
                if (bal == 1) return false;
                if (bal == 0) {
                    bal = -1;
                    continue;
                }
            }
        }
    }
    return true;
}



class Ellipse: public Shape {
protected:
    Point f1 = {0, 0}, f2 = {0, 0};
    double a = 1;
public:
    Ellipse() {}
    Ellipse(const Point& f1_, const Point& f2_, const double& dist): f1(f1_), f2(f2_), a(dist / 2) {}
    std::pair<Point,Point> focuses() const {return {f1, f2};}
    double eccentricity() const;
    std::pair<Line, Line> directrices() const;
    Point center() const {return {(f1.x + f2.x) / 2, (f1.y + f2.y) / 2};}
    virtual double perimeter() const override;
    virtual double area() const override;
    void rotate(const Point& center, const double& angle) override;
    void reflex(const Point& center) override;
    void reflex(const Line& axis) override;
    void scale(const Point& center, const double& coefficient) override;
    bool operator==(const Shape& another) const override;
    bool operator!=(const Shape& another) const override { return !(*this == another); };
    bool isCongruentTo(const Shape& another) const override;
    bool isSimilarTo(const Shape& another) const override;
    bool containsPoint(Point point) const;
    virtual ~Ellipse() override = default;
};


std::pair<Line, Line> Ellipse::directrices() const {
    double angle = std::atan((f2.y - f1.y) / (f2.x - f1.x));
    double e = eccentricity();
    Point A(a / e, 0);
    Point B(-a / e, 0);
    Point A_ = A + Point(0, 1);
    Point B_ = B + Point(0, 1);
    A.rotate(angle);
    A_.rotate(angle);
    B.rotate(angle);
    B_.rotate(angle);
    A = A + center();
    A_ = A_ + center();
    B = B + center();
    B_ = B_ + center();
    return {Line(A, A_), Line(B, B_)};
}


bool Ellipse::containsPoint(Point point) const {
    Point O = center();
    double c = a * eccentricity();
    double bb = a * a - c * c;
    double tmp = (point.x - O.x) * (point.x - O.x) / (a * a) + (point.y - O.y) * (point.y - O.y) / bb;
    if (tmp <= 1) return true;
    else return false;
}


bool Ellipse::isSimilarTo(const Shape& another) const {
    const Ellipse* anoter_ = dynamic_cast<const Ellipse*>(&another);
    if (anoter_ == nullptr) return false;
    const Ellipse& other = dynamic_cast<const Ellipse&>(another);

    if (abs(eccentricity() - other.eccentricity()) < eps) return true;
    else return false;
}


bool Ellipse::isCongruentTo(const Shape& another) const {
    const Ellipse* anoter_ = dynamic_cast<const Ellipse*>(&another);
    if (anoter_ == nullptr) return false;
    const Ellipse& other = dynamic_cast<const Ellipse&>(another);

    if (isSimilarTo(other) && perimeter() == other.perimeter()) return true;
    else return false;
}


bool Ellipse::operator==(const Shape& another) const {
    const Ellipse* anoter_ = dynamic_cast<const Ellipse*>(&another);
    if (anoter_ == nullptr) return false;
    const Ellipse& other = dynamic_cast<const Ellipse&>(another);

    if (((f1 == other.f1 && f2 == other.f2) || (f1 == other.f2 && f2 == other.f1)) && abs(a - other.a) < eps) return true;
    else return false;
}


void Ellipse::scale(const Point& center, const double& coefficient) {
    f1.scale(center, coefficient);
    f2.scale(center, coefficient);
    a *= coefficient;
}


void Ellipse::reflex(const Line& axis) {
    f1 = axis.reflex(f1);
    f2 = axis.reflex(f2);
}


void Ellipse::reflex(const Point& center) {
    f1.reflex(center);
    f2.reflex(center);
}


void Ellipse::rotate(const Point& center, const double& angle) {
    f1.rotate(center, angle);
    f2.rotate(center, angle);
}


double Ellipse::area() const {
    double c = f1.dist(f2) / 2;
    double b = std::sqrt(a * a - c * c);
    return PI * a * b;
}


double Ellipse::perimeter() const {
    double c = f1.dist(f2) / 2;
    double b = std::sqrt(a * a - c * c);
    return PI * (3 * (a + b) - std::sqrt((3 * a + b) * (a + 3 * b)));
}


double Ellipse::eccentricity() const {
    double c = f1.dist(f2) / 2;
    return c / a;
}


class Circle: public Ellipse {
public:
    Circle() {}
    Circle(const Point& o_, const double r_): Ellipse(o_, o_, r_ * 2) {}
    double radius() const {return a;}
    double perimeter() const override {return 2 * PI * a;}
    double area() const override {return PI * a * a;}
    ~Circle() override = default;
};


class Rectangle: public Polygon {
public:
    Rectangle() {}
    Rectangle(const Point& A_, const Point& C_, const double& tg_);
    std::pair<Line, Line> diagonals() const;
    Point center() const;
    double perimeter() const override {return 2 * (vertices[0].dist(vertices[1]) + vertices[0].dist(vertices[3]));}
    double area() const override {return vertices[0].dist(vertices[1]) * vertices[0].dist(vertices[3]);}
    virtual ~Rectangle() override = default;
};


Point Rectangle::center() const {
    auto diag = diagonals();
    return diag.first.intersec(diag.second);
}


std::pair<Line, Line> Rectangle::diagonals() const {
    Line AD(vertices[0], vertices[2]), BC(vertices[1], vertices[3]);
    return {AD, BC};
}


Rectangle::Rectangle(const Point& A_, const Point& C_, const double& tg_) {
    vertices.resize(4);
    vertices[0] = A_;
    vertices[2] = C_;
    double tg = tg_;
    if (tg < 1) tg = 1 / tg;
    double c = vertices[0].dist(vertices[2]);
    double a = c / (std::sqrt(1 + tg* tg));
    double b = a * tg;
    Point AC = vertices[2] - vertices[0];

    Point AB = AC * (a / c);
    AB.rotate(std::atan(tg));
    vertices[1] = vertices[0] + AB;

    Point AD = AC * (b / c);
    AD.rotate(std::atan(tg) + 3 * PI / 2);
    vertices[3] = vertices[0] + AD;
}


class Square: public Rectangle {
public:
    Square() {}
    Square(const Point& p1, const Point& p2): Rectangle(p1, p2, 1) {}
    Circle circumscribedCircle() const {return Circle(center(), vertices[0].dist(vertices[1]) / std::sqrt(2));}
    Circle circinscribedCircle() const {return Circle(center(), vertices[0].dist(vertices[1]) / 2);}
    ~Square() override = default;
};


class Triangle: public Polygon {
public:
    Triangle() {}
    Triangle(const Point& a, const Point& b, const Point& c);
    Circle circumscribedCircle() const;
    Circle inscribedCircle() const;
    Point centroid() const;
    Point orthocenter() const;
    Line EulerLine() const;
    Circle ninePointsCircle() const;
    ~Triangle() override = default;
};


Triangle::Triangle(const Point& a, const Point& b, const Point& c) {
    vertices.push_back(a);
    vertices.push_back(b);
    vertices.push_back(c);
}


Circle Triangle::ninePointsCircle() const {
    Point M = circumscribedCircle().center(), N = orthocenter();
    Point O = (M + N) / 2;
    Point ABm = (vertices[0] + vertices[1]) / 2;
    return Circle(O, O.dist(ABm));
}


Line Triangle::EulerLine() const {
    Point M = circumscribedCircle().center(), N = orthocenter();
    return Line(M, N);
}


Point Triangle::orthocenter() const {
    Line BC(vertices[1], vertices[2]);
    Line ha(vertices[0], Point(BC.a, BC.b), 1);
    Line CA(vertices[2], vertices[0]);
    Line hb(vertices[1], Point(CA.a, CA.b), 1);
    return ha.intersec(hb);
}


Point Triangle::centroid() const {
    Point M = vertices[0] + vertices[1] + vertices[2];
    return M / 3;
}


Circle Triangle::inscribedCircle() const {
    double a = vertices[1].dist(vertices[2]), b = vertices[2].dist(vertices[0]), c = vertices[0].dist(vertices[1]);
    Point M = (vertices[0] * a + vertices[1] * b + vertices[2] * c) / (a + b + c);
    return Circle(M, Line(vertices[0], vertices[1]).dist(M));
}


Circle Triangle::circumscribedCircle() const {
    double A = vertices[1].dist(vertices[2]), B = vertices[2].dist(vertices[0]), C = vertices[0].dist(vertices[1]);
    double S = area();
    double a = ((A * A) / (8 * S * S)) * ((vertices[0] - vertices[1]) * (vertices[0] - vertices[2]));
    double b = ((B * B) / (8 * S * S)) * ((vertices[1] - vertices[0]) * (vertices[1] - vertices[2]));
    double c = ((C * C) / (8 * S * S)) * ((vertices[2] - vertices[0]) * (vertices[2] - vertices[1]));
    Point M = vertices[0] * a + vertices[1] * b + vertices[2] * c;
    return Circle(M, M.dist(vertices[0]));
}
