#include <iostream>
#include <vector>
#include <string>
#include <complex>
#include <algorithm>

typedef unsigned char type;
typedef long long ll;
typedef std::complex<double> complex;
const long double PI = acosl(-1);
static const int DEF_BASE = 10;


class BigInteger {
private:
    static const int DIG_SZ = 100;

    int pos = 1;
    std::vector<type> val;

    void reduce() {
        while (val.size() > 1 && val.back() == 0) {
            val.pop_back();
        }
        if (val.empty())
            val.push_back(0);
    }

    static void fft(std::vector<complex>& poly, bool invert) {
        size_t sz = poly.size();

        size_t invi = 0;
        for (size_t i = 1; i < sz; ++i) {
            size_t bit = sz >> 1;

            while (invi >= bit) {
                invi -= bit;
                bit >>= 1;
            }
            invi += bit;

            if (i < invi)
                std::swap(poly[i], poly[invi]);
        }

        for (size_t len = 2; len <= sz; len <<= 1) {
            size_t mid = len >> 1;
            double angle = (2 * PI / len) * (invert ? -1 : 1);
            complex wf(cos(angle), sin(angle));

            for (size_t j = 0; j < sz; j += len) {
                complex ws = 1;

                for (size_t i = 0; i < mid; ++i, ws *= wf) {
                    complex first = poly[j + i];
                    complex second = poly[j + mid + i] * ws;

                    poly[j + i] = first + second;
                    poly[j + mid + i] = first - second;
                }
            }
        }

        if (invert) {
            for (size_t i = 0; i < sz; ++i) {
                poly[i] /= sz;
            }
        }
    }

    void elem_mul(type mul) {
        int ex = 0;
        int tmp = 0;
        for (type& u : val) {
            tmp = u;
            tmp *= mul;
            tmp += ex;

            u = tmp % DIG_SZ;
            ex = tmp / DIG_SZ;
        }
        if (ex)
            val.push_back(ex);
        reduce();
    }

    void elem_div(type div) {
        int ex = 0;
        for (int i = static_cast<int>(val.size()) - 1; i >= 0; --i) {
             int cur = val[i] + ex * DIG_SZ;
             val[i] = cur / div;
             ex = cur % div;
        }
        while (val.size() > 1 && val.back() == 0) {
            val.pop_back();
        }
        return;
    }

    type& operator[](const size_t ind) {
        while (val.size() < ind + 1) {
            val.push_back(0);
        }
        return val[ind];
    }

public:
    BigInteger();
    BigInteger(std::string s);
    BigInteger(int b);
    BigInteger(const BigInteger& bi) = default;
    BigInteger(const std::vector<unsigned char>& nval) : pos(true), val(nval) {}
    std::string toString() const;
    friend std::ostream& operator<<(std::ostream& out, const BigInteger& bi);
    friend std::istream& operator>>(std::istream& in, BigInteger& bi);
    ~BigInteger();

    bool operator==(const BigInteger& bi) const;
    bool operator!=(const BigInteger& bi) const;
    bool operator<(const BigInteger& bi) const;
    bool operator>(const BigInteger& bi) const;
    bool operator<=(const BigInteger& bi) const;
    bool operator>=(const BigInteger& bi) const;
    int abscomp(const BigInteger& other) const {
        if (val.size() == 1 && other.val.size() == 1 && val[0] == 0 && other.val[0] == 0)
            return 0;

        int comparator = static_cast<int>(val.size()) - other.val.size();
        if (comparator)
            return comparator;

        for (int i = static_cast<int>(val.size()) - 1; i >= 0; --i) {
            comparator = static_cast<int>(val[i]) - other.val[i];
            if (comparator)
                return comparator;
        }
        return 0;
    }

    BigInteger& operator=(const BigInteger&) = default;
    explicit operator bool() const;
    explicit operator int() const;
    BigInteger& operator+=(const BigInteger& bi);
    BigInteger& operator-=(const BigInteger& bi);
    BigInteger& operator++();
    BigInteger operator++(int);
    BigInteger& operator--();
    BigInteger operator--(int);
    BigInteger operator-() const {
        BigInteger ret = *this;
        if (*this != 0)
            ret.pos *= -1;
        return ret;
    }

    BigInteger& operator*=(const BigInteger& other);
    BigInteger& operator/=(const BigInteger& other);
    BigInteger& operator%=(const BigInteger& other);
    friend BigInteger abs(BigInteger);
};

BigInteger& BigInteger::operator%=(const BigInteger& other) {
    BigInteger temporary = *this;
    temporary /= other;
    temporary *= other;

    return *this -= temporary;
}

BigInteger& BigInteger::operator/=(const BigInteger& bi) {
    BigInteger self = *this;
    BigInteger other = bi;

    // It's Knuth's Algorithm
    // most names are similar to original ones
    // I tryed to make them as much logical as I can


    int newpos = self.pos * other.pos;
    other.pos = self.pos = 1;

    if (abscomp(other) < 0)
        return *this = 0;

    if (other.val.size() == 1) {
        self.pos = newpos;
        self.elem_div(other.val[0]);
        return *this = self;
    }

    size_t nsz = other.val.size();
    size_t msz = self.val.size() - nsz;

    BigInteger ex;
    ex.val.resize(nsz + 2, 0);
    ex.val[nsz + 1] = 1;

    BigInteger res;
    res.val.resize(msz + 1, 0);

    type div = DIG_SZ / (other[nsz - 1] + 1);
    self.elem_mul(div);
    other.elem_mul(div);

    for (int j = msz; j >= 0; --j) {
        int resdig = self[j + nsz] * DIG_SZ + self[j + nsz - 1];
        int resgidhelp = resdig % other[nsz - 1];
        resdig /= other[nsz - 1];

        while (resdig == DIG_SZ || resdig * other[nsz - 2] > DIG_SZ * resgidhelp + self[j + nsz - 2]) {
            --resdig;
            resgidhelp += other[nsz - 1];
            if (resgidhelp >= DIG_SZ)
                break;
        }

        BigInteger othercpy = other;
        othercpy.elem_mul(resdig);

        BigInteger selfcpy;
        selfcpy.val.clear();
        for (size_t i = j; i <= j + nsz; ++i) {
            selfcpy.val.push_back(self[i]);
        }

        selfcpy.reduce();
        selfcpy -= othercpy;

        bool flag = false;
        if (selfcpy < 0) {
            selfcpy += ex;
            flag = true;
        }

        res.val[j] = resdig;

        if (flag) {
            other[nsz] = 0;
            selfcpy += other;
            --res.val[j];

            if (selfcpy.val.size() == nsz + 2)
                selfcpy.val.pop_back();
        }

        for (size_t i = j; i <= j + nsz; ++i) {
            self[i] = selfcpy[i - j];
        }
    }

    self.val = res.val;
    self.reduce();
    self.pos = newpos;
    return *this = self;
}

BigInteger& BigInteger::operator*=(const BigInteger& bi) {
    pos *= bi.pos;

    if (val.size() == 1 || bi.val.size() == 1) {
        bool fl = val.size() == 1;

        type tomul = (fl ? val[0] : bi.val[0]);
        int newpos = pos;
        if (fl)
            *this = bi;

        elem_mul(tomul);
        pos = newpos;
        return *this;
    }

    std::vector<complex> first;
    std::vector<complex> second;


    for (type x : val) {
        first.push_back(x);
    }
    for (type x : bi.val) {
        second.push_back(x);
    }

    size_t sz = 1;
    while (sz <= std::max(first.size(), second.size())) {
        sz *= 2;
    }
    sz *= 2;
    first.resize(sz);
    second.resize(sz);

    fft(first, false);
    fft(second, false);
    for (size_t i = 0; i < sz; ++i) {
        first[i] *= second[i];
    }
    fft(first, true);

    val.clear();
    int mem = 0;
    for (size_t i = 0; mem != 0 || i < sz; ++i) {
        int tmp = mem;
        if (i < sz)
            tmp += static_cast<int>(first[i].real() + 0.5);
        mem = tmp / DIG_SZ;
        val.push_back(tmp - mem * DIG_SZ);
    }

    reduce();
    return *this;
}

BigInteger& BigInteger::operator-=(const BigInteger& bi) {
    BigInteger bii = bi;
    if (pos != bii.pos) {
        return *this += -bi;
    }
    if ((*this < bii) ^ (pos == -1)) {
        std::swap(*this, bii);
        pos *= -1;
    }
    val.push_back(0);
    int ex = 0;
    int fi = 0;
    int se = 0;
    for (size_t i = 0; i < val.size(); ++i) {
        fi = val[i];
        fi -= ex;
        ex = 0;
        if (i < bii.val.size()) se = bii.val[i];
        if (fi < se) {
            fi += DIG_SZ;
            ex = 1;
        }
        fi -= se;
        val[i] = fi;

        fi = se = 0;
    }

    reduce();
    return *this;
}


BigInteger& BigInteger::operator+=(const BigInteger& bi) {
    if (bi == 0) return *this;
    if (*this == 0) return *this = bi;
    BigInteger bii = bi;
    if (pos != bii.pos) {
        bii = -bi;
        return *this -= bii;
    }
    if ((*this < bii) ^ (pos == -1)) std::swap(*this, bii);

    bii.val.resize(val.size());
    int ex = 0;
    for (size_t i = 0; i < val.size(); ++i) {
        val[i] += bii.val[i] + ex;
        ex = val[i] / DIG_SZ;
        val[i] %= DIG_SZ;
        if (i == val.size() - 1 && ex > 0) {
            val.push_back(0);
            bii.val.push_back(0);
        }
    }
    reduce();
    return *this;
}

BigInteger& BigInteger::operator++() {
    return *this += BigInteger(1);
}

BigInteger& BigInteger::operator--() {
    return *this += BigInteger(-1);
}

BigInteger BigInteger::operator++(int) {
    const BigInteger ret = *this;
    ++*this;
    return ret;
}

BigInteger BigInteger::operator--(int) {
    BigInteger ret = *this;
    --*this;
    return ret;
}

BigInteger::operator int() const {
    int ret = 0;
    for (int i = val.size() - 1; i >= 0; --i) {
        ret *= DIG_SZ;
        ret += val[i];
    }
    if (!pos) ret *= -1;
    return ret;
}

BigInteger::operator bool() const {
    return val.size() > 1 || val[0] != 0;
}

BigInteger::~BigInteger() {
    val.clear();
}

std::ostream& operator<<(std::ostream& out, const BigInteger& bi) {
    out << bi.toString();
    return out;
}

std::istream& operator>>(std::istream& in, BigInteger& bi) {
    std::string str;
    in >> str;
    BigInteger bii(str);
    bi = bii;
    if (bi == 0) bi.pos = 1;
    return in;
}

std::string BigInteger::toString() const {
    std::string ret;
    if (pos == -1 && *this != 0) ret.push_back('-');
    auto dig = val.back();
    if (dig / DEF_BASE) ret.push_back(dig / DEF_BASE + '0');
    ret.push_back(dig % DEF_BASE + '0');
    for (int i = val.size() - 2; i >= 0; --i) {
        ret.push_back(val[i] / DEF_BASE + '0');
        ret.push_back(val[i] % DEF_BASE + '0');
    }
    return ret;
}

BigInteger::BigInteger(int nval): pos(nval < 0 ? -1 : 1) {
    nval = abs(nval);
    while (nval > 0) {
        unsigned long long oval = nval;
        nval /= DIG_SZ;
        val.push_back(oval - DIG_SZ * nval);
    }
    if (val.empty())
        val.push_back(0);
}

BigInteger::BigInteger(std::string str) {
    std::reverse(str.begin(), str.end());
    if (str.back() == '-') {
        pos = -1;
        str.pop_back();
    }
    while (str.back() == '0' && str.size() > 1) {
        str.pop_back();
    }
    for (size_t i = 0; i < str.size(); i += 2) {
        unsigned char dig = str[i] - '0';
        if (i != str.size() - 1) {
            dig += (str[i + 1] - '0') * DEF_BASE;
        }
        val.push_back(dig);
    }
}

BigInteger::BigInteger(): BigInteger(0) {}

BigInteger operator+(const BigInteger& bi, const BigInteger& bii) {
    BigInteger ret = bi;
    ret += bii;
    return ret;
}

BigInteger operator-(const BigInteger& bi, const BigInteger& bii) {
    BigInteger ret = bi;
    ret -= bii;
    return ret;
}

BigInteger operator*(const BigInteger& bi, const BigInteger& bii) {
    BigInteger ret = bi;
    ret *= bii;
    return ret;
}

BigInteger operator/(const BigInteger& bi, const BigInteger& bii) {
    BigInteger ret = bi;
    ret /= bii;
    return ret;
}

BigInteger operator%(const BigInteger& bi, const BigInteger& bii) {
    BigInteger ret = bi;
    ret %= bii;
    return ret;
}

bool BigInteger::operator>=(const BigInteger& bi) const {
    return !(*this < bi);
}

bool BigInteger::operator<=(const BigInteger& bi) const {
    return !(*this > bi);
}

bool BigInteger::operator>(const BigInteger& bi) const {
    return bi < *this;
}

bool BigInteger::operator<(const BigInteger& bi) const {
    if (pos < bi.pos) return true;
    if (pos > bi.pos) return false;
    bool fl = true;
    if (pos == -1 && bi.pos == -1) fl = false;
    if (val.size() < bi.val.size()) return fl;
    if (val.size() > bi.val.size()) return !fl;
    for (int i = val.size() - 1; i >= 0; --i) {
        if (val[i] < bi.val[i]) return fl;
        if (val[i] > bi.val[i]) return !fl;
    }
    return false;
}

bool BigInteger::operator!=(const BigInteger& bi) const {
    return !(*this == bi);
}

bool BigInteger::operator==(const BigInteger& bi) const {
    if (pos != bi.pos) return false;
    if (val.size() != bi.val.size()) return false;
    for (size_t i = 0; i < val.size(); ++i) {
        if (val[i] != bi.val[i]) return false;
    }
    return true;
}

BigInteger abs(BigInteger bi) {
    bi.pos = 1;
    return bi;
}

BigInteger GCD(BigInteger first, BigInteger second) {
    first = abs(first);
    second = abs(second);

    while (first > 0 && second > 0) {
        if (first < second) {
            second %= first;
        } else {
            first %= second;
        }
    }

    return first > 0 ? first : second;
}

BigInteger LCM(BigInteger first, BigInteger second) {
    return first / GCD(first, second) * second;
}


class Rational {
private:
    BigInteger nom = {1};
    BigInteger den = {1};
public:
    Rational();
    Rational(const BigInteger& bi);
    Rational(int x);
    ~Rational() {}
    Rational& operator=(const Rational&) = default;
    Rational& operator+=(const Rational& r);
    Rational& operator-=(const Rational& r);
    Rational& operator*=(const Rational& r);
    Rational& operator/=(const Rational& r);
    bool operator<(const Rational& r) const;
    bool operator==(const Rational& r) const;
    bool operator>(const Rational& r) const;
    bool operator<=(const Rational& r) const;
    bool operator>=(const Rational& r) const;
    bool operator!=(const Rational& r) const;
    void reduce();
    std::string toString() const;
    std::string asDecimal(size_t precision = 0) const;
    explicit operator double() const;
};

std::ostream& operator<<(std::ostream& out, const Rational& rat) {
    out << rat.toString();
    return out;
}


std::istream& operator>>(std::istream& in, Rational& rat) {
    BigInteger t;
    in >> t;
    rat = Rational(t);
    return in;
}

Rational operator-(const Rational& r) {
    if (r != 0) {
        Rational rt(r);
        rt *= -1;
        return rt;
    } else {
        return r;
    }
}

Rational::operator double() const {
    return std::stod(asDecimal(16));
}

std::string Rational::asDecimal(size_t precision) const {
    std::string ret;
    Rational rt = *this;
    if (rt.nom < 0) rt.nom *= -1;
    Rational rtt = rt;
    BigInteger bi = 1;
    for (size_t i = 0; i < precision; ++i) {
        bi *= DEF_BASE;
    }
    BigInteger exp = rt.nom;
    exp *= bi;
    exp /= rt.den;
    ret = exp.toString();
    if (precision) {
        if (ret.size() > precision) {
            ret = ret.substr(0, ret.size() - precision) + "." + ret.substr(ret.size() - precision);
        } else {
            std::string nul(precision - ret.size(), '0');
            ret = "0." + nul + ret;
        }
    }
    if (nom < 0 && ret != "0") ret = "-" + ret;
    return ret;
}

std::string Rational::toString() const {
    std::string ret;
    Rational rt = *this;
    rt.reduce();
    if (rt.nom % rt.den == 0) ret = (rt.nom / rt.den).toString();
    else ret = rt.nom.toString() + '/' + rt.den.toString();
    return ret;
}

void Rational::reduce() {
    Rational rt = *this;
    if (rt.nom < 0) rt.nom *= BigInteger(-1);
    BigInteger bi = GCD(rt.nom, rt.den);
    nom /= bi;
    den /= bi;
}

bool Rational::operator!=(const Rational& r) const {
    return nom * r.den != den * r.nom;
}

bool Rational::operator>=(const Rational& r) const {
    return nom * r.den >= den * r.nom;
}

bool Rational::operator<=(const Rational& r) const {
    return nom * r.den <= den * r.nom;
}

bool Rational::operator>(const Rational& r) const {
    return nom * r.den > den * r.nom;
}

bool Rational::operator==(const Rational& r) const {
    return nom * r.den == den * r.nom;
}

bool Rational::operator<(const Rational& r) const {
    return nom * r.den < den * r.nom;
}

Rational operator+(const Rational& r, const Rational& ri) {
    Rational cpy = r;
    return cpy += ri;
}

Rational operator-(const Rational& r, const Rational& ri) {
    Rational cpy = r;
    return cpy -= ri;
}

Rational operator*(const Rational& r, const Rational& ri) {
    Rational cpy = r;
    return cpy *= ri;
}

Rational operator/(const Rational& r, const Rational& ri) {
    Rational cpy = r;
    return cpy /= ri;
}

Rational& Rational::operator*=(const Rational& r) {
    nom *= r.nom;
    den *= r.den;
    reduce();
    return *this;
}

Rational& Rational::operator/=(const Rational& rt) {
    Rational rti;
    rti.nom = rt.den;
    rti.den = rt.nom;
    if (rti.den < 0) {
        rti.nom = -rti.nom;
        rti.den = -rti.den;
    }
    rti.reduce();
    return *this *= rti;
}

Rational& Rational::operator+=(const Rational& rt) {
    Rational tmp = rt;
    BigInteger nden = LCM(den, tmp.den);

    nom *= nden / den;
    tmp.nom *= nden / tmp.den;
    den = tmp.den = nden;
    nom += tmp.nom;

    reduce();
    return *this;
}

Rational& Rational::operator-=(const Rational& rt) {
    return *this += -rt;
}

Rational::Rational(int x) : nom(x) {}

Rational::Rational(const BigInteger& bi) : nom(bi) {}

Rational::Rational() : nom(0) {}
