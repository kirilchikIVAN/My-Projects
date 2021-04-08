#include <iostream>
#include <vector>
#include <string>
#include <complex>
#include <algorithm>

typedef unsigned char type;
typedef long long ll;
typedef std::complex<double> complex;
const long double PI = acosl(-1);


template <typename T>
T myabs(const T& a) {
    return a * T(-1) < a ? a : a * T(-1);
}


class BigInteger {
private:
    static const int DIG_SZ = 100;
    static const int DEF_BASE = 10;

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

public:
    BigInteger();
    BigInteger(std::string s);
    BigInteger(int b);
    BigInteger(const BigInteger& bi) = default;
    BigInteger(const std::vector<unsigned char>& val_) : pos(true), val(val_) {}
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

    type& operator[](const size_t ind) {
        while (val.size() < ind + 1) {
            val.push_back(0);
        }
        return val[ind];
    }
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
            fi += 100;
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
    return (val.size() > 1 || val[0] != 0) ? true : false;
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

BigInteger::BigInteger(int val_): pos(val_ < 0 ? -1 : 1) {
        val_ = abs(val_);
        while (val_ > 0) {
            unsigned long long oval = val_;
            val_ /= DIG_SZ;
            val.push_back(oval - DIG_SZ * val_);
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
        unsigned char dig = str[i] - 48;
        if (i != str.size() - 1) {
            dig += (str[i + 1] - 48) * 10;
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
    return (*this > bi || *this == bi) ? 1 : 0;
}

bool BigInteger::operator<=(const BigInteger& bi) const {
    return (*this < bi || *this == bi) ? 1 : 0;
}

bool BigInteger::operator>(const BigInteger& bi) const {
    return (*this < bi || *this == bi) ? 0 : 1;
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
    BigInteger nom, den = {1};
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

const Rational operator-(const Rational& r) {
    if (r != 0) {
        Rational rt(r);
        rt *= -1;
        return rt;
    } else return r;
}

Rational::operator double() const {
    std::string str = asDecimal(16);
    double ret = 0;
    int fl = 1, d = 10;
    for (size_t i = str[0] == '-' ? 1 : 0; i < str.size(); ++i) {
        if (str[i] == '.') {
            fl = 0;
            continue;
        }
        if (fl) {
            ret *= 10;
            ret += str[i] - 48;
        } else {
            ret += (double)(str[i] - 48) / d;
            d *= 10;
        }
    }
    if (str[0] == '-')
        ret *= -1;
    return ret;
}

std::string Rational::asDecimal(size_t precision) const {
    std::string ret;
    Rational rt = *this;
    if (rt.nom < 0) rt.nom *= -1;
    Rational rtt = rt;
    BigInteger bi = 1;
    for (size_t i = 0; i < precision; ++i) {
        bi *= 10;
    }
    BigInteger exp = rt.nom;
    exp *= bi;
    exp /= rt.den;
    ret = exp.toString();
    if (precision) {
        if (ret.size() > precision)
            ret = ret.substr(0, ret.size() - precision) + "." + ret.substr(ret.size() - precision);
        else {
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

Rational::Rational(int x) {
    nom = x;
}

Rational::Rational(const BigInteger& bi) {
    nom = bi;
}

Rational::Rational() {
    nom = {0};
}






















const unsigned MAXINT = 2147483647;


template <unsigned long long N, unsigned long long M, bool P>
struct is_prime_helper {
    static const bool flag = (M * M > N);
    static const unsigned long long newM = flag ? N + 2 : M + 2;
    static const unsigned minden = is_prime_helper<N, newM, flag || (N % M == 0)>::minden;
};

template <unsigned long long N, unsigned long long M>
struct is_prime_helper<N, M, true> {
    static const unsigned minden = M - 2;
};


template <unsigned N>
struct is_prime {
    static const unsigned minden = N % 2 == 0 ? 2 : is_prime_helper<N, 3, false>::minden;
    static const bool value = minden == N;
};

template <>
struct is_prime<2> {
    static const bool value = 1;
    static const unsigned minden = 2;
};

template <>
struct is_prime<1> {
    static const bool value = 0;
    static const unsigned minden = 1;
};

template <>
struct is_prime<0> {
    static const bool value = 0;
    static const unsigned minden = 1;
};


template <unsigned N>
const bool is_prime_v = is_prime<N>::value;

template <unsigned N>
const unsigned minden_v = is_prime<N>::minden;

/*---------------------------------*/

template <unsigned N, unsigned M, bool P>
struct is_deg {
    static const bool flag = N % M != 0 || M == 1 || N == 0;
    static const bool value = is_deg<N / M, M, flag>::value;
};

template <unsigned N, unsigned M>
struct is_deg<N, M, 1> {
    static const bool value = false;
};


template <unsigned N>
struct is_deg<N, N, 0> {
    static const bool value = true;
};


template <unsigned N, unsigned M>
const bool is_deg_v = is_deg<N, M, 0>::value;

template <unsigned N>
const bool is_deg_v<N, 0> = false;

/*-------------------------------*/

template <unsigned N>
struct has_primitive_root {
    static const bool flag = (N % 2 == 0);
    static const unsigned newN = N / (flag ? 2 : 1);
    static const unsigned newN_minden = minden_v<newN>;
    static const bool value = (N % 4 == 0) ?
                            false :
                            is_deg_v<newN, newN_minden>;
};

template <>
struct has_primitive_root<0> {
    static const bool value = false;
};

template <>
struct has_primitive_root<1> {
    static const bool value = false;
};

template <>
struct has_primitive_root<2> {
    static const bool value = true;
};

template <>
struct has_primitive_root<4> {
    static const bool value = true;
};


template <unsigned N>
const bool has_primitive_root_v = has_primitive_root<N>::value;


/*--------------------*/


unsigned EulerF(unsigned n) {
    unsigned euler = n;
    for (unsigned i = 2; i * i <= n; ++i) {
        if (n % i == 0) {
            while (n % i == 0) {
                n /= i;
            }
            euler -= euler / i;
        }
    }
    if (n > 1) euler -= euler / n;
    return euler;
}


unsigned GCD(unsigned a, unsigned b) {
    a = myabs(a);
    b = myabs(b);

    while (0 < a && 0 < b) {
        a < b ? b %= a : a %= b;
    }

    return 0 < a ? a : b;
}


template <bool fl>
struct err;

template<>
struct err<true> {};


template <unsigned N>
class Residue {
private:
    unsigned val = 0;
public:
    Residue() {}
    Residue(const unsigned& val_): val(val_ % N) {}
    Residue(const Residue<N>& res): val(res.val) {}
    Residue(const int& val_) {
        long long cpy = val_;
        unsigned big = MAXINT;
        int n = big > N ? big / N : 1;
        if (cpy < 0) cpy += n * N;
        val = cpy % N;
    }
    explicit operator int() const {
        return (int)val;
    }
    Residue<N>& operator=(const Residue<N>& res) {
        val = res.val;
        return *this;
    }
    Residue<N>& operator+=(const Residue<N>& res) {
        unsigned long long tmp = val;
        tmp += res.val;
        tmp %= N;
        val = tmp;
        return *this;
    }
    Residue<N>& operator-=(const Residue<N>& res) {
        unsigned long long tmp = val;
        if (tmp < res.val) {
            tmp += N;
        }
        tmp -= res.val;
        val = tmp;
        return *this;
    }
    Residue<N>& operator*=(const Residue<N>& res) {
        unsigned long long tmp = val;
        tmp *= res.val;
        tmp %= N;
        val = tmp;
        return *this;
    }
    Residue<N> pow(const unsigned& k) const;
    Residue<N> getInverse() const {
        err<is_prime_v<N>>();
        return pow(N - 2);
    }
    Residue<N>& operator/=(const Residue<N>& res) {
        err<is_prime_v<N>>();
        Residue<N> invres = res.getInverse();
        return *this *= invres;
    }
    unsigned order() const;
    static Residue<N> getPrimitiveRoot();
    bool operator<(const Residue<N>& res) const {
        return val < res.val;
    }

    bool operator==(const Residue<N>& res) const {
        return val == res.val;
    }
    bool operator!=(const Residue<N>& res) const {
        return !(*this == res);
    }

    template <unsigned M>
    friend std::ostream& operator<<(std::ostream& out, const Residue<M>& res);
    template <unsigned M>
    friend std::istream& operator>>(std::istream& in, Residue<M>& res);
};


template <unsigned N>
std::ostream& operator<<(std::ostream& out, const Residue<N>& res) {
    out << res.val;
    return out;
}

template <unsigned N>
std::istream& operator>>(std::istream& in, Residue<N>& res) {
    int tmp;
    in >> tmp;
    res = Residue<N>(tmp);
    return in;
}


template <unsigned N>
Residue<N> Residue<N>::getPrimitiveRoot() {
    err<has_primitive_root_v<N>>();

    unsigned euler = EulerF(N);

    std::vector<unsigned> fact;
    unsigned n = euler;
    for (unsigned i = 2; i * i <= n; ++i) {
        if (n % i == 0) {
            fact.push_back(i);
            while (n % i == 0) {
                n /= i;
            }
        }
    }
    if (n > 1) fact.push_back(n);

    for (unsigned tmp = 2; tmp <= N; ++tmp) {
        bool ok = true;
        if (GCD(tmp, N) != 1) {
            ok = false;
        }
        for (size_t i = 0; i < fact.size() && ok; ++i) {
            ok &= Residue<N>(tmp).pow(euler / fact[i]).val != 1;
        }
        if (ok)  return Residue<N>(tmp);
    }
    return 0;
}


template <unsigned N>
unsigned Residue<N>::order() const {
    if (GCD(val, N) != 1) return 0;

    unsigned euler = EulerF(N);

    unsigned karm = euler / (has_primitive_root_v<N> ? 1 : 2);
    unsigned ret = karm;
    for (unsigned i = 1; i * i <= karm; ++i) {
        if (karm % i == 0) {
            if (pow(i).val == 1) ret = std::min(ret, i);
            if (pow(karm / i).val == 1) ret = std::min(ret, karm / i);
        }
    }
    return ret;
}


template <unsigned N>
Residue<N> Residue<N>::pow(const unsigned& k) const {
    if (k == 0) return Residue<N>(1);
    if (k == 1) return *this;
    unsigned kk = k;
    std::vector<int> tmp;
    while (kk != 0) {
        tmp.push_back(1 & kk);
        kk >>= 1;
    }
    Residue<N> ret = *this;
    for (size_t i = 2; i <= tmp.size(); ++i) {
        ret *= ret;
        if (tmp[tmp.size() - i] == 1) ret *= *this;
    }
    return ret;
}


template <unsigned N>
Residue<N> operator+(const Residue<N>& res1, const Residue<N>& res2) {
    Residue<N> ret(res1);
    ret += res2;
    return ret;
}


template <unsigned N>
Residue<N> operator-(const Residue<N>& res1, const Residue<N>& res2) {
    Residue<N> ret(res1);
    ret -= res2;
    return ret;
}


template <unsigned N>
Residue<N> operator*(const Residue<N>& res1, const Residue<N>& res2) {
    Residue<N> ret(res1);
    ret *= res2;
    return ret;
}


template <unsigned N>
Residue<N> operator/(const Residue<N>& res1, const Residue<N>& res2) {
    Residue<N> ret(res1);
    ret /= res2;
    return ret;
}


template<unsigned N>
Residue<N> myabs(const Residue<N>& res) {
    return res;
}

























template <typename T>
void swap(T& a, T& b) {
    T t = a;
    a = b;
    b = t;
}


template <unsigned N, unsigned M, typename Field = Rational>
class Matrix {
public:
    std::vector<std::vector<Field>> val;
public:
    Matrix(): val({N, std::vector<Field>(M)}) {}
    Matrix(const std::initializer_list<std::vector<Field>>& val_): val(val_) {}
    Matrix(const Matrix<N, M, Field>& mat): val(mat.val) {}
    Matrix(std::vector<std::vector<Field>> val_): val(val_) {}
    Matrix<N, M, Field>& operator=(const Matrix<N, M, Field>& mat) = default;
    ~Matrix() {}

    Matrix<N, M, Field>& operator+=(const Matrix<N, M, Field>& mat);
    Matrix<N, M, Field>& operator-=(const Matrix<N, M, Field>& mat);
    Matrix<N, M, Field>& operator*=(const Field& fi);

    Field det() const;
    Matrix<M, N, Field> transposed() const;
    unsigned rank() const;
    Matrix<N, M, Field> inverted() const;
    Matrix<N, M, Field>& invert();

    std::vector<Field> getRow(const unsigned& row) const;
    std::vector<Field> getColumn(const unsigned& col) const;
    const std::vector<Field>& operator[](const unsigned& row) const { return val[row]; }
    std::vector<Field>& operator[](const unsigned& row) { return val[row]; }

    unsigned col_max(const unsigned& col) const;
    unsigned triangulation();

    Field trace() const;

    bool operator==(const Matrix<N, M, Field>& mat) const;
    bool operator!=(const Matrix<N, M, Field>& mat) const { return !(*this == mat); }

    Matrix<N, M, Field>& operator*=(const Matrix<N, M, Field>& mat);
};


template <unsigned N, unsigned M, typename Field>
bool Matrix<N, M, Field>::operator==(const Matrix<N, M, Field>& mat) const {
    for (unsigned i = 0; i < N; ++i) {
        for (unsigned j = 0; j < M; ++j) {
            if (val[i][j] != mat.val[i][j]) return false;
        }
    }
    return true;
}


template <unsigned N, unsigned M, typename Field>
Field Matrix<N, M, Field>::trace() const {
    err<N == M>();

    Field ret(0);
    for (unsigned i = 0; i < N; ++i) {
        ret += val[i][i];
    }
    return ret;
}


template <unsigned N, unsigned M, typename Field>
std::vector<Field> Matrix<N, M, Field>::getColumn(const unsigned& col) const {
    std::vector<Field> ret(N);
    for (unsigned i = 0; i < N; ++i) {
        ret.push_back(val[i][col]);
    }
    return ret;
}


template <unsigned N, unsigned M, typename Field>
std::vector<Field> Matrix<N, M, Field>::getRow(const unsigned& row) const {
    std::vector<Field> ret(val[row]);
    return ret;
}


template <unsigned N, unsigned M, typename Field>
Matrix<M, N, Field> Matrix<N, M, Field>::transposed() const {
    Matrix<M, N, Field> res;
    for (unsigned i = 0; i < N; ++i) {
        for (unsigned j = 0; j < M; ++j) {
            res.val[j][i] = val[i][j];
        }
    }
    return res;
}


template <unsigned N, unsigned M, typename Field>
Matrix<N, M, Field>& Matrix<N, M, Field>::operator*=(const Field& fi) {
    for (unsigned i = 0; i < N; ++i){
        for (unsigned j = 0; j < M; ++j) {
            val[i][j] *= fi;
        }
    }
    return *this;
}


template <unsigned N, unsigned M, typename Field>
Matrix<N, M, Field>& Matrix<N, M, Field>::operator-=(const Matrix<N, M, Field>& mat) {
    for (unsigned i = 0; i < N; ++i){
        for (unsigned j = 0; j < M; ++j) {
            val[i][j] -= mat.val[i][j];
        }
    }
    return *this;
}


template <unsigned N, unsigned M, typename Field>
Matrix<N, M, Field>& Matrix<N, M, Field>::operator+=(const Matrix<N, M, Field>& mat) {
    for (unsigned i = 0; i < N; ++i){
        for (unsigned j = 0; j < M; ++j) {
            val[i][j] += mat.val[i][j];
        }
    }
    return *this;
}


template <unsigned N, typename Field = Rational>
using SquareMatrix = Matrix<N, N, Field>;


template <unsigned N, unsigned M, typename Field>
Matrix<N, M, Field> operator+(const Matrix<N, M, Field>& m1, const Matrix<N, M, Field>& m2) {
    Matrix<N, M, Field> cpy = m1;
    cpy += m2;
    return cpy;
}


template <unsigned N, unsigned M, typename Field>
Matrix<N, M, Field> operator-(const Matrix<N, M, Field>& m1, const Matrix<N, M, Field>& m2) {
    Matrix<N, M, Field> cpy = m1;
    cpy -= m2;
    return cpy;
}


template <unsigned N, unsigned M, typename Field>
Matrix<N, M, Field> operator*(const Matrix<N, M, Field>& mat, const Field& fi) {
    Matrix<N, M, Field> cpy = mat;
    cpy *= fi;
    return cpy;
}


template <unsigned N, unsigned M, typename Field>
Matrix<N, M, Field> operator*(const Field& fi, const Matrix<N, M, Field>& mat) {
    Matrix<N, M, Field> cpy = mat;
    cpy *= fi;
    return cpy;
}


template <unsigned N, unsigned M, unsigned K, typename Field>
Matrix<N, K, Field> operator*(const Matrix<N, M, Field>& m1, const Matrix<M, K, Field>& m2) {
    Matrix<N, K, Field> res;
    for (unsigned i = 0; i < N; ++i) {
        for (unsigned j = 0; j < K; ++j) {
            res[i][j] = 0;
            for (unsigned t = 0; t < M; ++t) {
                res[i][j] += m1.val[i][t] * m2.val[t][j];
            }
        }
    }
    return res;
}


template <unsigned N, unsigned M, typename Field>
unsigned Matrix<N, M, Field>::triangulation() {
    unsigned swapCount = 0;

    const unsigned n_col = M;
    for (unsigned i = 0; i < N - 1; ++i) {
        unsigned imax = col_max(i);
        if (i != imax) {
            swap(val[i], val[imax]);
            ++swapCount;
        }
        if (val[i][i] != Field(0)) {
            for (unsigned j = i + 1; j < N; ++j) {
                Field coef = val[j][i] / val[i][i];
                for (unsigned k = i; k < n_col; ++k) {
                    if (i == 14) {
                    }
                    val[j][k] -= val[i][k] * coef;
                }
            }
        }
    }
    return swapCount;
}


template <unsigned N, unsigned M, typename Field>
unsigned Matrix<N, M, Field>::col_max(const unsigned& col) const {
    for (unsigned i = col; i < N; ++i) {
        if (val[i][col] != 0) {
            return i;
        }
    }
    return col;
}


template <unsigned N, unsigned M, typename Field>
Field Matrix<N, M, Field>::det() const {
    err<N == M>();
    Matrix<N, M, Field> cpy = *this;
    unsigned swapCount = cpy.triangulation();

    long long n = std::min(N, M);
    Field ret = 1;
    for (unsigned i = 0; i < n; ++i) {
        if (cpy.val[i][i] == Field(0)) return 0;
        ret *= cpy.val[i][i];
    }

    if (swapCount & 1) ret *= Field(-1);
    return ret;
}


template <unsigned N, unsigned M, typename Field>
unsigned Matrix<N, M, Field>::rank() const {
    Matrix<N, M, Field> cpy = *this;
    cpy.triangulation();

    unsigned ret = 0;
    unsigned i = 0, j = 0;
    for (; i < N && j < M; ) {
        if (cpy.val[i][j] != Field(0)) {
            ++ret;
            ++i;
        }
        ++j;
    }
    return ret;
}


template <unsigned N, unsigned M, typename Field>
Matrix<N, M, Field>& Matrix<N, M, Field>::invert() {
    err<M == N>();
    Matrix<M, N, Field> tmp = inverted();
    for (unsigned i = 0; i < N; ++i) {
        for (unsigned j = 0; j < M; ++j) {
            val[i][j] = tmp.val[i][j];
        }
    }
    return *this;
}


template <unsigned N, unsigned M, typename Field>
Matrix<N, M, Field> Matrix<N, M, Field>::inverted() const {
    Matrix<N, M * 2, Field> tmp;

    for (unsigned i = 0; i < N; ++i) {
        for (unsigned j = 0; j < M; ++j) {
            tmp.val[i][j] = val[i][j];
        }
    }
    for (unsigned i = 0; i < N; ++i) {
        for (unsigned j = 0; j < M; ++j) {
            tmp.val[i][j + M] = i == j ? Field(1) : Field(0);
        }
    }

    tmp.triangulation();

    for (unsigned i = N - 1; ; --i) {
        if (i > 0) {
            for (unsigned ii = i - 1; ; --ii) {
                Field coef = tmp.val[ii][i] / tmp.val[i][i];
                for (unsigned j = i; j < 2 * M; ++j) {
                    tmp.val[ii][j] -= tmp.val[i][j] * coef;
                }
                if (ii == 0) break;
            }
        }
        for (unsigned j = i + 1; j < 2 * M; ++j) {
            tmp.val[i][j] /= tmp.val[i][i];
        }
        if (i == 0) break;
    }

    Matrix<M, N, Field> res;
    for (unsigned i = 0; i < N; i++) {
        for (unsigned j = 0; j < M; j++) {
            res.val[i][j] = tmp.val[i][j + M];
        }
    }

    return res;
}


template <unsigned N, unsigned M, typename Field>
Matrix<N, M, Field>& Matrix<N, M, Field>::operator*=(const Matrix<N, M, Field>& mat) {
    err<N == M>();
    *this = *this * mat;
    return *this;
}
