#include <iostream>
#include <vector>


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
    static const bool value = (N % 4 == 0) ? false : is_deg_v<newN, newN_minden>;
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
    if (a < b) std::swap(a, b);
    if (b == 0) return a;
    return GCD(b, a % b);
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
    explicit Residue(const int& val_) {
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
};


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
