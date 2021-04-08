#include <iostream>
#include <cstring>

class String {
private:
    size_t len;
    size_t cap;
    char* str;
    void swap(String tmp) {
        std::swap(len, tmp.len);
        std::swap(cap, tmp.cap);
        std::swap(str, tmp.str);
    }
public:
    ~String();
    String();
    String(char c);
    String(const char* ptr);
    String(size_t n, char c);
    String(const String& tmp);
    String(const String& tmp, size_t s, size_t n);
    String& operator=(const String& tmp);
    const char& operator[](size_t n) const;
    char& operator[](size_t n);
    bool operator==(const String& tmp) const;
    size_t length() const;
    void resize(size_t n);
    void push_back(char c);
    void pop_back();
    char& front();
    char& back();
    const char& front() const;
    const char& back() const;
    String& operator+=(const String& tmp);
    void reverse();
    size_t find(const String& tmp) const;
    size_t rfind(const String& tmp) const;
    String substr(size_t s, size_t n) const;
    void clear();
    bool empty() const;
    friend std::istream& operator>>(std::istream& in, String& tmp);
    char* getstr() const {
        return str;
    }
};


std::istream& operator>>(std::istream& in, String& tmp) {
    tmp.clear();
    char c;
    while ((c = in.get()) != -1 && isspace(c)) {}
    tmp.push_back(c);
    while ((c = in.get()) != -1 && !isspace(c)) {
        tmp.push_back(c);
    }
    return in;
}


std::ostream& operator<<(std::ostream& out, const String& tmp) {
    out << tmp.getstr();
    return out;
}


bool String::empty() const {
    return len == 0;
}


void String::clear() {
    resize(0);
}


String String::substr(size_t s, size_t n) const {
    return String(this->str, s, n);
}


size_t String::rfind(const String& tmp) const {
    String self(*this);
    String tmpp(tmp);
    self.reverse();
    tmpp.reverse();
    if (strstr(self.str, tmpp.str) == nullptr) return len;
    return len - 1 - (strstr(self.str, tmpp.str) - self.str) - tmpp.len + 1;
}


size_t String::find(const String& tmp) const {
    if (strstr(str, tmp.str) == nullptr) return len;
    return strstr(str, tmp.str) - str;
}


void String::reverse() {
    for (size_t i = 0; i < len / 2; ++i) {
        std::swap(str[i], str[len - 1 - i]);
    }
}


String operator+(const String& tmp, const String& tmpp) {
    String ret = tmp;
    ret += tmpp;
    return ret;
}


String& String::operator+=(const String& tmp) {
    resize(len + tmp.len);
    memcpy(str + len - tmp.len, tmp.str, tmp.len);
    return *this;
}


const char& String::back() const {
    return str[len - 1];
}


const char& String::front() const {
    return *str;
}


char& String::back() {
    return str[len - 1];
}


char& String::front() {
    return *str;
}


void String::pop_back() {
    if (len) resize(len - 1);
}


void String::push_back(char c) {
    resize(len + 1);
    str[len - 1] = c;
    str[len] = '\0';
}


void String::resize(size_t n) {
    if (cap <= n) {
        cap = n * 2;
        char* tmp = new char[cap];
        memcpy(tmp, str, len);
        delete [] str;
        str = tmp;
    }
    len = n;
    str[len] = '\0';
}


size_t String::length() const {
    return len;
}


bool String::operator==(const String& tmp) const {
    return !strcmp(str, tmp.str);
}


char& String::operator[](size_t n) {
    return str[n];
}


const char& String::operator[](size_t n) const {
    return str[n];
}


String& String::operator=(const String& tmp) {
    String cpy(tmp);
    swap(cpy);
    return *this;
}


String::String(const String& tmp, size_t s, size_t n)
        : len(std::min(n, tmp.len - s)), cap(2 * len + 1), str(new char[cap]) {
    String cpy(tmp);
    memcpy(str, cpy.str + s, len);
    str[len] = '\0';
}


String::String(const String& tmp) : String(tmp.str) {}


String::String(size_t n, char c)
        : len(n), cap(len * 2 + 1), str(new char[cap]) {
    memset(str, c, n);
    str[n] = '\0';
}


String::String(const char* ptr)
        : len(strlen(ptr)), cap(len * 2 + 1), str(new char[cap]) {
    memcpy(str, ptr, len);
    str[len] = '\0';
}


String::String(char c) : len(1), cap(len * 2 + 1), str(new char[cap]) {
    str[0] = c;
    str[1] = '\0';
}


String::String() : len(0), cap(len * 2 + 1), str(new char[cap]) {
    *str = '\0';
}


String::~String() {
    delete [] str;
}
