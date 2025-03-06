#ifndef TINYNURBS_INTEGRAL_H
#define TINYNURBS_INTEGRAL_H

template <typename T> T integral(T (*f)(T), T a, T b, int N) {
    T h = (b - a) / N;
    T s = f(a) + f(b);
    for (int i = 1; i <= N - 1; ++i) {
        T x = a + h * i;
        s += f(x) * ((i & 1) ? 4 : 2);
    }
    s *= h / 3;
    return s;
}

template <typename T> T integral(T *f, T a, T b, int N) {
    T h = (b - a) / N;
    T s = f[0] + f[N];
    for (int i = 1; i <= N - 1; ++i) {
        s += f[i] * ((i & 1) ? 4 : 2);
    }
    s *= h / 3;
    return s;
}

#endif // TINYNURBS_INTEGRAL_H
