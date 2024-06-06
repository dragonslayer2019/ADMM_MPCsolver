#pragma once

#include "BlockMatrix.h"

template<typename T>
class Interval{
public:
    T m, M;

    inline void Prescaling(T x) {
        m *= x;
        M *= x;
    }

    inline T Project(T val) {
        if(val <= m) return m;
        if(val >= M) return M;
        return val;
    }

    T dis(T x) {
        if(x <= m) return m - x;
        if(x >= M) return x - M;
        return 0;
    }
};

template<typename T>
class QuadraticFunction{
public:
    T a, b, c; // coefficients
    T L, R; // Interval;
    
    QuadraticFunction() {}

    QuadraticFunction(T _a, T _b, T _c = 0, T _L = -1000, T _R = 1000) 
        : a(_a), b(_b), c(_c), L(_L), R(_R) {}

    inline T minimal() {
        if(a == 0) { //seems that we do not need to consider this case
            return (b >= 0) ? L : R;
        }
        return (-b) / (2 * a);
    }

    inline void Prescaling(T lamb) {
        a = a / (lamb * lamb);
        b = b / lamb;
        L *= lamb;
        R *= lamb;
    }

    QuadraticFunction<T> operator + (QuadraticFunction<T> other) {
        other.a += a; other.b += b; other.c += c;
        other.L = L; other.R = R;
        return other;
    }

    inline T Cost(T x) {
        return a * x * x + b * x + c;
    }
};

template<typename T>
class QuadraticSet{
public:
    int QuadraticNum;
    std::vector< QuadraticFunction<T> > q;

    inline void Initial(int N) {QuadraticNum = N; q.clear();}

    inline void Prescaling(T lamb) {
        for(int i = 0; i < QuadraticNum; ++i) {
            q[i].Prescaling(lamb);
        }
    }

    inline T Minimal(QuadraticFunction<T> h) {
        for(int i = 0; i < QuadraticNum; ++i) {
            T res = (q[i] + h).minimal();
            if(res <= q[i].R) return std::max(q[i].L, res);
        }
        return q[QuadraticNum - 1].R;
    }

    inline T Cost(T x) {
        for(int i = 0; i < QuadraticNum; ++i) {
            if(x >= q[i].L && x <= q[i].R) return q[i].Cost(x);
        }
        return 1e5;
    }

};

template<typename T>
class FunctionG {
public:
    bool IndicatorFlag;
    bool QuadraticFlag;
    Interval<T> I;
    QuadraticSet<T> Q;
    FunctionG () {IndicatorFlag = QuadraticFlag = 0;}

    inline void AddIndicator(T L, T R) {
        IndicatorFlag = 1;
        I.m = L, I.M = R;
    }

    inline void AddQuadratic(int QuadraticNum, std::vector<T> a, std::vector<T> b, std::vector<T> c, std::vector<T> p) {
        QuadraticFlag = 1;
        Q.Initial(QuadraticNum);
        for(int i = 0; i < QuadraticNum; ++i) {
            QuadraticFunction<T> qi;
            qi.a = a[i]; qi.b = b[i]; qi.c = c[i];
            qi.L = p[i]; qi.R = p[i + 1];
            Q.q.push_back(qi);
        }
    }

    inline void Prescaling(T lamb) {
        if(IndicatorFlag) I.Prescaling(lamb);
        if(QuadraticFlag) Q.Prescaling(lamb);
    }

    T DistanceOfIndicatorPart(T x) {
        if(IndicatorFlag) return I.dis(x);
        return 0;
    }

    T CostOfQuadraticPart(T x) {
        if(QuadraticFlag) return Q.Cost(x);
        return 0;
    }

    T Minimizer(T w, T rho) {
        T zstar = w;
        if(QuadraticFlag) {
            QuadraticFunction<T> h(rho/2, -1 * rho * w, 0);
            zstar = Q.Minimal(h);
        }
        if(IndicatorFlag) {
            zstar = I.Project(zstar);
        }
        return zstar;
    }
};
