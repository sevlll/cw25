#include <bits/stdc++.h>

const long double PI = acosl(-1);

class Interval {
public:
    long double left = 0.;
    long double right = 0.;

    Interval() = default;

    Interval(long double l, long double r) : left(l), right(r) {
        if (l > r) {
            throw std::invalid_argument("GG l > r");
        }
    }

    Interval operator+(const Interval &other) const {
        return Interval(left + other.left, right + other.right);
    }

    Interval operator-(const Interval &other) const {
        return Interval(left - other.right, right - other.left);
    }

    Interval operator*(const Interval &other) const {
        long double p1 = left * other.left;
        long double p2 = left * other.right;
        long double p3 = right * other.left;
        long double p4 = right * other.right;
        long double new_lo = std::min({p1, p2, p3, p4});
        long double new_hi = std::max({p1, p2, p3, p4});
        return Interval(new_lo, new_hi);
    }

    Interval operator/(const Interval &other) const {
        if (other.left <= 0 && other.right >= 0) {
            throw std::invalid_argument("Ошибка: деление на интервал, содержащий 0, не определено.");
        }
        long double r1 = 1.0L / other.left;
        long double r2 = 1.0L / other.right;
        Interval reciprocal(std::min(r1, r2), std::max(r1, r2));
        return (*this) * reciprocal;
    }

    static long double normalizeAngle(long double angle) {
        long double two_pi = 2 * PI;
        angle = fmodl(angle, two_pi);
        if (angle < 0) {
            angle += two_pi;
        }
        return angle;
    }

    Interval sqrt() const {
        if (left < 0) {
            throw std::domain_error("GANDON");
        }
        long double lo = sqrtl(left);
        long double hi = sqrtl(right);
        return Interval(lo, hi);
    }


    [[nodiscard]] Interval sin() const {
        long double two_pi = 2 * PI;
        if (right - left >= two_pi)
            return Interval(-1., 1.);

        long double a = normalizeAngle(left);
        long double b = a + (right - left);
        if (b >= two_pi)
            return Interval(-1, 1);

        long double sin_a = sinl(a);
        long double sin_b = sinl(b);
        long double min_val = std::min(sin_a, sin_b);
        long double max_val = std::max(sin_a, sin_b);

        if (a <= PI / 2 && b >= PI / 2)
            max_val = 1;
        if (a <= (3 * PI) / 2 && b >= (3 * PI) / 2)
            min_val = -1;
        return Interval(min_val, max_val);
    }

    Interval cos() const {
        long double two_pi = 2 * PI;
        if (right - left >= two_pi)
            return Interval(-1, 1);

        long double a = normalizeAngle(left);
        long double b = a + (right - left);
        if (b >= two_pi)
            return Interval(-1, 1);

        long double cos_a = cosl(a);
        long double cos_b = cosl(b);
        long double min_val = std::min(cos_a, cos_b);
        long double max_val = std::max(cos_a, cos_b);

        if (a == 0 || b == two_pi)
            max_val = 1;
        if (a <= PI && b >= PI)
            min_val = -1;
        return Interval(min_val, max_val);
    }

    Interval operator-() const {
        return Interval(-right, -left);
    }

    Interval tan() const {
        long double period = PI;
        if (right - left >= period)
            return Interval(-std::numeric_limits<long double>::infinity(),
                            std::numeric_limits<long double>::infinity());
        long double a = fmodl(left, period);
        if (a < 0)
            a += period;
        long double b = a + (right - left);
        if (b >= period)
            return Interval(-std::numeric_limits<long double>::infinity(),
                            std::numeric_limits<long double>::infinity());
        if (a <= period / 2 && b >= period / 2)
            return Interval(-std::numeric_limits<long double>::infinity(),
                            std::numeric_limits<long double>::infinity());
        long double tan_a = tanl(a);
        long double tan_b = tanl(b);
        long double min_val = std::min(tan_a, tan_b);
        long double max_val = std::max(tan_a, tan_b);
        return Interval(min_val, max_val);
    }

    Interval arccos() const {
        if (left < -1 || right > 1) {
            throw std::domain_error("Ошибка: arccos определен только на интервале [-1, 1].");
        }
        return Interval(acosl(right), acosl(left));
    }
};

std::ostream &operator<<(std::ostream &os, const Interval &iv) {
    os << "[" << iv.left << ", " << iv.right << "]";
    return os;
}

Interval angle(const Interval &x, const Interval &y, const Interval &z) {
    Interval numerator = (y * y) + (x * y) + (y * z) - (x * z);

    Interval denominator = (x + y) * (y + z);

    Interval fraction = numerator / denominator;

    return fraction.arccos();
}

const Interval O = Interval(1., 1.);
const Interval R = Interval(0.71331, 0.71332);
const Interval S = Interval(0.62746, 0.62747);
const Interval T = Interval(0.55623, 0.55624);


Interval check1() {
    Interval one(1.0, 1.0), two(2.0, 2.0), three(3.0, 3.0), four(4.0, 4.0);

    Interval R2 = R * R;
    Interval T2 = T * T;
    Interval RT = R * T;

    Interval Rp1 = R + one;
    Interval Tp1 = T + one;
    Interval Tp1_sq = Tp1 * Tp1;

    Interval big = R2 * T2
                   + two * R2 * T
                   + R2
                   + two * R * T2
                   + four * R * T
                   + two * R
                   + T2
                   + two * T
                   + one;

    Interval term1 = -three.sqrt()
                     * R.sqrt()
                     * T.sqrt()
                     * (R + T + one).sqrt()
                     / (Tp1 * big.sqrt());

    Interval term2 = -R.sqrt()
                     * T
                     * (T + two).sqrt()
                     * (R + T + one).sqrt()
                     / ((T2 + two * T + one).sqrt() * (big).sqrt());

    Interval term3 = -(three).sqrt()
                     * (T).sqrt()
                     * (T + two).sqrt()
                     * (-RT + R + T + one)
                     / (two * Rp1 * Tp1 * (T2 + two * T + one).sqrt());

    Interval term4 = (
            one / two
            + (-R2 + two * R + one)
              / (two * Rp1 * Rp1)
    ).sqrt();

    Interval term5 = (-RT + R + T + one)
                     / (two * Rp1 * Tp1_sq);

    return term1 + term2 + term3 + term4 + term5;
}

Interval check2() {
    Interval one(1.0L, 1.0L),
            two(2.0L, 2.0L),
            four(4.0L, 4.0L);

    Interval R2 = R * R;
    Interval R3 = R2 * R;
    Interval R4 = R2 * R2;
    Interval S2 = S * S;
    Interval T2 = T * T;

    Interval Rp1 = R + one;
    Interval RpT = R + T;
    Interval RpS = R + S;
    Interval RS = R * S;
    Interval RT = R * T;

    Interval A = R2 + two * R + one;

    Interval big1 = R4
                    + two * R3 * S
                    + two * R3 * T
                    + R2 * S2
                    + four * R2 * S * T
                    + R2 * T2
                    + two * R * S2 * T
                    + two * R * S * T2
                    + S2 * T2;

    Interval big2 = R4
                    + two * R3 * T
                    + two * R3
                    + R2 * T2
                    + four * R2 * T
                    + R2
                    + two * R * T2
                    + two * R * T
                    + T2;

    Interval frac = two * T / (Rp1 * RpT);
    Interval factor1 = frac - one;

    Interval term1 = two
                     * R.sqrt()
                     * S.sqrt()
                     * T.sqrt()
                     * (two * R + one).sqrt()
                     * factor1
                     * (R + S + T).sqrt()
                     / (A.sqrt() * big1.sqrt());

    Interval numer2 = two
                      * R.sqrt()
                      * T.sqrt()
                      * (two * R + one).sqrt()
                      * (RpT + one).sqrt()
                      * (R2 + R * (S + T) - S * T);
    Interval term2 = -numer2 / (RpS * RpT * A.sqrt() * big2.sqrt());

    Interval numer3 = four
                      * R2
                      * S.sqrt()
                      * T
                      * (R + S + T).sqrt()
                      * (RpT + one).sqrt();
    Interval term3 = -numer3 / (Rp1 * big1.sqrt() * big2.sqrt());

    Interval numer4 = R * factor1 * (R2 + R * (S + T) - S * T);
    Interval term4 = -numer4 / (Rp1 * RpS * RpT);

    Interval term5 = one;

    return term1 + term2 + term3 + term4 + term5;
}

Interval check3() {
    Interval one(1.0L, 1.0L),
            two(2.0L, 2.0L),
            four(4.0L, 4.0L);

    Interval R2 = R * R;
    Interval S2 = S * S;
    Interval T2 = T * T;
    Interval T3 = T2 * T;
    Interval T4 = T2 * T2;

    Interval RS = R * S;
    Interval RpT = R + T;
    Interval SpT = S + T;
    Interval Tp1 = T + one;

    Interval big1 = R2 * S2
                    + two * R2 * S * T
                    + R2 * T2
                    + two * R * S2 * T
                    + four * R * S * T2
                    + two * R * T3
                    + S2 * T2
                    + two * S * T3
                    + T4;

    Interval big2 = R2 * T2
                    + two * R2 * T
                    + R2
                    + two * R * T3
                    + four * R * T2
                    + two * R * T
                    + T4
                    + two * T3
                    + T2;

    Interval term1 = -four
                     * R
                     * S.sqrt()
                     * T
                     * (R + S + T).sqrt()
                     * (R + T + one).sqrt()
                     / (big1.sqrt() * big2.sqrt());

    Interval term2 = (one - one / (Tp1 * Tp1)).sqrt();

    Interval frac3 = -two * R / (RpT * Tp1) + one;
    Interval numer3 = -RS + T2 + T * (R + S);
    Interval denom3 = RpT * SpT;
    Interval term3 = frac3 * numer3 / denom3;

    return term1 + term2 + term3;
}


int main() {
    std::cout << std::fixed << std::setprecision(20);
    std::cout << check1() << '\n';
    std::cout << check2() << '\n';
    std::cout << check3() << '\n';

    return 0;
}
