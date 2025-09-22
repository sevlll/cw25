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

    [[nodiscard]] Interval sin() const {
        long double two_pi = 2 * PI;
        if (right - left >= two_pi)
            return Interval(-1., 1.);

        long double a = normalizeAngle(left);
        long double b = a + (right - left);
        if (b >= two_pi)
            return Interval(-1, 1); // safe

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
const Interval R = Interval(0.71330, 0.71331);
const Interval S = Interval(0.62745, 0.62746);
const Interval T = Interval(0.55622, 0.55623);


void one_coronas() {
    for (int len = 3; len <= 11; len++) {
        int MASK = (1 << (2 * len));
        for (int mask = 0; mask < MASK; mask++) {
            std::vector<Interval> seq(len);
            for (int i = 0; i < len; i++) {
                seq[i] = (std::vector<Interval>) {O, R, S, T}[(mask >> (2 * i)) & 3];
            }
            Interval circ_angle(0., 0.);
            for (int i = 0; i < len; i++) {
                circ_angle = circ_angle + angle(seq[i], O, seq[(i + 1) % len]);
            }
            if (circ_angle.left <= 2 * PI && 2 * PI <= circ_angle.right) {
                for (int i = 0; i < len; i++) {
                    std::cout << "1RST"[(mask >> (2 * i)) & 3];
                }
                std::cout << std::endl;
            }
        }
    }
}

void r_coronas() {
    for (int len = 3; len <= 11; len++) {
        int MASK = (1 << (2 * len));
        for (int mask = 0; mask < MASK; mask++) {
            std::vector<Interval> seq(len);
            for (int i = 0; i < len; i++) {
                seq[i] = (std::vector<Interval>) {O, R, S, T}[(mask >> (2 * i)) & 3];
            }
            Interval circ_angle(0., 0.);
            for (int i = 0; i < len; i++) {
                circ_angle = circ_angle + angle(seq[i], R, seq[(i + 1) % len]);
            }
            if (circ_angle.left <= 2 * PI && 2 * PI <= circ_angle.right) {
                for (int i = 0; i < len; i++) {
                    std::cout << "1RST"[(mask >> (2 * i)) & 3];
                }
                std::cout << std::endl;
            }
        }
    }
}

void s_coronas() {
    for (int len = 3; len <= 11; len++) {
        int MASK = (1 << (2 * len));
        for (int mask = 0; mask < MASK; mask++) {
            std::vector<Interval> seq(len);
            for (int i = 0; i < len; i++) {
                seq[i] = (std::vector<Interval>) {O, R, S, T}[(mask >> (2 * i)) & 3];
            }
            Interval circ_angle(0., 0.);
            for (int i = 0; i < len; i++) {
                circ_angle = circ_angle + angle(seq[i], S, seq[(i + 1) % len]);
            }
            if (circ_angle.left <= 2 * PI && 2 * PI <= circ_angle.right) {
                for (int i = 0; i < len; i++) {
                    std::cout << "1RST"[(mask >> (2 * i)) & 3];
                }
                std::cout << std::endl;
            }
        }
    }
}

void t_coronas() {
    for (int len = 3; len <= 11; len++) {
        int MASK = (1 << (2 * len));
        for (int mask = 0; mask < MASK; mask++) {
            std::vector<Interval> seq(len);
            for (int i = 0; i < len; i++) {
                seq[i] = (std::vector<Interval>) {O, R, S, T}[(mask >> (2 * i)) & 3];
            }
            Interval circ_angle(0., 0.);
            for (int i = 0; i < len; i++) {
                circ_angle = circ_angle + angle(seq[i], T, seq[(i + 1) % len]);
            }
            if (circ_angle.left <= 2 * PI && 2 * PI <= circ_angle.right) {
                for (int i = 0; i < len; i++) {
                    std::cout << "1RST"[(mask >> (2 * i)) & 3];
                }
                std::cout << std::endl;
            }
        }
    }
}

int main() {
    one_coronas();
    std::cout << std::endl;
    r_coronas();
    std::cout << std::endl;
    s_coronas();
    std::cout << std::endl;
    t_coronas();
    std::cout << std::endl;
    return 0;
}
