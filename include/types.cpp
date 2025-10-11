#include "types.h"
#include <algorithm>

namespace packing {

    const Interval PI(boost::math::constants::pi<double>());
    const Interval ONE{1.0, 1.0};

    bool IntervalTrio::gap_norm() {
        if (std::max(T_LB, t.lower()) > t.upper()) return false;
        Interval tnew{std::max(T_LB, t.lower()), t.upper()};

        if (std::max(tnew.lower() + GAP, s.lower()) > s.upper()) return false;
        Interval snew{std::max(tnew.lower() + GAP, s.lower()), s.upper()};

        if (std::max(snew.lower() + GAP, r.lower()) > std::min(r.upper(), 1.0 - GAP)) return false;
        Interval rnew{std::max(snew.lower() + GAP, r.lower()), std::min(r.upper(), 1.0 - GAP)};

        r = rnew;
        s = snew;
        t = tnew;

        if (s.lower() > std::min(s.upper(), r.upper() - GAP)) return false;
        snew = Interval{s.lower(), std::min(s.upper(), r.upper() - GAP)};

        if (t.lower() > std::min(t.upper(), snew.upper() - GAP)) return false;
        tnew = Interval{t.lower(), std::min(t.upper(), snew.upper() - GAP)};

        s = snew;
        t = tnew;
        return true;
    }

    double IntervalTrio::area() const {
        return (r.upper() - r.lower()) * (s.upper() - s.lower()) * (t.upper() - t.lower());
    }

}