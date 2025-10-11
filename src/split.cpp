#include "split.h"

namespace packing {

    std::vector<Interval> split_in_parts(const Interval &q, int k) {
        std::vector<Interval> parts;
        parts.reserve(k);
        Interval step = (q.upper() - q.lower()) / static_cast<double>(k);
        for (int i = 0; i < k; ++i) {
            Interval left = q.lower() + step * static_cast<double>(i);
            Interval right = q.lower() + step * static_cast<double>(i + 1);
            parts.emplace_back(left.lower(), right.upper());
        }
        return parts;
    }

    std::vector<IntervalTrio> smart_split(IntervalTrio &rst) {
        std::vector<IntervalTrio> out;
        auto r = rst.r, s = rst.s, t = rst.t;
        double dr = r.upper() - r.lower();
        double ds = s.upper() - s.lower();
        double dt = t.upper() - t.lower();

        if (dr >= ds && dr >= dt) {
            double mid = (r.lower() + r.upper()) / 2.0;
            out.emplace_back(Interval(r.lower(), mid), s, t);
            out.emplace_back(Interval(mid, r.upper()), s, t);
        } else if (ds >= dr && ds >= dt) {
            double mid = (s.lower() + s.upper()) / 2.0;
            out.emplace_back(r, Interval(s.lower(), mid), t);
            out.emplace_back(r, Interval(mid, s.upper()), t);
        } else {
            double mid = (t.lower() + t.upper()) / 2.0;
            out.emplace_back(r, s, Interval(t.lower(), mid));
            out.emplace_back(r, s, Interval(mid, t.upper()));
        }
        return out;
    }

}