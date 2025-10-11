#include "angle.h"
#include <algorithm>

namespace packing {

    Interval angle(const Interval &x, const Interval &y, const Interval &z) {
        const Interval num = Interval{2.0, 2.0} * x * z;
        const Interval den = y * y + y * (x + z) + x * z;
        Interval q = Interval{1.0, 1.0} - num / den;
        double lo = std::max(q.lower(), -1.0);
        double hi = std::min(q.upper(), 1.0);
        const Interval q_cap(lo, hi);
        return acos(q_cap);
    }

    void sort_merge_adjacent_filter(std::vector<Interval> &q) {
        std::sort(q.begin(), q.end(), [](const Interval &a, const Interval &b) {
            return a.lower() < b.lower();
        });
        std::vector<Interval> merged;
        bool first = true;
        Interval cur{0.0, 0.0};
        for (const auto &intr: q) {
            if (first) {
                cur = intr;
                first = false;
                continue;
            }
            if (intr.lower() > cur.upper()) {
                if (cur.lower() <= 2 * PI.lower()) {
                    merged.push_back(cur);
                }
                cur = intr;
            } else {
                cur = Interval(cur.lower(), std::max(cur.upper(), intr.upper()));
            }
        }
        if (!first && cur.lower() <= 2 * PI.lower()) {
            merged.push_back(cur);
        }
        q.swap(merged);
    }

    static constexpr int LOCAL_MAX_DISKS = 7;

    bool first_order_cut(IntervalTrio &rst, int mid_id) {
        auto r = rst.r, s = rst.s, t = rst.t;
        std::vector<Interval> v = {ONE, r, s, t};
        Interval mid = v[mid_id];

        for (int start = 0; start < 4; ++start) {
            std::vector<std::vector<std::vector<Interval>>> total(LOCAL_MAX_DISKS + 1,
                                                                  std::vector<std::vector<Interval>>(4));
            total[0][start].push_back(Interval{0.0, 0.0});

            for (int level = 0; level < LOCAL_MAX_DISKS; ++level) {
                for (int j = start; j < 4; ++j) {
                    for (int i = start; i < 4; ++i) {
                        if (start == mid_id && level == 0 && j == mid_id) {
                            continue;
                        }
                        Interval w = angle(v[i], mid, v[j]);
                        for (const auto &prev: total[level][i]) {
                            total[level + 1][j].push_back(prev + w);
                        }
                    }
                    sort_merge_adjacent_filter(total[level + 1][j]);
                }
                for (const auto &intr: total[level + 1][start]) {
                    if (intr.lower() <= 2 * PI.lower() && 2 * PI.lower() <= intr.upper()) {
                        return true;
                    }
                }
            }
        }
        return false;
    }

    bool cheap_check_if_done(IntervalTrio &rst) {
        if (!first_order_cut(rst, 0)) return true;
        if (!first_order_cut(rst, 1)) return true;
        if (!first_order_cut(rst, 2)) return true;
        if (!first_order_cut(rst, 3)) return true;
        return false;
    }

}
