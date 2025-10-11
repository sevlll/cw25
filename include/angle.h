#pragma once

#include "types.h"
#include <vector>

namespace packing {

    Interval angle(const Interval &x, const Interval &y, const Interval &z);

    void sort_merge_adjacent_filter(std::vector<Interval> &qs);

    bool first_order_cut(IntervalTrio &rst, int mid_id);

    bool cheap_check_if_done(IntervalTrio &rst);

}
