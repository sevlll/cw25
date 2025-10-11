#include "types.h"
#include "coronas.h"
#include "ub_search.h"
#include <iostream>
#include <vector>
#include <cstdlib>

using namespace packing;

int main(int argc, char *argv[]) {
    if (argc < 3) {
        return 1;
    }

    int i = std::atoi(argv[1]);
    int j = std::atoi(argv[2]);

    fill_all_valid_distributions();

    std::vector<Interval> fixed_intervals;
    for (double low = T_LB; low + GAP < 1.0 - GAP / 2.0; low += GAP) {
        fixed_intervals.push_back({low, low + GAP});
    }

    if (i < 0 || i >= static_cast<int>(fixed_intervals.size()) ||
        j < 0 || j >= static_cast<int>(fixed_intervals.size()) || i >= j) {
        return 2;
    }

    IntervalTrio W{
            Interval{fixed_intervals[j].upper(), 1.0 - GAP},
            fixed_intervals[j],
            fixed_intervals[i]
    };

    int ub = find_ub(W);
    std::cout << ub << "\n";
    return 0;
}
