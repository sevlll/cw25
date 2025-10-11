#include "coronas.h"
#include "angle.h"
#include <algorithm>

namespace packing {

    std::vector<std::vector<std::vector<std::vector<int>>>> ALL_VALID_DISTRIBUTIONS;

    static bool is_connected(std::vector<std::vector<int>> &adj) {
        std::vector<bool> used(4, false);
        int components = 0;
        for (int v = 0; v < 4; ++v) {
            int deg = 0;
            for (int i = 0; i < 4; ++i) deg += adj[std::min(i, v)][std::max(i, v)];
            if (deg == 0) continue;
            if (used[v]) continue;
            used[v] = true;
            std::vector<int> q{v};
            while (!q.empty()) {
                int u = q.back();
                q.pop_back();
                for (int w = 0; w < 4; ++w) {
                    if (!used[w] && adj[std::min(u, w)][std::max(u, w)]) {
                        used[w] = true;
                        q.push_back(w);
                    }
                }
            }
            components++;
        }
        return components <= 1;
    }

    static void brute(int i, int j, int sum, std::vector<std::vector<int>> &a) {
        if (i >= 4) {
            if (sum >= 3 && is_connected(a)) {
                ALL_VALID_DISTRIBUTIONS[sum].push_back(a);
            }
            return;
        }
        for (int add = 0; sum + add <= MAX_DISKS; ++add) {
            a[i][j] = add;
            int ni = i, nj = j;
            if (nj == 3) {
                ni++;
                nj = ni;
                int deg_i = 0;
                for (int prev = 0; prev < i; ++prev) deg_i += a[prev][i];
                for (int nxt = i + 1; nxt < 4; ++nxt) deg_i += a[i][nxt];
                if (deg_i % 2 == 0) {
                    brute(ni, nj, sum + add, a);
                }
            } else {
                nj++;
                brute(ni, nj, sum + add, a);
            }
        }
        a[i][j] = 0;
    }

    void fill_all_valid_distributions() {
        ALL_VALID_DISTRIBUTIONS.assign(MAX_DISKS + 1, {});
        std::vector<std::vector<int>> a(4, std::vector<int>(4, 0));
        brute(0, 0, 0, a);
    }

    int look_for_coronas(IntervalTrio &rst, int center_id,
                         std::vector<std::vector<std::pair<int, int>>> &where_to_look) {
        if (where_to_look.empty()) {
            where_to_look.resize(4);
            for (int d = 0; d < 4; ++d) {
                for (int disks = 3; disks <= MAX_DISKS; ++disks) {
                    const auto &bucket = ALL_VALID_DISTRIBUTIONS[disks];
                    for (int id = 0; id < static_cast<int>(bucket.size()); ++id) {
                        where_to_look[d].push_back({disks, id});
                    }
                }
            }
        }

        auto r = rst.r, s = rst.s, t = rst.t;
        std::vector<Interval> v = {ONE, r, s, t};
        Interval center = v[center_id];

        int solutions = 0;
        std::vector<std::pair<int, int>> filtered;

        for (auto [disks, id]: where_to_look[center_id]) {
            const auto &adj = ALL_VALID_DISTRIBUTIONS[disks][id];
            Interval sum{0.0, 0.0};
            for (int i = 0; i < 4; ++i) {
                for (int j = i; j < 4; ++j) {
                    int cnt = adj[i][j];
                    if (cnt == 0) continue;
                    sum += Interval{static_cast<double>(cnt), static_cast<double>(cnt)} *
                           angle(v[i], center, v[j]);
                }
            }
            if (sum.lower() <= 2 * PI.lower() && 2 * PI.lower() <= sum.upper()) {
                solutions++;
                filtered.push_back({disks, id});
            }
        }

        where_to_look[center_id] = std::move(filtered);
        return solutions;
    }

}
