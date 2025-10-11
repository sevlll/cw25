#include "ub_search.h"
#include "angle.h"
#include "coronas.h"
#include "split.h"
#include <iostream>
#include <set>
#include <map>
#include <algorithm>

namespace packing {

    struct CompareTrios {
        bool operator()(const std::pair<IntervalTrio, std::vector<std::vector<std::pair<int, int>>>> &a,
                        const std::pair<IntervalTrio, std::vector<std::vector<std::pair<int, int>>>> &b) const {
            return a.first.area() > b.first.area();
        }
    };

    static constexpr double AREA_LIM = 1e-14;

    int find_ub(IntervalTrio &rst_init) {
        std::multiset<std::pair<IntervalTrio, std::vector<std::vector<std::pair<int, int>>>>, CompareTrios> Q;
        Q.insert({rst_init, {}});

        int ub = 0;

        while (!Q.empty()) {
            auto cur = *Q.begin();
            Q.erase(Q.begin());

            IntervalTrio rst = cur.first;
            auto where = cur.second;

            if (!rst.gap_norm()) continue;
            if (cheap_check_if_done(rst)) continue;
            if (rst.area() < AREA_LIM) break;

            int sols1 = look_for_coronas(rst, 0, where);
            int solsR = look_for_coronas(rst, 1, where);
            int solsS = look_for_coronas(rst, 2, where);
            int solsT = look_for_coronas(rst, 3, where);

            if (std::min({sols1, solsR, solsS, solsT}) == 1) {
                continue;
            }

            auto children = smart_split(rst);
            for (const auto &ch: children) {
                Q.insert({ch, where});
            }
        }

        std::map<std::vector<std::vector<std::vector<std::vector<int>>>>, IntervalTrio> all_map;

        for (const auto &item: Q) {
            IntervalTrio rst = item.first;
            auto where = item.second;

            std::vector<std::vector<std::vector<std::vector<int>>>> all_coronas(4);
            for (int x = 1; x <= 3; ++x) {
                for (auto [disks, id]: where[x]) {
                    const auto &adj = ALL_VALID_DISTRIBUTIONS[disks][id];
                    if (adj[x][x] == 6) continue;
                    all_coronas[x].push_back(adj);
                }
            }

            int Rn = static_cast<int>(all_coronas[1].size());
            int Sn = static_cast<int>(all_coronas[2].size());
            int Tn = static_cast<int>(all_coronas[3].size());

            for (int mr = 1; mr < (1 << Rn); ++mr) {
                for (int ms = 1; ms < (1 << Sn); ++ms) {
                    for (int mt = 1; mt < (1 << Tn); ++mt) {
                        std::vector<std::vector<std::vector<std::vector<int>>>> cor(4);
                        for (int i = 0; i < Rn; ++i) if ((mr >> i) & 1) cor[1].push_back(all_coronas[1][i]);
                        for (int i = 0; i < Sn; ++i) if ((ms >> i) & 1) cor[2].push_back(all_coronas[2][i]);
                        for (int i = 0; i < Tn; ++i) if ((mt >> i) & 1) cor[3].push_back(all_coronas[3][i]);

                        bool bad = false;
                        for (int fmask = 1; fmask <= 7; ++fmask) {
                            bool fund = false;
                            for (int i = 1; i < 4 && !fund; ++i) {
                                if ((fmask >> (i - 1)) & 1) {
                                    for (int j = 0; j < 4 && !fund; ++j) {
                                        if (j == 0 || !((fmask >> (j - 1)) & 1)) {
                                            for (const auto &adj: cor[i]) {
                                                if (adj[std::min(i, j)][std::max(i, j)]) {
                                                    fund = true;
                                                    break;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            if (!fund) {
                                bad = true;
                                break;
                            }
                        }
                        if (bad) continue;

                        std::sort(cor[1].begin(), cor[1].end());
                        std::sort(cor[2].begin(), cor[2].end());
                        std::sort(cor[3].begin(), cor[3].end());

                        all_map[{cor[1], cor[2], cor[3]}] = rst;
                    }
                }
            }
        }

        ub += static_cast<int>(all_map.size());
        return ub;
    }

} 
