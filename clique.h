#ifndef CLIQUE
#define CLIQUE

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

class CliqueCenter {
    public:
     int vertex_num;
     int k;
     std::vector<double> weight; // index starts from 0.
     std::vector<double> edge; // edge[i] means the length of v_i to v_((i + 1) % vertex_num).
     std::vector<double> loc; // location of v_i on circle starting form v_0 as 0.
    private:
     double perimeter;
     std::vector<double> center;
     std::vector<std::pair<int, int>> center_loc;
     double radius, radius_l, radius_r;
     double LQ_num;
     std::vector<bool> L;
     std::vector<std::pair<double, double>> intervals;
     std::vector<int> GD;
     std::vector<std::vector<int>> LQs;
     std::vector<int> next_vertex;
    public:
     CliqueCenter(const std::vector<double> &w, const std::vector<double> &e, int _k) {
        vertex_num = w.size();
        weight = w;
        edge = e;
        k = _k;
        loc.push_back(0);
        for (int i = 0; i < vertex_num - 1; i++)
            loc.push_back(loc.back() + edge[i]);
        perimeter = 0;
        for (auto len: edge)
            perimeter += len;
        radius_l = 0; radius_r = *std::max_element(weight.begin(), weight.end()) * (perimeter / 2);
        L.resize(vertex_num); next_vertex.resize(vertex_num);
     }
     void printCenter() {
        getCenter();
        std::cout << std::fixed;
        std::cout.precision(3);
        std::cout << "Radius: " << radius << "\n";
        for (int i = 0; i < k; i++) {
            std::cout << "Center " << i + 1 << " is between v" << center_loc[i].first << " and v" << center_loc[i].second 
                      << ", " << center[i] - loc[center_loc[i].first] << " away from v" << center_loc[i].first << ".\n";
        }
     }
    private:
     bool accuracyCheck(double l, double r, double epsilon) {
        return std::fabs(l - r) <= epsilon;
     }
     bool isOverlap(int idx1, int idx2) {
        std::pair<double, double> i1 = intervals[idx1], i2 = intervals[idx2];
        if (i1.first > i2.first) std::swap(i1, i2);
        return (i2.first < i1.second) || (i1.first + perimeter < i2.second);
     }
     double scale(double here) {
        if (here < 0) return here + perimeter;
        else if (here > perimeter) return here - perimeter;
        return here;
     }
     bool isBetween(double loc1, double loc_target, double loc2) {
        loc1 = scale(loc1);
        loc_target = scale(loc_target);
        loc2 = scale(loc2);
        if (loc1 > loc2) return (loc1 <= loc_target) || (loc_target <= loc2);
        else return (loc1 <= loc_target) && (loc_target <= loc2);
     }
     void makeIntervals() {
        intervals.clear();
        radius = (radius_l + radius_r) / 2;
        for (int i = 0; i < vertex_num; i++) {
            double left_end = std::max(loc[i] - radius / weight[i], loc[i] - perimeter / 2);
            double right_end = std::min(loc[i] + radius / weight[i], loc[i] + perimeter / 2);
            if (!intervals.empty() && left_end < intervals.back().first)
                continue;
            while (!intervals.empty() && right_end < intervals.back().second)
                intervals.pop_back();
            intervals.push_back({left_end, right_end});
        }
     }
     int NEXT(int idx) {
        int next_idx = idx + 1;
        while (isBetween(intervals[idx].second, intervals[next_idx % intervals.size()].first, intervals[next_idx % intervals.size()].second) == false) {
            next_idx++;
            if (next_idx % intervals.size() == idx)
                return idx;
        }
        return next_idx % intervals.size();
     }
     void GREEDY() {
        //1. Set L[i] = 0, for each arc i in F. Arbitrarily select an arc i in F;
        std::fill(L.begin(), L.end(), false);
        int i = 0, cur;

        //2. While L[i] = 0 do {L[i] <- 1; i <- NEXT(i)};
        while (L[i] == false) {
            L[i] = true;
            i = next_vertex[i] = NEXT(i);
        }

        //3. Find GD(i).
        cur = i; LQ_num = 0;
        GD.clear();
        GD.push_back(cur); LQ_num++;
        while (1) {
            cur = next_vertex[cur];
            if (isOverlap(cur, i)) 
                break;
            GD.push_back(cur); 
            LQ_num++; 
        }
        if (next_vertex[GD.back()] != i) 
            LQ_num++;
     }
     void makeLQ() {
        if (next_vertex[GD.back()] != GD.front()) 
            GD.push_back(next_vertex[GD.back()]);
        for (int i: GD) {
            std::vector<int> LQ;
            LQ.push_back(i);
            int idx = i + 1;
            while (isBetween(intervals[idx % intervals.size()].first, intervals[i].second, intervals[idx % intervals.size()].second))
                LQ.push_back(idx++ % intervals.size());
            LQs.push_back(LQ);
        }
     }
     std::pair<double, double> getOverlap(std::vector<int> &LQ) {
        std::vector<std::pair<double, double>> overlaps;
        overlaps.push_back({0, perimeter});
        for (auto LQ_idx: LQ) {
            std::vector<std::pair<double, double>> LQ_intervals;
            double left = scale(intervals[LQ_idx].first);
            double right = scale(intervals[LQ_idx].second);
            if (left > right) {
                LQ_intervals.push_back({left, perimeter});
                LQ_intervals.push_back({0, right});
            }
            else LQ_intervals.push_back({left, right});
            std::vector<std::pair<double, double>> tmp_overlaps = overlaps;
            overlaps.clear();
            for (auto overlap: tmp_overlaps) {
                for (auto LQ_interval: LQ_intervals) {
                    double oleft = std::max(LQ_interval.first, overlap.first);
                    double oright = std::min(LQ_interval.second, overlap.second);
                    if (oleft < oright) overlaps.push_back({oleft, oright});
                }
            }
        }
        if (overlaps.size() == 1) return overlaps[0];
        else return std::make_pair(std::max(overlaps[0].first, overlaps[1].first),
                                   std::min(overlaps[0].second, overlaps[1].second));
     }
     void FindCenterLoc() {
        auto vertex_itr = loc.begin();
        for (auto LQ: LQs) {
            std::pair<double, double> range = getOverlap(LQ);
            double _center = (range.first + range.second) / 2;
            center.push_back(_center);
            while (vertex_itr != loc.end() && *vertex_itr < _center)
                vertex_itr++;
            int idx = vertex_itr - loc.begin();
            center_loc.push_back({idx - 1, idx % vertex_num}); 
        }
     }
     void getCenter() {
        while (LQ_num > k || accuracyCheck(radius_l, radius_r, 0.0001) == false) {
            makeIntervals();
            GREEDY();
            if (LQ_num <= k) radius_r = radius;
            else radius_l = radius;
        }
        makeLQ();
        FindCenterLoc();
     }   
};

#endif
