// Tuned batch Collapse (q=0.75, heap finish at need <= 16384).
// Extracted from lex_stages.cpp so drivers can share it. Not part of lex1d.
#pragma once
#include "detail.hpp"
#include <algorithm>
#include <string>
#include <vector>

struct CollapseStats {
    int passes = 0;
    std::vector<int> merges_per_pass; // one entry per batch or heap-finish pass
    int merges_total() const {
        int s = 0;
        for (int m : merges_per_pass) s += m;
        return s;
    }
    std::string merges_csv() const {
        std::string o;
        for (std::size_t i = 0; i < merges_per_pass.size(); ++i) {
            if (i) o += ",";
            o += std::to_string(merges_per_pass[i]);
        }
        return o.empty() ? "0" : o;
    }
};

// Directional binary min-heap used only for the one-at-a-time finish of tuned batch.
class MinHeapDir {
public:
    explicit MinHeapDir(int capacity) : key_(capacity), heap_(capacity), pos_(capacity, -1) {}
    void append(int id, double k) { put(id, k, size_++); }
    void heapify() { for (int i = size_ / 2 - 1; i >= 0; --i) down(i); }
    int pop() {
        const int top = heap_[0];
        if (--size_ > 0) { put(heap_[size_], key_[size_], 0); down(0); }
        pos_[top] = -1;
        return top;
    }
    bool has(int id) const { return id >= 0 && id < static_cast<int>(pos_.size()) && pos_[id] >= 0; }
    int size() const { return size_; }
    void update(int id, double k) {
        const int i = pos_[id];
        const double old = key_[i];
        key_[i] = k;
        if (k < old) up(i);
        else if (k > old) down(i);
    }
private:
    std::vector<double> key_;
    std::vector<int> heap_;
    std::vector<int> pos_;
    int size_ = 0;
    void put(int id, double k, int i) { heap_[i] = id; key_[i] = k; pos_[id] = i; }
    static bool less(double ka, int a, double kb, int b) { return ka < kb || (ka == kb && a < b); }
    void up(int i) {
        const int id = heap_[i];
        const double k = key_[i];
        while (i > 0 && less(k, id, key_[(i - 1) / 2], heap_[(i - 1) / 2])) {
            const int p = (i - 1) / 2;
            put(heap_[p], key_[p], i);
            i = p;
        }
        put(id, k, i);
    }
    void down(int i) {
        const int id = heap_[i];
        const double k = key_[i];
        for (;;) {
            const int l = 2 * i + 1;
            const int c = l + 1 < size_ && less(key_[l + 1], heap_[l + 1], key_[l], heap_[l]) ? l + 1 : l;
            if (c < size_ && less(key_[c], heap_[c], k, id)) { put(heap_[c], key_[c], i); i = c; }
            else break;
        }
        put(id, k, i);
    }
};

// Tuned batch: 75th-percentile threshold, no extra need-th nth_element, reused
// pass buffers, and a real heap finish once (bins-target) <= 16384.
static std::vector<int> tuned_collapse(const lex::detail::Atoms& a, int target, CollapseStats* stats = nullptr) {
    const int d = a.d;
    std::vector<int> prev(d), next(d), ids;
    std::vector<double> costs, tmp;
    ids.reserve(static_cast<std::size_t>(d));
    costs.reserve(static_cast<std::size_t>(d));
    tmp.reserve(static_cast<std::size_t>(d));
    for (int j = 0; j < d; ++j) { prev[j] = j - 1; next[j] = j + 1; }
    auto cost = [&](int x) { return a.merge_cost(a.start[prev[x]], a.start[x], a.start[next[x]]); };

    auto heap_finish = [&](int& bins) {
        MinHeapDir heap(d);
        for (int x = next[0]; x < d; x = next[x]) heap.append(x, cost(x));
        heap.heapify();
        int merged = 0;
        while (bins > target && heap.size() > 0) {
            const int x = heap.pop();
            if (prev[x] < 0 || next[prev[x]] != x) continue;
            const int p = prev[x], q = next[x];
            next[p] = q;
            if (q < d) prev[q] = p;
            --bins;
            ++merged;
            if (p > 0 && heap.has(p)) heap.update(p, cost(p));
            if (q < d && heap.has(q)) heap.update(q, cost(q));
        }
        if (stats) { stats->passes += 1; stats->merges_per_pass.push_back(merged); }
    };

    int bins = d;
    const int cutover = 16384;
    const double q = 0.75;
    while (bins > target) {
        const int need = bins - target;
        if (need <= cutover) { heap_finish(bins); break; }
        ids.clear();
        costs.clear();
        for (int x = next[0]; x < d; x = next[x]) {
            ids.push_back(x);
            costs.push_back(cost(x));
        }
        if (ids.empty()) break;
        tmp.assign(costs.begin(), costs.end());
        const std::size_t qi = std::min(tmp.size() - 1,
            static_cast<std::size_t>(q * static_cast<double>(tmp.size())));
        std::nth_element(tmp.begin(), tmp.begin() + qi, tmp.end());
        const double thr = tmp[qi];
        int merged = 0;
        bool skip = false;
        for (std::size_t i = 0; i < ids.size() && bins > target; ++i) {
            if (skip) { skip = false; continue; }
            const int x = ids[i];
            if (prev[x] < 0 || next[prev[x]] != x) continue;
            if (costs[i] <= thr) {
                const int p = prev[x], qq = next[x];
                next[p] = qq;
                if (qq < d) prev[qq] = p;
                --bins;
                ++merged;
                skip = true;
            }
        }
        if (merged == 0) { heap_finish(bins); break; }
        if (stats) { stats->passes += 1; stats->merges_per_pass.push_back(merged); }
    }

    std::vector<int> b(target + 1);
    int j = 0;
    for (int bin = 0; bin < d && j < target; bin = next[bin]) b[j++] = a.start[bin];
    b[target] = a.n;
    return b;
}

