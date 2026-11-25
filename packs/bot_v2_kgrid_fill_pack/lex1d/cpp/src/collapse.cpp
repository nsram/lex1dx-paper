// Greedy adjacent merging (Ward / PNN restricted to neighbours): repeatedly merge the adjacent pair
// whose merge raises SSE least, until k bins remain. Equal costs merge the leftmost pair: the heap
// orders by (cost, id), and a tie means exactly equal doubles (no tolerance). This rule only fixes
// the start. Refine runs afterwards, and its tie sweeps move bounds to the largest cuts (as R's
// Ckmeans.1d.dp does), so a leftmost choice here does not decide where final ties land.
// O(d log d) with an indexed 4-ary heap keyed by each bin's cost of merging into its left neighbour.
#include "detail.hpp"

namespace lex::detail {
namespace {

/// 4-ary min-heap over ids 0 .. capacity-1, ordered by (key, id), with key updates. Each node
/// keeps its key and id together (16 bytes) and a node's 4 children are adjacent (64 bytes), which
/// halves the depth. The pop order is the (key, id) minimum, so it does not depend on the arity.
class MinHeap {
public:
    explicit MinHeap(int capacity) : nodes_(capacity), pos_(capacity, -1) {}

    /// Adds an entry without restoring heap order; call heapify() after the last one.
    void append(int id, double k) { put(id, k, size_++); }
    void heapify() {
        if (size_ < 2) return; // nothing to order (and (size_ - 2) / 4 would round to 0)
        for (int i = (size_ - 2) / 4; i >= 0; --i) down(i);
    }

    int pop() {
        const int top = nodes_[0].id;
        if (--size_ > 0) {
            put(nodes_[size_].id, nodes_[size_].key, 0);
            down(0);
        }
        pos_[top] = -1;
        return top;
    }

    void update(int id, double k) {
        const int i = pos_[id];
        const double old_k = nodes_[i].key;
        nodes_[i].key = k;
        if (k < old_k) up(i);
        else if (k > old_k) down(i);
        else if (k != k || old_k != old_k) { up(i); down(pos_[id]); } // NaN key: both directions
    }

private:
    struct Node { double key; int id; };
    std::vector<Node> nodes_; // key and id adjacent in the same cache line
    std::vector<int> pos_;    // heap position of id
    int size_ = 0;

    void put(int id, double k, int i) { nodes_[i] = {k, id}; pos_[id] = i; }
    static bool less(double ka, int a, double kb, int b) {
        return ka < kb || (ka == kb && a < b);
    }

    void up(int i) {
        const int id = nodes_[i].id;
        const double k = nodes_[i].key;
        while (i > 0) {
            const int p = (i - 1) >> 2;
            if (!less(k, id, nodes_[p].key, nodes_[p].id)) break;
            put(nodes_[p].id, nodes_[p].key, i);
            i = p;
        }
        put(id, k, i);
    }

    void down(int i) {
        const int id = nodes_[i].id;
        const double k = nodes_[i].key;
        for (;;) {
            const int l = (i << 2) + 1;
            if (l >= size_) break;
            int c = l;
            double ck = nodes_[l].key;
            int cid = nodes_[l].id;
            const int end = l + 4 < size_ ? l + 4 : size_;
            for (int j = l + 1; j < end; ++j) {
                if (less(nodes_[j].key, nodes_[j].id, ck, cid)) {
                    c = j; ck = nodes_[j].key; cid = nodes_[j].id;
                }
            }
            if (less(ck, cid, k, id)) {
                put(cid, ck, i);
                i = c;
            } else break;
        }
        put(id, k, i);
    }
};

} // namespace

std::vector<int> collapse(const Atoms& a, int k) {
    const int d = a.d;
    std::vector<int> prev(d), next(d); // next == d marks the end
    for (int j = 0; j < d; ++j) { prev[j] = j - 1; next[j] = j + 1; }
    auto cost = [&](int x) { return a.merge_cost(a.start[prev[x]], a.start[x], a.start[next[x]]); };

    MinHeap heap(d);
    for (int id = 1; id < d; ++id) heap.append(id, cost(id));
    heap.heapify();
    for (int bins = d; bins > k; --bins) {
        const int x = heap.pop(); // merge bin x into its left neighbour p
        const int p = prev[x], q = next[x];
        next[p] = q;
        if (q < d) prev[q] = p;
        if (p > 0) heap.update(p, cost(p));
        if (q < d) heap.update(q, cost(q));
    }

    std::vector<int> b(k + 1);
    int j = 0;
    for (int bin = 0; bin < d; bin = next[bin]) b[j++] = a.start[bin];
    b[k] = a.n;
    return b;
}

} // namespace lex::detail
