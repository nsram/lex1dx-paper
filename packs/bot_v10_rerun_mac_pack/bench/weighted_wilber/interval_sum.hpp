#ifndef INTERVAL_SUM_GUARD
#define INTERVAL_SUM_GUARD

#include <cstdlib>
#include <vector>

template<typename T>
class interval_sum {
public:
    interval_sum();
    interval_sum(const std::vector<T> &points);
    // Weighted items: per-item weight W, sum S1, sum of squares S2 (same length).
    interval_sum(const std::vector<T> &w, const std::vector<T> &wx, const std::vector<T> &wx2);
    T query_w(size_t i, size_t j) const;
    T query_sq(size_t i, size_t j) const;
    T query(size_t i, size_t j) const;
    T cost_interval_l2(size_t i, size_t j) const;
    T mean(size_t i, size_t j) const;
    ~interval_sum();
private:
    std::vector<T> prefix_sum;
    std::vector<T> prefix_sum_of_squares;
    std::vector<T> prefix_w;
};

#endif
