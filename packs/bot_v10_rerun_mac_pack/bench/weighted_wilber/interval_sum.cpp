#include "interval_sum.hpp"
#include <vector>
#include <cstdint>
#include <iostream>
#include <cassert>

template<typename T>
interval_sum<T>::interval_sum(const std::vector<T> &points) : prefix_sum(points.size() + 1),
                                                              prefix_sum_of_squares(points.size() + 1),
                                                              prefix_w(points.size() + 1) {
    T prefix_sum_val = 0;
    T prefix_square_val = 0;
    size_t n = points.size();
    this->prefix_sum[0] = 0;
    this->prefix_sum_of_squares[0] = 0;
    this->prefix_w[0] = 0;
    for (size_t i = 1; i <= n; ++i) {
        prefix_sum_val += points[i-1];
        prefix_square_val += points[i-1]*points[i-1];
        this->prefix_sum[i] = prefix_sum_val;
        this->prefix_sum_of_squares[i] = prefix_square_val;
        this->prefix_w[i] = static_cast<T>(i);
    }
}

template<typename T>
interval_sum<T>::interval_sum(const std::vector<T> &w, const std::vector<T> &wx, const std::vector<T> &wx2)
    : prefix_sum(w.size() + 1), prefix_sum_of_squares(w.size() + 1), prefix_w(w.size() + 1) {
    T W = 0, S = 0, Q = 0;
    size_t n = w.size();
    this->prefix_sum[0] = 0;
    this->prefix_sum_of_squares[0] = 0;
    this->prefix_w[0] = 0;
    for (size_t i = 1; i <= n; ++i) {
        W += w[i-1];
        S += wx[i-1];
        Q += wx2[i-1];
        this->prefix_w[i] = W;
        this->prefix_sum[i] = S;
        this->prefix_sum_of_squares[i] = Q;
    }
}

template<typename T>
interval_sum<T>::interval_sum() {}

template<typename T>
interval_sum<T>::~interval_sum() {}

/**
 * @returns sum from i to j of points[l]^2 where i included, j not included.
 */
template<typename T>
T interval_sum<T>::query_sq(size_t i, size_t j) const {
    return this->prefix_sum_of_squares[j] - this->prefix_sum_of_squares[i];
}

/**
 * @returns sum from i to j of points[l] where i included, j not included.
 */
template<typename T>
T interval_sum<T>::query(size_t i, size_t j) const {
    return this->prefix_sum[j] - this->prefix_sum[i];
}

template<typename T>
T interval_sum<T>::query_w(size_t i, size_t j) const {
    return this->prefix_w[j] - this->prefix_w[i];
}

/**
 * returns the cost of clustering points [start, i] into one cluster with L_2
 * norm.
 */
template<typename T>
T interval_sum<T>::cost_interval_l2(size_t start, size_t i) const {
  if(i <= start){ return 0.0;}
  T S = this->query(start, i+1);
  T S2 = this->query_sq(start, i+1);
  T W = this->query_w(start, i+1);
  if (W == 0) return 0.0;
  T interval_cost = S2 - (S * S) / W;
  if(interval_cost < 0){
    // numerical issues
    return 0.0;
  }
  return interval_cost;
}

/**
 * Returns the mean of points[i..j], both included.
 */
template<typename T>
T interval_sum<T>::mean(size_t i, size_t j) const {
    T W = this->query_w(i, j+1);
    if (W == 0) return 0;
    return this->query(i, j+1) / W;
}

template class interval_sum<double>;
template class interval_sum<std::int64_t>;
