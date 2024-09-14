// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2024 Luca Perju

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef ACCELERATED_BILLIARD_WALK_UTILS_HPP
#define ACCELERATED_BILLIARD_WALK_UTILS_HPP

#include <Eigen/Eigen>
#include <vector>

const double eps = 1e-10;

// data structure which maintains the values of (b - Ar)/Av, and can extract the minimum positive value and the facet associated with it
// vec[i].first contains the value of (b(i) - Ar(i))/Av(i) + moved_dist, where moved_dist is the total distance that the point has travelled so far
// The heap will only contain the values from vec which are greater than moved_dist (so that they are positive)
template<typename NT>
class BoundaryOracleHeap {
public:
    int n, heap_size;
    std::vector<std::pair<NT, int>> heap;
    std::vector<std::pair<NT, int>> vec;

private:
    int siftDown(int index) {
        while((index << 1) + 1 < heap_size) {
            int child = (index << 1) + 1;
            if(child + 1 < heap_size && heap[child + 1].first < heap[child].first - eps) {
                child += 1;
            }
            if(heap[child].first < heap[index].first - eps)
            {
                std::swap(heap[child], heap[index]);
                std::swap(vec[heap[child].second].second, vec[heap[index].second].second);
                index = child;
            } else {
                return index;
            }
        }
        return index;
    }

    int siftUp(int index) {
        while(index > 0 && heap[(index - 1) >> 1].first - eps > heap[index].first) {
            std::swap(heap[(index - 1) >> 1], heap[index]);
            std::swap(vec[heap[(index - 1) >> 1].second].second, vec[heap[index].second].second);
            index = (index - 1) >> 1;
        }
        return index;
    }

    // takes the index of a facet, and (in case it is in the heap) removes it from the heap.
    void remove (int index) {
        index = vec[index].second;
        if(index == -1) {
            return;
        }
        std::swap(heap[heap_size - 1], heap[index]);
        std::swap(vec[heap[heap_size - 1].second].second, vec[heap[index].second].second);
        vec[heap[heap_size - 1].second].second = -1;
        heap_size -= 1;
        index = siftDown(index);
        siftUp(index);
    }

    // inserts a new value into the heap, with its associated facet
    void insert (const std::pair<NT, int> val) {
        vec[val.second].second = heap_size;
        vec[val.second].first = val.first;
        heap[heap_size++] = val;
        siftUp(heap_size - 1);
    }

public:
    BoundaryOracleHeap() {}

    BoundaryOracleHeap(int n) : n(n), heap_size(0) {
        heap.resize(n);
        vec.resize(n);
    }

    // rebuilds the heap with the existing values from vec
    // O(n)
    void rebuild (const NT &moved_dist) {
        heap_size = 0;
        for(int i = 0; i < n; ++i) {
            vec[i].second = -1;
            if(vec[i].first - eps > moved_dist) {
                vec[i].second = heap_size;
                heap[heap_size++] = {vec[i].first, i};
            }
        }
        for(int i = heap_size - 1; i >= 0; --i) {
            siftDown(i);
        }
    }

    // returns (b(i) - Ar(i))/Av(i) + moved_dist
    // O(1)
    NT get_val (const int &index) {
        return vec[index].first;
    }

    // returns the nearest facet
    // O(1)
    std::pair<NT, int> get_min () {
        return heap[0];
    }

    // changes the stored value for a given facet, and updates the heap accordingly
    // O(logn)
    void change_val(const int& index, const NT& new_val, const NT& moved_dist) {
        if(new_val < moved_dist - eps) {
            vec[index].first = new_val;
            remove(index);
        } else {
            if(vec[index].second == -1) {
                insert({new_val, index});
            } else {
                heap[vec[index].second].first = new_val;
                vec[index].first = new_val;
                siftDown(vec[index].second);
                siftUp(vec[index].second);
            }
        }
    }
};


#endif
