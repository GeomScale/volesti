// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2020 Apostolos Chalkis

// Contributed and/or modified by Repouskos Panagiotis, as part of Google Summer of Code 2019 program.
// Contributed and/or modified by Korakitis Angelos, as part of Google Summer of Code 2025 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef VOLESTI_SLIDING_WINDOW_HPP
#define VOLESTI_SLIDING_WINDOW_HPP

#include <list>
#include <cmath>

/// Sliding window for tracking convergence based on relative error
/// Maintains a fixed-size window of recent values to detect when optimization
/// has converged by comparing the oldest and newest entries.
/// \tparam NT Numeric type for stored values
template<typename NT>
class SlidingWindow {
public:
    /// Stored values (newest at front, oldest at back)
    std::list<NT> values;
    
    /// Maximum number of values to store in the window
    int windowSize;
    
    /// Current number of entries in the window
    int numEntries;
    
    /// Construct a sliding window with specified size
    /// \param[in] windowSize Maximum number of values to store
    SlidingWindow(int windowSize) : windowSize(windowSize), numEntries(0) {}
    
    /// Add a new value to the window
    /// Adds value to the front of the window. If window is full,
    /// removes the oldest value from the back.
    /// \param[in] value The new value to add
    void push(NT value) {
        // If window is full, remove the oldest value
        if (numEntries >= windowSize) {
            values.pop_back();
        } else {
            numEntries++;
        }
        values.push_front(value);
    }
    
    /// Calculate relative error between newest and oldest values
    /// Computes |oldest - newest| / |oldest| to measure convergence.
    /// Returns 1.0 if window not full or to avoid division by zero.
    /// \return Relative error in range [0, infinity), or 1.0 if not converged
    NT getRelativeError() const {
        if (numEntries < windowSize) {
            return NT(1); // Not converged yet
        }
        NT newest = values.front();
        NT oldest = values.back();
        if (std::abs(oldest) < NT(1e-10)) {
            return NT(1); // Avoid division by zero
        }
        return std::abs((oldest - newest) / oldest);
    }
    
    /// Check if window is full
    /// \return True if window has collected windowSize entries, false otherwise
    bool isFull() const {
        return numEntries >= windowSize;
    }
    
    /// Clear all values from the window
    /// Resets the window to empty state.
    void clear() {
        values.clear();
        numEntries = 0;
    }
};

#endif // VOLESTI_SLIDING_WINDOW_HPP