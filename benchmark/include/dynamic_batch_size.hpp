#pragma once

#include <algorithm>

// Calculates the optimal next batch size based on the current ESS
inline unsigned int compute_next_batch_size(
    unsigned int target_ESS, 
    unsigned int current_ESS, 
    size_t total_samples, 
    unsigned int current_batch_size)
{
    double ess_per_sample = static_cast<double>(current_ESS) / static_cast<double>(total_samples);
    
    if (ess_per_sample > 1e-6) {
        unsigned int remaining_ESS = target_ESS - current_ESS;
        
        if (current_ESS < (target_ESS * 0.85)) {
            // We are less than 85% to the target.
            // Since ESS generation usually improves over time, we intentionally 
            // ask for ONLY 80% of what the math predicts. This prevents massive overshoots.
            unsigned int samples_needed = static_cast<unsigned int>((remaining_ESS / ess_per_sample) * 0.80);
            return std::max(samples_needed, 500u);
        } else {
            // We are in the final 15%.
            // Now the ess_per_sample is highly accurate. We add a small 5% buffer 
            // to ensure we cross the finish line and don't waste time on tiny micro-batches.
            unsigned int samples_needed = static_cast<unsigned int>((remaining_ESS / ess_per_sample) * 1.05);
            return std::max(samples_needed, 500u);
        }
    } 
    
    // Fallback if ess_per_sample is virtually zero
    return static_cast<unsigned int>(current_batch_size * 1.5);
}