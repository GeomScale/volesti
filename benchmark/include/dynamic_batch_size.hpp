#pragma once

#include <algorithm>

// Calculates the optimal next batch size based on the current ESS
inline unsigned int compute_next_batch_size(
    unsigned int target_ESS, 
    unsigned int current_ESS, 
    size_t total_samples, 
    unsigned int current_batch_size)
{
    const unsigned int MAX_BATCH_SIZE = 500000;

    double ess_per_sample = static_cast<double>(current_ESS) / static_cast<double>(total_samples);
    unsigned int next_batch_size = 0;
    
    // The main idea is to ask for samples based on how many samples we need for 1 ESS.
    if (ess_per_sample > 1e-6) {
        unsigned int remaining_ESS = target_ESS - current_ESS;

        // we always ask for enough points to generate at least ~15 ESS.
        unsigned int min_viable_batch = static_cast<unsigned int>(15.0 / ess_per_sample);

        // Still, better keep an absolute minimum just in case ess_per_sample is oddly high
        min_viable_batch = std::max(min_viable_batch, 2500u);
        
        if (current_ESS < (target_ESS * 0.85)) {
            // We are less than 85% to the target.
            // Since ESS generation usually improves over time, we intentionally 
            // ask for only 80% of what the math predicts. This prevents massive overshoots.
            unsigned int samples_needed = static_cast<unsigned int>((remaining_ESS / ess_per_sample) * 0.80);
            next_batch_size = std::max(samples_needed, min_viable_batch);
        } else {
            // We are in the final 15%.
            // Now the ess_per_sample is highly accurate. We add a small 5% buffer 
            // to ensure we cross the finish line and don't waste time on tiny micro-batches.
            unsigned int samples_needed = static_cast<unsigned int>((remaining_ESS / ess_per_sample) * 1.05);
            next_batch_size = std::max(samples_needed, min_viable_batch);
        }
    } 
    else {
        // Fallback if ess_per_sample is virtually zero (terrible mixing)
        // back off the aggressive multiplier to prevent rapid spiraling
        next_batch_size = static_cast<unsigned int>(current_batch_size * 1.25);
    }
    
    // Fallback if ess_per_sample is virtually zero
    return std::min(next_batch_size, MAX_BATCH_SIZE);
}