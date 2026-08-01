#pragma once

#include <algorithm>

// Checks if a number is "7-smooth" (only divisible by 2, 3, 5, or 7)
inline bool is_fft_friendly(unsigned int n) {
    if (n == 0) return false;
    while (n % 2 == 0) n /= 2;
    while (n % 3 == 0) n /= 3;
    while (n % 5 == 0) n /= 5;
    while (n % 7 == 0) n /= 7;
    return n == 1; // If it reduces exactly to 1, the FFT will process it instantly.
}

// Finds the nearest FFT-friendly number by bumping it up by 1 until it fits
inline unsigned int get_next_fft_friendly_size(unsigned int target) {
    unsigned int n = target;
    while (!is_fft_friendly(n)) {
        n++;
    }
    return n;
}

// Calculates the optimal next batch size based on the current ESS
inline unsigned int compute_next_batch_size(
    unsigned int target_ESS, 
    unsigned int current_ESS,
    unsigned int previous_ESS, 
    size_t total_samples, 
    unsigned int current_batch_size,
    unsigned int dimension,
    double remaining_time_sec,
    double samples_per_sec)
{
    unsigned int base_cap = std::max(500000u, target_ESS * 10); 
    unsigned int penalty = 1000 * dimension;
    unsigned int MAX_BATCH_SIZE = (base_cap > penalty) ? (base_cap - penalty) : 1000;
    MAX_BATCH_SIZE = std::max(1000u, std::min(MAX_BATCH_SIZE, 2000000u));

    double ess_per_sample = static_cast<double>(current_ESS) / static_cast<double>(total_samples);

    // The next lines attempt to catch a stuck sampler by checking the ESS efficiency between batches
    // calculate how efficiently this specific batch generated ESS
    unsigned int ess_gained = (current_ESS > previous_ESS) ? (current_ESS - previous_ESS) : 0;
    double marginal_ess_per_sample = static_cast<double>(ess_gained) / static_cast<double>(current_batch_size);
    size_t previous_samples = (total_samples > current_batch_size) ? (total_samples - current_batch_size) : 0;
    double previous_ess_per_sample = (previous_samples > 0) ? (static_cast<double>(previous_ESS) / static_cast<double>(previous_samples)) : 0.0;

    // warning check
    if (previous_samples > 0) {
        // Condition 1: efficiency dropped by 90% or more compared to historical average
        if (previous_ess_per_sample > 1e-6 && marginal_ess_per_sample < (previous_ess_per_sample * 0.1)) {
            std::cerr << " | [WARNING] Sampler might be stuck! ESS yield dropped heavily." << std::flush;
        } 
        // Condition 2: absolute terrible mixing (>10,000 samples per 1 ESS)
        else if (marginal_ess_per_sample < 1e-4) {
            std::cerr << " | [WARNING] Terrible mixing detected. Sampler may be stuck in a corner." << std::flush;
        }
    }

    unsigned int next_batch_size = 0;
    unsigned int min_viable_batch = 2500u; 

    // The main idea is to ask for samples based on how many samples we need for 1 ESS.
    if (ess_per_sample > 1e-6) {
        unsigned int remaining_ESS = target_ESS - current_ESS;

        // we always ask for enough points to generate at least ~15 ESS.
        min_viable_batch = static_cast<unsigned int>(15.0 / ess_per_sample);

        // Keep an absolute minimum, but never let the minimum exceed the calculated maximum
        min_viable_batch = std::max(min_viable_batch, 2500u);
        min_viable_batch = std::min(min_viable_batch, MAX_BATCH_SIZE); 
        
        if (current_ESS < (target_ESS * 0.85)) {
            // we are less than 85% to the target.
            unsigned int samples_needed = static_cast<unsigned int>((remaining_ESS / ess_per_sample) * 0.80);
            next_batch_size = std::max(samples_needed, min_viable_batch);
        } else {
            // we are in the final 15%.
            unsigned int samples_needed = static_cast<unsigned int>((remaining_ESS / ess_per_sample) * 1.10);
            next_batch_size = std::max(samples_needed, min_viable_batch);
        }
    } 
    else {
        // Fallback if ess_per_sample is virtually zero (terrible mixing)
        // Back off the aggressive multiplier. 
        next_batch_size = static_cast<unsigned int>(current_batch_size * 1.25);
        
        // If we are stuck, do not let it grow all the way to MAX_BATCH_SIZE.
        // Cap the runaway growth at 25% of our maximum (with a floor of 1000).
        unsigned int safe_stuck_cap = std::max(1000u, MAX_BATCH_SIZE / 4);
        next_batch_size = std::min(next_batch_size, safe_stuck_cap);
    }
    
    if (remaining_time_sec > 0 && samples_per_sec > 0) {
        unsigned int time_budget_batch = static_cast<unsigned int>(remaining_time_sec * samples_per_sec * 1.2);
        next_batch_size = std::min(next_batch_size, std::max(time_budget_batch, min_viable_batch));
    }

    // Calculate the raw batch size bounded by our maximums
    unsigned int raw_next_batch = std::min(next_batch_size, MAX_BATCH_SIZE);
    
    if (raw_next_batch == 0) return 0;

    // We calculate what the absolute total number of samples will be after this batch
    unsigned int target_total_samples = total_samples + raw_next_batch;
    
    // We bump the total up slightly (usually < 100 points) to the next Smooth Number
    unsigned int optimal_total_samples = get_next_fft_friendly_size(target_total_samples);
    
    // Return the adjusted batch size
    return optimal_total_samples - total_samples;
}