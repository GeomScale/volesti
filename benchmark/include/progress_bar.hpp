#pragma once

#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <algorithm>

#if defined(_WIN32)
#include <windows.h>
#else
#include <sys/ioctl.h>
#include <unistd.h>
#endif

inline unsigned int get_terminal_width(unsigned int fallback = 80) {
#if defined(_WIN32)
    CONSOLE_SCREEN_BUFFER_INFO csbi;
    if (GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &csbi)) {
        return static_cast<unsigned int>(csbi.srWindow.Right - csbi.srWindow.Left + 1);
    }
    return fallback;
#else
    struct winsize w;
    if (ioctl(STDOUT_FILENO, TIOCGWINSZ, &w) == 0 && w.ws_col > 0) {
        return w.ws_col;
    }
    return fallback;
#endif
}

inline void draw_progress_bar(const std::string& walk_name, 
    unsigned int current, 
    unsigned int total, 
    unsigned int total_samples_so_far,
    unsigned int current_ESS,
    unsigned int walk_len,
    double elapsed_seconds,
    double estimated_remaining_seconds,
    int bar_width = 25) 
{
    if (total == 0) return; 

    float progress = static_cast<float>(current) / total;
    int pos = static_cast<int>(bar_width * progress);

    std::ostringstream oss;

    oss << "[" << walk_name << "] Batch: [";
    for (int i = 0; i < bar_width; ++i) {
        if (i < pos) oss << "=";
        else if (i == pos) oss << ">";
        else oss << " ";
    }
    oss << "] " << static_cast<int>(progress * 100.0) << "% (" 
        << current << "/" << total << ")";

    if (current_ESS > 0) {
        double live_mixing_ratio = static_cast<double>(total_samples_so_far * walk_len) / current_ESS;
        oss << " | Mix Ratio: " << std::fixed << std::setprecision(2) << live_mixing_ratio;
    } else {
        oss << " | Mix Ratio: N/A";
    }

    unsigned int e_total_secs = static_cast<unsigned int>(elapsed_seconds);
    unsigned int e_hours = e_total_secs / 3600;
    unsigned int e_minutes = (e_total_secs % 3600) / 60;
    unsigned int e_seconds = e_total_secs % 60;

    oss << " | Elapsed: ";
    if (e_hours > 0) oss << e_hours << "h " << e_minutes << "m";
    else if (e_minutes > 0) oss << e_minutes << "m " << e_seconds << "s";
    else oss << e_seconds << "s";

    if (estimated_remaining_seconds >= 0.0) {
        unsigned int total_secs = static_cast<unsigned int>(estimated_remaining_seconds);
        unsigned int hours = total_secs / 3600;
        unsigned int minutes = (total_secs % 3600) / 60;
        unsigned int seconds = total_secs % 60;

        oss << " | ETA: ";
        if (hours > 0) oss << hours << "h " << minutes << "m";
        else if (minutes > 0) oss << minutes << "m " << seconds << "s";
        else oss << seconds << "s";
    } else {
        oss << " | ETA: Calculating...";
    }

    std::string line = oss.str();

    unsigned int term_width = get_terminal_width();
    unsigned int max_width = (term_width > 1) ? term_width - 1 : term_width;

    if (line.size() > max_width) {
        line = line.substr(0, max_width);
    }

    std::cout << "\r" << std::string(term_width, ' ') << "\r" << line << std::flush;
}