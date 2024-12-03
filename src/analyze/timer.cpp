#include <iostream>
#include <chrono>

class Timer {
public:
    Timer() : start_(std::chrono::high_resolution_clock::now()), elapsed_(0) {}

    // Start or restart the timer
    void start() {
        start_ = std::chrono::high_resolution_clock::now();
    }

    // Stop the timer and add to elapsed time
    void stop() {
        auto end = std::chrono::high_resolution_clock::now();
        elapsed_ += std::chrono::duration_cast<std::chrono::microseconds>(end - start_).count();
    }

    // Reset the timer
    void reset() {
        elapsed_ = 0;
    }

    // Get the total elapsed time in microseconds
    long long elapsedMicroseconds() const {
        return elapsed_;
    }

    // Get the total elapsed time in milliseconds
    double elapsedMilliseconds() const {
        return elapsed_ / 1000.0;
    }

    // Get the total elapsed time in seconds
    double elapsedSeconds() const {
        return elapsed_ / 1e6;
    }

private:
    std::chrono::time_point<std::chrono::high_resolution_clock> start_;
    long long elapsed_; // Total elapsed time in microseconds
};
