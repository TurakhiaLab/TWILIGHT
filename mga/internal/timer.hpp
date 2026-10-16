#pragma once

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <tbb/spin_rw_mutex.h>

class Timer {
private:
    struct TimerSlot {
        std::chrono::steady_clock::time_point startTime;
        std::chrono::nanoseconds accumulatedTime{0};
        bool active = false;
    };

    std::unordered_map<std::string, TimerSlot> slots_;
    bool isGlobalActive = true; 
    mutable tbb::spin_rw_mutex mutex_;

    std::pair<double, std::string> convertDuration(std::chrono::nanoseconds totalNs, const std::string& unit) const {
        if (unit == "s")  return { totalNs.count() / 1e9, "s" };
        if (unit == "us") return { totalNs.count() / 1e3, "us" };
        if (unit == "ns") return { static_cast<double>(totalNs.count()), "ns" };
        return { totalNs.count() / 1e6, "ms" };
    }

    void printSlot(const std::string& name, const TimerSlot& slot, const std::string& unit) const {
        auto [value, unitStr] = convertDuration(slot.accumulatedTime, unit);
        std::cout << "  ├─ " << std::left << std::setw(25) << name << " : " 
                  << std::right << std::setw(10) << std::fixed << std::setprecision(3) << value 
                  << " " << unitStr;
        if (slot.active) {
            std::cout << " [RUNNING...]";
        }
        std::cout << "\n";
    }

public:
    Timer() = default;

    void init() {
        clear();
    }

    void activate() {
        tbb::spin_rw_mutex::scoped_lock lock(mutex_, /*write=*/true);
        isGlobalActive = true;
    }

    void deactivate() {
        tbb::spin_rw_mutex::scoped_lock lock(mutex_, /*write=*/true);
        isGlobalActive = false;
    }

    void start(const std::string& name) {
        tbb::spin_rw_mutex::scoped_lock lock(mutex_, /*write=*/true);
        if (!isGlobalActive) return;

        auto& slot = slots_[name]; 
        slot.startTime = std::chrono::steady_clock::now();
        slot.active = true;
    }

    void stop(const std::string& name) {
        auto now = std::chrono::steady_clock::now();
        tbb::spin_rw_mutex::scoped_lock lock(mutex_, /*write=*/true);
        if (!isGlobalActive) return;

        auto it = slots_.find(name);
        if (it != slots_.end() && it->second.active) {
            auto& slot = it->second;
            auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(now - slot.startTime);
            slot.accumulatedTime += duration;
            slot.active = false;
        }
    }

    double getTime(const std::string& name, const std::string& unit = "ms") const {
        tbb::spin_rw_mutex::scoped_lock lock(mutex_, /*write=*/false);
        auto it = slots_.find(name);
        if (it != slots_.end()) {
            return convertDuration(it->second.accumulatedTime, unit).first;
        }
        return 0.0;
    }

    // Add accumulated time directly with nanoseconds duration (thread-safe)
    void add(const std::string& name, std::chrono::nanoseconds duration) {
        tbb::spin_rw_mutex::scoped_lock lock(mutex_, /*write=*/true);
        if (!isGlobalActive) return;
        slots_[name].accumulatedTime += duration;
    }

    // Add accumulated time directly with milliseconds duration (thread-safe)
    void add(const std::string& name, double duration_ms) {
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::duration<double, std::milli>(duration_ms));
        add(name, duration);
    }

    void addTime(const std::string& name, std::chrono::nanoseconds duration) {
        add(name, duration);
    }

    void addTime(const std::string& name, double duration_ms) {
        add(name, duration_ms);
    }

    // Merge another timer's accumulated records into this timer (thread-safe)
    void merge(const Timer& other) {
        tbb::spin_rw_mutex::scoped_lock lock_other(other.mutex_, /*write=*/false);
        tbb::spin_rw_mutex::scoped_lock lock_self(mutex_, /*write=*/true);
        if (!isGlobalActive) return;
        for (const auto& [name, slot] : other.slots_) {
            slots_[name].accumulatedTime += slot.accumulatedTime;
        }
    }

    void print(const std::string& name = "", const std::string& unit = "ms") const {
        tbb::spin_rw_mutex::scoped_lock lock(mutex_, /*write=*/false);
        std::cout << "\n==================== PERFORMANCE REPORT ====================\n";
        if (!name.empty()) {
            auto it = slots_.find(name);
            if (it != slots_.end()) {
                printSlot(name, it->second, unit);
            } else {
                std::cout << "  ⚠️ Timer '" << name << "' does not exist.\n";
            }
        } else {
            if (slots_.empty()) {
                std::cout << "  (No timers recorded)\n";
            } else {
                std::vector<std::pair<std::string, TimerSlot>> sorted_slots(slots_.begin(), slots_.end());
                std::sort(sorted_slots.begin(), sorted_slots.end(),
                          [](const auto& a, const auto& b) {
                              return a.first < b.first;
                          });
                for (const auto& [timerName, slot] : sorted_slots) {
                    printSlot(timerName, slot, unit);
                }
            }
        }
        std::cout << "============================================================\n\n";
    }

    void clear() {
        tbb::spin_rw_mutex::scoped_lock lock(mutex_, /*write=*/true);
        slots_.clear();
    }
};

inline Timer global_timer;