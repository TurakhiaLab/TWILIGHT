#pragma once

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <chrono>
#include <iomanip>

class Timer {
private:
    struct TimerSlot {
        std::chrono::steady_clock::time_point startTime;
        std::chrono::nanoseconds accumulatedTime{0};
        bool active = false;
    };

    std::unordered_map<std::string, TimerSlot> slots_;
    bool isGlobalActive = true; 

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
        isGlobalActive = true;
    }

    void deactivate() {
        isGlobalActive = false;
    }

    void start(const std::string& name) {
        if (!isGlobalActive) return;

        auto& slot = slots_[name]; 
        slot.startTime = std::chrono::steady_clock::now();
        slot.active = true;
    }

    void stop(const std::string& name) {
        auto now = std::chrono::steady_clock::now();
        if (!isGlobalActive) return;

        auto it = slots_.find(name);
        if (it != slots_.end() && it->second.active) {
            auto& slot = it->second;
            auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(now - slot.startTime);
            slot.accumulatedTime += duration;
            slot.active = false;
        }
    }

    void print(const std::string& name = "", const std::string& unit = "ms") const {
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
        slots_.clear();
    }
};

inline Timer global_timer;