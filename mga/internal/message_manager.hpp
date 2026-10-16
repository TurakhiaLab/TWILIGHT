#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <queue>
#include <unordered_map>
#include <mutex>
#include <condition_variable>
#include <thread>
#include <atomic>
#include <sstream>
#include <memory>
#include <chrono>
#include <cstdarg>
#include <cstdio>
#include <ctime>
#include <iomanip>

namespace mga {

/**
 * @brief Output target stream mode.
 */
enum class OutputMode : int {
    Cout = 1,
    Cerr = 2
};

inline OutputMode toOutputMode(int mode) {
    return (mode == 2) ? OutputMode::Cerr : OutputMode::Cout;
}

/**
 * @brief Severity level for logging messages.
 */
enum class LogLevel : int {
    DEBUG = 0,
    INFO = 1,
    WARNING = 2,
    ERROR = 3,
    NONE = 4
};

namespace color {
    constexpr const char* RESET   = "\033[0m";
    constexpr const char* RED     = "\033[31m";
    constexpr const char* GREEN   = "\033[32m";
    constexpr const char* YELLOW  = "\033[33m";
    constexpr const char* BLUE    = "\033[34m";
    constexpr const char* MAGENTA = "\033[35m";
    constexpr const char* CYAN    = "\033[36m";
    constexpr const char* GRAY    = "\033[90m";
    constexpr const char* BOLD    = "\033[1m";
}

/**
 * @brief Strip ANSI escape codes from string for clean log file writing.
 */
inline std::string stripAnsiCodes(const std::string& input) {
    std::string result;
    result.reserve(input.size());
    bool inEscape = false;
    for (size_t i = 0; i < input.size(); ++i) {
        if (input[i] == '\033') {
            inEscape = true;
        } else if (inEscape) {
            if ((input[i] >= 'A' && input[i] <= 'Z') || (input[i] >= 'a' && input[i] <= 'z')) {
                inEscape = false;
            }
        } else {
            result += input[i];
        }
    }
    return result;
}

/**
 * @brief Thread-safe message controller and logger.
 * 
 * Features:
 * - Dedicated background worker thread with FIFO queue for serialized terminal output.
 * - Prevents interleaving of std::cout and std::cerr across multiple threads.
 * - Immediate output with print() and category buffering with push() / flush().
 * - Atomic batch flushing: flush(name) prints all buffered messages consecutively
 *   without being interleaved by real-time print() calls from other threads.
 * - Dual output (tee) to log file without ANSI color codes.
 * - Graceful shutdown and optional exit on ERROR.
 */
class MessageManager {
public:
    struct OutputBatch {
        std::vector<std::string> messages;
        OutputMode mode;
    };

    /**
     * @brief Meyers' Singleton accessor (thread-safe lazy initialization).
     */
    static MessageManager& getInstance() {
        static MessageManager instance;
        return instance;
    }

    // Disable copy and assignment
    MessageManager(const MessageManager&) = delete;
    MessageManager& operator=(const MessageManager&) = delete;

    ~MessageManager() {
        shutdown();
    }

    // ==========================================
    // 1. Log Level & Mode Control (Verbose / Quiet)
    // ==========================================
    void setLogLevel(LogLevel level) {
        minLogLevel_.store(level);
    }

    LogLevel getLogLevel() const {
        return minLogLevel_.load();
    }

    void setVerbose(bool enable = true) {
        verbose_.store(enable);
    }

    bool isVerbose() const {
        return verbose_.load();
    }

    void setQuiet(bool enable = true) {
        quiet_.store(enable);
    }

    bool isQuiet() const {
        return quiet_.load();
    }

    bool shouldLog(LogLevel level) const {
        // ERROR is always logged regardless of verbose/quiet settings
        if (level == LogLevel::ERROR) {
            return true;
        }
        // In quiet mode, only ERROR messages are output
        if (quiet_.load()) {
            return false;
        }
        // INFO is always logged unless quiet mode is active
        if (level == LogLevel::INFO) {
            return true;
        }
        // WARNING and DEBUG are only logged when verbose is enabled
        if (verbose_.load()) {
            return true;
        }
        return false;
    }

    // ==========================================
    // 2. Prefix Configuration (Timestamp, Thread ID, Level Tag, Color)
    // ==========================================
    void setEnableTimestamp(bool enable) { enableTimestamp_.store(enable); }
    void setEnableThreadId(bool enable)  { enableThreadId_.store(enable); }
    void setEnableLevelTag(bool enable)  { enableLevelTag_.store(enable); }
    void setEnableColor(bool enable)     { enableColor_.store(enable); }

    bool isTimestampEnabled() const { return enableTimestamp_.load(); }
    bool isThreadIdEnabled()  const { return enableThreadId_.load(); }
    bool isLevelTagEnabled()  const { return enableLevelTag_.load(); }
    bool isColorEnabled()     const { return enableColor_.load(); }

    // ==========================================
    // 3. File Logging (Tee mode)
    // ==========================================
    bool setLogFile(const std::string& filePath) {
        std::lock_guard<std::mutex> lock(fileMutex_);
        if (fileStream_.is_open()) {
            fileStream_.close();
        }
        fileStream_.open(filePath, std::ios::out | std::ios::app);
        if (!fileStream_.is_open()) {
            std::cerr << "[MessageManager] Failed to open log file: " << filePath << std::endl;
            return false;
        }
        logFilePath_ = filePath;
        return true;
    }

    void closeLogFile() {
        std::lock_guard<std::mutex> lock(fileMutex_);
        if (fileStream_.is_open()) {
            fileStream_.flush();
            fileStream_.close();
        }
        logFilePath_.clear();
    }

    // ==========================================
    // 6. Category Buffer Capacity Protection
    // ==========================================
    void setCategoryLimit(size_t limit) {
        categoryLimit_.store(limit);
    }

    size_t getCategoryLimit() const {
        return categoryLimit_.load();
    }

    // ==========================================
    // 8. Error Handling (Exit on Error)
    // ==========================================
    void setExitOnError(bool enable) {
        exitOnError_.store(enable);
    }

    bool isExitOnError() const {
        return exitOnError_.load();
    }

    // ==========================================
    // Immediate Output
    // ==========================================
    void print(const std::string& msg, OutputMode mode = OutputMode::Cout) {
        {
            std::lock_guard<std::mutex> lock(mainMutex_);
            mainQueue_.push(OutputBatch{{msg}, mode});
        }
        cv_.notify_one();
    }

    void print(const std::string& msg, int mode) {
        print(msg, toOutputMode(mode));
    }

    // ==========================================
    // Category Buffering (push & flush)
    // ==========================================
    void push(const std::string& name, const std::string& msg) {
        bool needAutoFlush = false;
        {
            std::lock_guard<std::mutex> lock(categoryMutex_);
            auto& buf = categoryBuffers_[name];
            buf.push_back(msg);

            // Check if buffer capacity exceeds limit to prevent unbounded memory growth
            size_t limit = categoryLimit_.load();
            if (limit > 0 && buf.size() >= limit) {
                needAutoFlush = true;
            }
        }

        if (needAutoFlush) {
            std::string warnMsg = formatPrefix(LogLevel::WARNING) +
                                  "Warning: Category buffer '" + name +
                                  "' reached capacity limit (" + std::to_string(categoryLimit_.load()) +
                                  "), auto-flushing to terminal.\n";
            print(warnMsg, OutputMode::Cerr);
            flush(name, OutputMode::Cout);
        }
    }

    void flush(const std::string& name, OutputMode mode = OutputMode::Cout) {
        std::vector<std::string> toFlush;
        {
            std::lock_guard<std::mutex> lock(categoryMutex_);
            auto it = categoryBuffers_.find(name);
            if (it != categoryBuffers_.end() && !it->second.empty()) {
                toFlush = std::move(it->second);
                categoryBuffers_.erase(it);
            }
        }

        if (!toFlush.empty()) {
            {
                std::lock_guard<std::mutex> lock(mainMutex_);
                mainQueue_.push(OutputBatch{std::move(toFlush), mode});
            }
            cv_.notify_one();
        }
    }

    void flush(const std::string& name, int mode) {
        flush(name, toOutputMode(mode));
    }

    void clear(const std::string& name) {
        std::lock_guard<std::mutex> lock(categoryMutex_);
        categoryBuffers_.erase(name);
    }

    void clearAllCategories() {
        std::lock_guard<std::mutex> lock(categoryMutex_);
        categoryBuffers_.clear();
    }

    /**
     * @brief Synchronously wait until all queued messages have been printed.
     */
    void waitEmpty() {
        std::unique_lock<std::mutex> lock(mainMutex_);
        emptyCv_.wait(lock, [this]() {
            return mainQueue_.empty() && !isProcessing_;
        });
    }

    // ==========================================
    // 7. Formatted Output Core Logic
    // ==========================================
    std::string formatPrefix(LogLevel level) const {
        std::ostringstream ss;
        bool color = enableColor_.load();

        // 1. Timestamp [YYYY-MM-DD HH:MM:SS.mmm]
        if (enableTimestamp_.load()) {
            auto now = std::chrono::system_clock::now();
            auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;
            std::time_t timer = std::chrono::system_clock::to_time_t(now);
            std::tm bt{};
            localtime_r(&timer, &bt);

            if (color) ss << color::GRAY;
            char timeBuf[32];
            std::strftime(timeBuf, sizeof(timeBuf), "[%Y-%m-%d %H:%M:%S", &bt);
            ss << timeBuf << "." << std::setfill('0') << std::setw(3) << ms.count() << "] ";
            if (color) ss << color::RESET;
        }

        // 2. Thread ID [TID: ...]
        if (enableThreadId_.load()) {
            if (color) ss << color::CYAN;
            ss << "[TID:" << std::this_thread::get_id() << "] ";
            if (color) ss << color::RESET;
        }

        // 3. Log Level Tag [INFO], [WARN], etc. (with color)
        if (enableLevelTag_.load()) {
            switch (level) {
                case LogLevel::DEBUG:
                    if (color) ss << color::GRAY << "[DEBUG] " << color::RESET;
                    else ss << "[DEBUG] ";
                    break;
                case LogLevel::INFO:
                    if (color) ss << color::GREEN << "[INFO]  " << color::RESET;
                    else ss << "[INFO]  ";
                    break;
                case LogLevel::WARNING:
                    if (color) ss << color::YELLOW << "[WARN]  " << color::RESET;
                    else ss << "[WARN]  ";
                    break;
                case LogLevel::ERROR:
                    if (color) ss << color::RED << "[ERROR] " << color::RESET;
                    else ss << "[ERROR] ";
                    break;
                default:
                    break;
            }
        }

        return ss.str();
    }

    void logDirect(LogLevel level, const std::string& msg, OutputMode mode) {
        if (!shouldLog(level)) return;
        std::string fullMsg = formatPrefix(level) + msg;
        print(fullMsg, mode);

        // When logging an ERROR and exitOnError is enabled, flush all and terminate process
        if (level == LogLevel::ERROR && exitOnError_.load()) {
            waitEmpty();
            std::exit(1);
        }
    }

    void logPush(const std::string& name, LogLevel level, const std::string& msg) {
        if (!shouldLog(level)) return;
        std::string fullMsg = formatPrefix(level) + msg;
        push(name, fullMsg);
    }

private:
    MessageManager()
        : minLogLevel_(LogLevel::INFO),
          enableTimestamp_(false),
          enableThreadId_(false),
          enableLevelTag_(true),
          enableColor_(true),
          categoryLimit_(10000),
          exitOnError_(true),
          verbose_(false),
          quiet_(false),
          stopFlag_(false),
          isProcessing_(false) {
        workerThread_ = std::thread(&MessageManager::workerLoop, this);
    }

    void shutdown() {
        {
            std::lock_guard<std::mutex> lock(mainMutex_);
            if (stopFlag_) return;
            stopFlag_ = true;
        }
        cv_.notify_all();

        if (workerThread_.joinable()) {
            workerThread_.join();
        }

        closeLogFile();
    }

    void workerLoop() {
        while (true) {
            OutputBatch batch;
            {
                std::unique_lock<std::mutex> lock(mainMutex_);
                cv_.wait(lock, [this]() {
                    return stopFlag_ || !mainQueue_.empty();
                });

                if (stopFlag_ && mainQueue_.empty()) {
                    break;
                }

                if (!mainQueue_.empty()) {
                    isProcessing_ = true;
                    batch = std::move(mainQueue_.front());
                    mainQueue_.pop();
                } else {
                    continue;
                }
            }

            // If log file is open, output to log file instead of terminal
            bool wroteToLog = false;
            {
                std::lock_guard<std::mutex> lock(fileMutex_);
                if (fileStream_.is_open()) {
                    for (const auto& line : batch.messages) {
                        fileStream_ << stripAnsiCodes(line);
                    }
                    fileStream_.flush();
                    wroteToLog = true;
                }
            }

            // Output to terminal only if no log file is active
            if (!wroteToLog) {
                std::ostream& os = (batch.mode == OutputMode::Cerr) ? std::cerr : std::cout;
                for (const auto& line : batch.messages) {
                    os << line;
                }
                os.flush();
            }

            {
                std::lock_guard<std::mutex> lock(mainMutex_);
                isProcessing_ = false;
                if (mainQueue_.empty()) {
                    emptyCv_.notify_all();
                }
            }
        }
    }

    // Configuration flags
    std::atomic<LogLevel> minLogLevel_;
    std::atomic<bool> enableTimestamp_;
    std::atomic<bool> enableThreadId_;
    std::atomic<bool> enableLevelTag_;
    std::atomic<bool> enableColor_;
    std::atomic<size_t> categoryLimit_;
    std::atomic<bool> exitOnError_;
    std::atomic<bool> verbose_;
    std::atomic<bool> quiet_;

    // File logging members
    std::string logFilePath_;
    std::ofstream fileStream_;
    std::mutex fileMutex_;

    // Main FIFO queue members
    std::queue<OutputBatch> mainQueue_;
    std::mutex mainMutex_;
    std::condition_variable cv_;
    std::condition_variable emptyCv_;
    std::atomic<bool> stopFlag_;
    std::atomic<bool> isProcessing_;
    std::thread workerThread_;

    // Category buffer members
    std::unordered_map<std::string, std::vector<std::string>> categoryBuffers_;
    std::mutex categoryMutex_;
};

// ==========================================
// Formatting Helper Function
// ==========================================
inline std::string vformatHelper(const char* fmt, va_list args) {
    va_list argsCopy;
    va_copy(argsCopy, args);
    int size = std::vsnprintf(nullptr, 0, fmt, argsCopy);
    va_end(argsCopy);

    if (size <= 0) return std::string();

    std::string buf(size, '\0');
    std::vsnprintf(&buf[0], size + 1, fmt, args);
    return buf;
}

// ==========================================
// Stream Builder Class (operator<< syntax)
// ==========================================
class MessageStream {
public:
    enum class Target { DirectPrint, PushCategory };

    MessageStream(LogLevel level, OutputMode mode = OutputMode::Cout)
        : target_(Target::DirectPrint), level_(level), mode_(mode) {}

    MessageStream(const std::string& category, LogLevel level)
        : target_(Target::PushCategory), categoryName_(category), level_(level), mode_(OutputMode::Cout) {}

    template <typename T>
    MessageStream& operator<<(const T& value) {
        oss_ << value;
        return *this;
    }

    ~MessageStream() {
        std::string content = oss_.str();
        if (content.empty()) return;

        if (target_ == Target::DirectPrint) {
            MessageManager::getInstance().logDirect(level_, content, mode_);
        } else {
            MessageManager::getInstance().logPush(categoryName_, level_, content);
        }
    }

private:
    Target target_;
    std::string categoryName_;
    LogLevel level_;
    OutputMode mode_;
    std::ostringstream oss_;
};

/**
 * @brief RAII Scoped category flush helper.
 * Automatically flushes the category when leaving the scope.
 */
class ScopedCategoryFlush {
public:
    ScopedCategoryFlush(std::string name, OutputMode mode = OutputMode::Cout)
        : name_(std::move(name)), mode_(mode), active_(true) {}

    ScopedCategoryFlush(std::string name, int mode)
        : name_(std::move(name)), mode_(toOutputMode(mode)), active_(true) {}

    ~ScopedCategoryFlush() {
        if (active_) {
            MessageManager::getInstance().flush(name_, mode_);
        }
    }

    void cancel() { active_ = false; }

private:
    std::string name_;
    OutputMode mode_;
    bool active_;
};

// ==========================================
// Convenience Global Wrappers (Namespace: mga::msg)
// ==========================================
namespace msg {

inline void setLogLevel(LogLevel level) { MessageManager::getInstance().setLogLevel(level); }
inline void setVerbose(bool enable = true) { MessageManager::getInstance().setVerbose(enable); }
inline bool isVerbose() { return MessageManager::getInstance().isVerbose(); }
inline void setQuiet(bool enable = true) { MessageManager::getInstance().setQuiet(enable); }
inline bool isQuiet() { return MessageManager::getInstance().isQuiet(); }
inline void init(bool verbose, bool quiet, const std::string& logFile = "") {
    MessageManager::getInstance().setVerbose(verbose);
    MessageManager::getInstance().setQuiet(quiet);
    if (!logFile.empty()) {
        MessageManager::getInstance().setLogFile(logFile);
    }
}
inline void setLogFile(const std::string& path) { MessageManager::getInstance().setLogFile(path); }
inline void enableTimestamp(bool enable = true) { MessageManager::getInstance().setEnableTimestamp(enable); }
inline void enableThreadId(bool enable = true)  { MessageManager::getInstance().setEnableThreadId(enable); }
inline void enableLevelTag(bool enable = true)  { MessageManager::getInstance().setEnableLevelTag(enable); }
inline void enableColor(bool enable = true)     { MessageManager::getInstance().setEnableColor(enable); }
inline void setCategoryLimit(size_t limit)      { MessageManager::getInstance().setCategoryLimit(limit); }
inline void setExitOnError(bool enable = true)  { MessageManager::getInstance().setExitOnError(enable); }

inline void print(const std::string& msg, OutputMode mode = OutputMode::Cout) {
    MessageManager::getInstance().print(msg, mode);
}

inline void print(const std::string& msg, int mode) {
    MessageManager::getInstance().print(msg, mode);
}

inline void push(const std::string& name, const std::string& msg) {
    MessageManager::getInstance().push(name, msg);
}

inline void flush(const std::string& name, OutputMode mode = OutputMode::Cout) {
    MessageManager::getInstance().flush(name, mode);
}

inline void flush(const std::string& name, int mode) {
    MessageManager::getInstance().flush(name, mode);
}

inline void clear(const std::string& name) {
    MessageManager::getInstance().clear(name);
}

inline void waitEmpty() {
    MessageManager::getInstance().waitEmpty();
}

// Formatted logging APIs (immediate print with prefix and color)
inline void logf(LogLevel level, OutputMode mode, const char* fmt, ...) {
    if (!MessageManager::getInstance().shouldLog(level)) return;
    va_list args;
    va_start(args, fmt);
    std::string msg = vformatHelper(fmt, args);
    va_end(args);
    MessageManager::getInstance().logDirect(level, msg, mode);
}

inline void debug(const char* fmt, ...) {
    if (!MessageManager::getInstance().shouldLog(LogLevel::DEBUG)) return;
    va_list args;
    va_start(args, fmt);
    std::string msg = vformatHelper(fmt, args);
    va_end(args);
    MessageManager::getInstance().logDirect(LogLevel::DEBUG, msg, OutputMode::Cout);
}

inline void info(const char* fmt, ...) {
    if (!MessageManager::getInstance().shouldLog(LogLevel::INFO)) return;
    va_list args;
    va_start(args, fmt);
    std::string msg = vformatHelper(fmt, args);
    va_end(args);
    MessageManager::getInstance().logDirect(LogLevel::INFO, msg, OutputMode::Cout);
}

inline void warn(const char* fmt, ...) {
    if (!MessageManager::getInstance().shouldLog(LogLevel::WARNING)) return;
    va_list args;
    va_start(args, fmt);
    std::string msg = vformatHelper(fmt, args);
    va_end(args);
    MessageManager::getInstance().logDirect(LogLevel::WARNING, msg, OutputMode::Cerr);
}

inline void error(const char* fmt, ...) {
    if (!MessageManager::getInstance().shouldLog(LogLevel::ERROR)) return;
    va_list args;
    va_start(args, fmt);
    std::string msg = vformatHelper(fmt, args);
    va_end(args);
    MessageManager::getInstance().logDirect(LogLevel::ERROR, msg, OutputMode::Cerr);
}

// Formatted logging APIs (buffered into category)
inline void pushf(const std::string& name, LogLevel level, const char* fmt, ...) {
    if (!MessageManager::getInstance().shouldLog(level)) return;
    va_list args;
    va_start(args, fmt);
    std::string msg = vformatHelper(fmt, args);
    va_end(args);
    MessageManager::getInstance().logPush(name, level, msg);
}

// Stream builders
inline MessageStream stream(LogLevel level = LogLevel::INFO, OutputMode mode = OutputMode::Cout) {
    return MessageStream(level, mode);
}

inline MessageStream categoryStream(const std::string& name, LogLevel level = LogLevel::INFO) {
    return MessageStream(name, level);
}

} // namespace msg
} // namespace mga
