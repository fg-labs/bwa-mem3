// test/framework/scoped_env.h
//
// RAII override of an environment variable for unit tests that flip the
// BWA3_RESCUE_* kernel toggles in-process.

#ifndef BWA_TESTS_SCOPED_ENV_H
#define BWA_TESTS_SCOPED_ENV_H

#include <cstdlib>
#include <string>

namespace bwa_tests {

// Set on construction, restore the prior value (or unset if it was absent) on
// destruction, so doctest run order cannot leak an override into a sibling
// case that assumes the defaults.
class ScopedEnv {
public:
    ScopedEnv(const char *name, const char *value) : name_(name) {
        const char *cur = getenv(name);
        had_ = cur != nullptr;
        if (had_) saved_ = cur;
        setenv(name, value, 1);
    }
    ~ScopedEnv() {
        if (had_) setenv(name_.c_str(), saved_.c_str(), 1);
        else      unsetenv(name_.c_str());
    }
    ScopedEnv(const ScopedEnv &) = delete;
    ScopedEnv &operator=(const ScopedEnv &) = delete;
private:
    std::string name_;
    std::string saved_;
    bool had_ = false;
};

} // namespace bwa_tests

#endif // BWA_TESTS_SCOPED_ENV_H
