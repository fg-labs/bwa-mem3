// test/framework/scoped_rescue_env.h
//
// Scoped settings for the per-call kswv rescue toggles, and for the extension-dedup mode,
// so a unit test's kernel body does not depend on what the caller's environment set.

#ifndef BWA_TESTS_SCOPED_RESCUE_ENV_H
#define BWA_TESTS_SCOPED_RESCUE_ENV_H

#include <initializer_list>
#include <memory>
#include <string>
#include <vector>

#include "doctest/doctest.h"
#include "bwamem.h"       // mem_dedup_configure
#include "rescue_env.h"   // rescue_env_on
#include "scoped_env.h"

namespace bwa_tests {

// Every BWA3_RESCUE_* on/off toggle kswv reads on each call through rescue_env_on
// (src/kswv.cpp). Each defaults on; "0" turns it off.
constexpr const char *kKswvRescueToggles[] = {
    "BWA3_RESCUE_FSCAN",   "BWA3_RESCUE_USQADD",     "BWA3_RESCUE_ROWPAIR",     "BWA3_RESCUE_LAZYQE",
    "BWA3_RESCUE_PADCOLS", "BWA3_RESCUE_QE_REVSCAN", "BWA3_RESCUE_SCORE2_BITS",
};

// Sets every kswv toggle: the named ones to the given values, the rest to "1" (their
// default). Restores the caller's values on destruction. The constructor REQUIREs that
// rescue_env_on reads back exactly the scoped tuple.
class ScopedKswvRescueEnv {
public:
    struct Setting { const char *name, *value; };
    explicit ScopedKswvRescueEnv(std::initializer_list<Setting> overrides = {})
        : ScopedKswvRescueEnv(std::vector<Setting>(overrides)) {}
    explicit ScopedKswvRescueEnv(const std::vector<Setting> &overrides) {
        for (const char *name : kKswvRescueToggles) {
            const char *value = "1";
            for (const Setting &o : overrides) if (std::string(o.name) == name) value = o.value;
            envs_.emplace_back(new ScopedEnv(name, value));
            REQUIRE(rescue_env_on(name) == (value[0] != '0'));
        }
        for (const Setting &o : overrides) {
            bool known = false;
            for (const char *name : kKswvRescueToggles) known = known || std::string(o.name) == name;
            REQUIRE(known);   // a misspelt override would silently scope nothing
        }
    }
private:
    std::vector<std::unique_ptr<ScopedEnv>> envs_;
};

// Sets the extension-dedup mode for the scope (mem_dedup_configure: "off" | "on" | "auto")
// and on destruction resolves it again from BWAMEM3_DEDUP / the default, as an entry point
// that never configured it would.
// Sets the extension dedup mode for one scope. On exit it re-resolves the mode from the
// environment (BWAMEM3_DEDUP, or the default), which is the mode every test starts from;
// there is no getter, so a mode set by an enclosing scope is not restored.
class ScopedDedupMode {
public:
    explicit ScopedDedupMode(const char *mode) { mem_dedup_configure(mode); }
    ~ScopedDedupMode() { mem_dedup_configure(nullptr); }
    ScopedDedupMode(const ScopedDedupMode &) = delete;
    ScopedDedupMode &operator=(const ScopedDedupMode &) = delete;
};

} // namespace bwa_tests

#endif // BWA_TESTS_SCOPED_RESCUE_ENV_H
