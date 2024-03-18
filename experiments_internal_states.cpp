#include "experiments_internal_states.h"
#include <rfl.hpp>
#include <rfl/json.hpp>

namespace exp_states {
#define DEFINE_EXP_STATE_TO_JSON_FUNCTION(Type)    \
  auto to_json(const Type& obj) noexcept -> json { \
    /* reflect-cpp -> yyjson -> nlohmann::json */  \
    auto json_str = rfl::json::write(obj);         \
    return json::parse(json_str);                  \
  }

EXP_STATES_TO_JSON_FUNCTION_AUTO_GENERATED(DEFINE_EXP_STATE_TO_JSON_FUNCTION)
} // namespace exp_states

#undef DEFINE_EXP_STATE_TO_JSON_FUNCTION
