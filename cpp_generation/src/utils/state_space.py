def get_member_type(val):
  if isinstance(val, int):
    return "int"
  elif isinstance(val[0], int):
    return f'std::array<int, {len(val)}>'
  elif isinstance(val[0], float):
    return f'std::array<double, {len(val)}>'
  else:
    return f'std::array<{get_member_type(val[0])}, {len(val)}>'


def get_member_value(val):
  if isinstance(val, (int, float)):
    return val
  if isinstance(val[0], (int, float)):
    els = ", ".join(str(x) for x in val)
    return f'{{ {els} }}'
  else:
    els = ", ".join(str(get_member_value(x)) for x in val)
    return f'{{{{ {els} }}}}'


def extract_state_space_info(dat):
  base = dat["global_state_space"]
  configs_and_overrides = []
  for config in dat["configs"]:
    ss = config["state_space"]
    default = ss["default"]
    config_or_override = {
      "condition": f'ModelVariant::{config["enable_if"]}',
      "name": config["name"],
      "vars": default
    }
    configs_and_overrides.append(config_or_override)

    overrides = ss.get("overrides")
    if not overrides: continue
    for index, override in enumerate(overrides):
      config_or_override = {
        "condition": f'ModelVariant::{config["enable_if"]} && ModelVariant::{override["applies_if"]}',
        "name": f'{config["name"]}Override{index}',
        "vars": override["vars"]
      }
      configs_and_overrides.append(config_or_override)

  # we need to reverse it because the template iteration with inheritance gives
  # higher overload priority to the start structs rather than the end ones
  configs_and_overrides.reverse()

  return { "base": base, "configs_and_overrides": configs_and_overrides }
