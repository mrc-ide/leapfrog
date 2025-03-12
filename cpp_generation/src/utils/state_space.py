def get_template_specialisation(name, ss):
  if ss.get("requires"):
    collapsed_requires = " && ".join([f'ModelVariant::{req}' for req in ss["requires"]])
    return (
      f'requires({collapsed_requires})\n'
      f'struct {name}SS<ModelVariant>'
    )
  else:
    return f'struct {name}SS'


def get_member_type(val):
  if isinstance(val, int):
    return "int"
  else:
    return f'std::array<int, {len(val)}>'


def get_member_value(val):
  if isinstance(val, int):
    return val
  else:
    els = ", ".join(str(x) for x in val)
    return f'{{ {els} }}'
