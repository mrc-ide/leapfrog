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
