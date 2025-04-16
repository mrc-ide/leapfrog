def get_member_type(val):
  if isinstance(val, int):
    return "int"
  elif isinstance(val[0], int):
    return f'std::array<int, {len(val)}>'
  else:
    return f'std::array<double, {len(val)}>'


def get_member_value(val):
  if isinstance(val, int):
    return val
  else:
    els = ", ".join(str(x) for x in val)
    return f'{{ {els} }}'
