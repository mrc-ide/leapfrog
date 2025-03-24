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
