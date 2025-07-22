def collapse_dims(cfg, sep = ", "):
  if cfg.get("dims"):
    return sep.join(cfg["dims"])
  else:
    return ""


def collapse_dims_with_trailing_sep(cfg, sep = ", "):
  if cfg.get("dims"):
    dims_str = sep.join(cfg["dims"])
    return dims_str + sep
  else:
    return ""


def collapse(l, sep = ", "):
  return sep.join(l)


def nda_all_chip(cfg):
  return ", ".join(["nda::_"] * (len(cfg["dims"]) - 1))


def dim_len(cfg):
  if not cfg.get("dims"):
    return 0
  return len(cfg["dims"])


def get_member_type(cfg):
  if cfg["type"] == "TM":
    return f'TM{dim_len(cfg)}<{cfg["num_type"]}>'
  elif cfg['type'] == "TFS":
    return f'TFS<{cfg["num_type"]}, {collapse_dims(cfg)}>'
  else:
    return f'{cfg["num_type"]}'


def get_shape_dim_types(cfg):
  all_static_dims = True
  shape_dim_types = []
  for index, dim in enumerate(cfg["dims"]):
    if not all_static_dims:
      shape_dim_types.append("nda::dim<0, nda::dynamic, nda::dynamic>")
    elif "opts." in dim or "output_years" in dim:
      stride = "1" if index == 0  else " * ".join([f"({d})" for d in cfg["dims"][:index]])
      shape_dim_types.append(f"nda::dim<0, nda::dynamic, {stride}>")
      all_static_dims = False
    else:
      stride = "1" if index == 0 else " * ".join([f"({d})" for d in cfg["dims"][:index]])
      shape_dim_types.append(f"nda::dim<0, {dim}, {stride}>")
  return shape_dim_types


def get_py_parse_data(name, cfg):
  # TODO: MANTRA remove r alias everywhere
  r_alias = cfg["alias"]["r"]
  if cfg["type"] == "scalar":
    return f'nb::cast<{cfg["num_type"]}>(data["{r_alias}"])'
  else:
    dim_types = []
    for index, dim in enumerate(cfg["dims"]):
      prod_prev_dim = "1" if index == 0 else " * ".join([f"({d})" for d in cfg["dims"][:index]])
      dim_types.append(f'nda::dim<>(0, {dim}, {prod_prev_dim})')
    return f'parse_data<{cfg["num_type"]}, {len(cfg["dims"])}>(data, "{r_alias}", {{ {", ".join(dim_types)} }})'


def get_py_initial_state(cfg, name):
  if cfg["type"] == "scalar":
      return f'state.{name} = nb::cast<{cfg["num_type"]}>(data["{name}"])'
  else:
    shape_path = f'Config::State::shape_{name}'
    return f'fill_initial_state<{cfg["num_type"]}, typename {shape_path}>(data, "{name}", state.{name})'


def get_r_internal_data_pointer(cfg):
  return "INTEGER" if cfg["num_type"] == "int" else "REAL"


def get_r_parse_data(cfg):
  r_alias = cfg["alias"]["r"]
  if cfg["type"] == "scalar":
    return f'Rcpp::as<{cfg["num_type"]}>(data["{r_alias}"])'
  else:
    dim_types = []
    for index, dim in enumerate(cfg["dims"]):
      prod_prev_dim = "1" if index == 0 else " * ".join([f"({d})" for d in cfg["dims"][:index]])
      dim_types.append(f'nda::dim<>(0, {dim}, {prod_prev_dim})')
    return f'parse_data<{cfg["num_type"]}, {len(cfg["dims"])}>(data, "{r_alias}", {{ {", ".join(dim_types)} }})'


def get_r_initial_state(cfg, name):
  if cfg["type"] == "scalar":
    return f'state.{name} = Rcpp::as<{cfg["num_type"]}>(data["{name}"])'
  else:
    shape_path = f'Config::State::shape_{name}'
    return f'fill_initial_state<{cfg["num_type"]}, typename {shape_path}>(data, "{name}", state.{name})'


def get_parse_pars_cpp(name, cfg):
  r_alias = cfg["alias"]["r"]
  if cfg["type"] == "scalar":
    return f'read_data_scalar<{cfg["num_type"]}>(params_file, "{r_alias}")'
  else:
    dim_types = []
    for index, dim in enumerate(cfg["dims"]):
      prod_prev_dim = "1" if index == 0 else " * ".join([f"({d})" for d in cfg["dims"][:index]])
      dim_types.append(f'nda::dim<>(0, {dim}, {prod_prev_dim})')
    shape = f'typename Pars::shape_{name}'
    return f'read_data<{cfg["num_type"]}, {shape}>(params_file, "{r_alias}", {{ {", ".join(dim_types)} }})'


def get_pars_cpp(name, cfg, ns):
  el = f'owned_pars.{ns}.{name}'
  if cfg["type"] == "scalar":
    return el
  else:
    return f'{{ {el}.data(), {el}.shape() }}'


def get_cpp_read_data(name, cfg):
  if cfg["type"] == "scalar":
    return f'read_data<{cfg["num_type"]}>(input_dir, "{name}")'
  else:
    return f'read_data<{cfg["num_type"]}>(input_dir, "{name}", {collapse_dims(cfg)})'


def get_c_read_data(config_name, name, cfg):
  config_name = config_name.lower()
  if cfg["type"] == "scalar":
    return f'params.{config_name}->{name}'
  else:
    dim_types = []
    for index, dim in enumerate(cfg["dims"]):
      prod_prev_dim = "1" if index == 0 else " * ".join([f"({d})" for d in cfg["dims"][:index]])
      dim_types.append(f'nda::dim<>(0, {dim}, {prod_prev_dim})')
    return f'read_data<{cfg["num_type"]}, {dim_len(cfg)}>(params.{config_name}->{name}, params.{config_name}->{name}_length, "{name}", {{ {", ".join(dim_types)} }})'


def get_c_initial_state(config_name, name, cfg):
  config_name = config_name.lower()
  if cfg["type"] == "scalar":
    return f'initial_state.{name} = *(state.{config_name}->{name})'
  else:
    shape_path = f'Config::State::shape_{name}'
    return f'fill_initial_state<{cfg["num_type"]}, typename {shape_path}>(state.{config_name}->{name}, state.{config_name}->{name}_length, "{name}", initial_state.{name})'


def get_c_write_data(config_name, name, cfg):
  config_name = config_name.lower()
  shape_path = f'Config::OutputState::shape_{name}'
  return f'write_data<{cfg["num_type"]}, typename {shape_path}>(state.{name}, out.{config_name}->{name}, out.{config_name}->{name}_length, "{name}");'


def get_c_write_data_single_year(config_name, name, cfg):
  config_name = config_name.lower()
  if dim_len(cfg) == 0:
    return f'*(out.{config_name}->{name}) = state.{name};'
  else:
    shape_path = f'Config::State::shape_{name}'
    return f'write_data<{cfg["num_type"]}, typename {shape_path}>(state.{name}, out.{config_name}->{name}, out.{config_name}->{name}_length, "{name}");'


def get_reset_value(cfg):
  return cfg.get("default") or 0


def get_output_state_member_type(cfg):
  n_dim = dim_len(cfg) + 1 if cfg.get("dims") else 1
  return f'T{n_dim}<{cfg["num_type"]}>'


def get_output_state_chip(name, cfg):
  if cfg["type"] == "scalar":
    return f'{name}(i)'
  else:
    return f'{name}.chip(i, {name}.NumDimensions - 1)'


def get_num_type(cfg):
  return cfg["num_type"]
