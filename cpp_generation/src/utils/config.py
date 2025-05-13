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


def get_r_internal_data_pointer(cfg):
  return "INTEGER" if cfg["num_type"] == "int" else "REAL"


def get_r_parse_data(cfg, alias = None):
  r_alias = alias or cfg["alias"]["r"]
  if cfg["type"] == "scalar":
    return f'Rcpp::as<{cfg["num_type"]}>(data["{r_alias}"])'
  else:
    return f'parse_data<{cfg["num_type"]}>(data, "{r_alias}", {collapse_dims(cfg)})'


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
    return f'read_data<{cfg["num_type"]}>(params.{config_name}->{name}, params.{config_name}->{name}_length, "{name}", {collapse_dims(cfg)})'


def get_c_initial_state(config_name, name, cfg):
  config_name = config_name.lower()
  if cfg["type"] == "scalar":
    return f'*(state.{config_name}->{name})'
  else:
    return f'read_data<{cfg["num_type"]}>(state.{config_name}->{name}, state.{config_name}->{name}_length, "{name}", {collapse_dims(cfg)})'


def get_c_write_data(config_name, name, cfg):
  config_name = config_name.lower()
  return f'write_data<{cfg["num_type"]}, {dim_len(cfg) + 1}>(state.{name}, out.{config_name}->{name}, out.{config_name}->{name}_length, "{name}");'


def get_c_write_data_single_year(config_name, name, cfg):
  config_name = config_name.lower()
  if dim_len(cfg) == 0:
    return f'*(out.{config_name}->{name}) = state.{name};'
  else:
    return f'write_data<{cfg["num_type"]}, {dim_len(cfg)}>(state.{name}, out.{config_name}->{name}, out.{config_name}->{name}_length, "{name}");'


def get_delphi_setter(name, cfg):
  name = to_camel_case(name)
  return f'Set{ name }(var in{ name }: TGBFixedArray<{get_delphi_num_type(cfg)}>);'


def get_delphi_param_setter(config_name):
  return f'Set{ config_name }Params(var { config_name.lower() }Params: Leapfrog{ config_name }Params)'


def get_delphi_state_setter(config_name):
  return f'Set{ config_name }State(var { config_name.lower() }State: Leapfrog{ config_name }State)'


def get_reset(name, cfg):
  if cfg["type"] == "scalar":
    return f'{name} = {cfg.get("default") or 0}'
  else:
    return f'{name}.setZero()'


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


def get_num_type_delphi(cfg):
  n_dim = dim_len(cfg)
  if n_dim == 0:
    return cfg["num_type"]
  else:
    return f'{cfg["num_type"]}*'

def dict_len(dict):
  return len(dict.keys())


def to_camel_case(snake_str):
  return "".join(x.capitalize() for x in snake_str.lower().split("_"))


def to_lower_camel_case(snake_str):
  camel_string = to_camel_case(snake_str)
  return snake_str[0].lower() + camel_string[1:]


def get_delphi_ptr_type(cfg):
  return 'P' + get_delphi_num_type(cfg)


def get_delphi_num_type(cfg):
  if cfg["num_type"] == 'real_type':
    return 'Double'
  elif cfg["num_type"] == 'int':
    return 'Integer'
  else:
    raise Exception(f'Unknown type {cfg["num_type"]}, must be one of "real_type" or "int".')
