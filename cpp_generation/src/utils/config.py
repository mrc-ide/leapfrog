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


def get_r_parse_data(cfg):
  r_alias = cfg["alias"]["r"]
  if cfg["type"] == "scalar":
    return f'Rcpp::as<{cfg["num_type"]}>(data["{r_alias}"])'
  else:
    return f'parse_data<{cfg["num_type"]}>(data, "{r_alias}", {collapse_dims(cfg)})'


def get_cpp_read_data(name, cfg):
  if cfg["type"] == "scalar":
    return f'read_data<{cfg["num_type"]}>(input_dir, "{name}")'
  else:
    return f'read_data<{cfg["num_type"]}>(input_dir, "{name}", {collapse_dims(cfg)})'


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


def dict_len(dict):
  return len(dict.keys())
