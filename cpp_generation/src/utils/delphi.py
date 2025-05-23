from utils.config import dim_len
from utils.general import to_camel_case

def get_num_type_c(cfg):
  n_dim = dim_len(cfg)
  if n_dim == 0:
    return cfg["num_type"]
  else:
    return f'{cfg["num_type"]}*'


def get_delphi_param_type(cfg):
  n_dim = dim_len(cfg)
  if n_dim == 0:
    return get_delphi_num_type(cfg)
  else:
    return f'TGBFixedArray<{ get_delphi_num_type(cfg) }>'


def get_delphi_param_view_type(cfg):
  n_dim = dim_len(cfg)
  if n_dim == 0:
    return get_delphi_num_type(cfg)
  else:
    return 'P' + get_delphi_num_type(cfg)


def get_delphi_state_type(cfg):
  return f'TGBFixedArray<{ get_delphi_num_type(cfg) }>'


def get_delphi_state_view_type(cfg):
  return 'P' + get_delphi_num_type(cfg)


def get_delphi_num_type(cfg):
  if cfg["num_type"] == 'real_type':
    return 'Double'
  elif cfg["num_type"] == 'int':
    return 'Integer'
  else:
    raise Exception(f'Unknown type {cfg["num_type"]}, must be one of "real_type" or "int".')


def get_delphi_setter(name, cfg):
  name = to_camel_case(name)
  return f'Set{ name }(const in{ name }: TGBFixedArray<{get_delphi_num_type(cfg)}>);'


def get_delphi_param_setter(config_name):
  return f'Set{ config_name }Params(const { config_name.lower() }Params: Leapfrog{ config_name }ParamsView);'


def get_delphi_state_setter(config_name):
  return f'Set{ config_name }State(const { config_name.lower() }State: Leapfrog{ config_name }StateView);'
