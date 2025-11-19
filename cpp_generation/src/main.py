import json
import os
import copy

from jinja2 import Environment, FileSystemLoader, select_autoescape

import utils.concepts
import utils.config
import utils.state_space
import utils.delphi
import utils.general


def relative_file_path(*paths):
  dirname = os.path.dirname(__file__)
  return os.path.join(dirname, *paths)


def load_json(*paths):
  with open(relative_file_path(*paths)) as f:
    return json.load(f)


def load_children_model_schemas(paths):
  if isinstance(paths, str):
    return load_json("..", "modelSchemas", paths)
  else:
    return [load_json("..", "modelSchemas", p) for p in paths]


def generate(template_path, dest_path, *args, **kwargs):
  template = env.get_template(template_path)
  output = template.render(*args, **kwargs)
  with open(dest_path, "w") as f:
    f.write(output)

supported_languages = ["r", "python"]
supported_num_types = ["real_type", "int"]

def process_var_config(name, cfg, is_par = False):
  if not cfg.get("num_type"):
    unsupported_num_type_err_msg = f'num_type for {name} should be one of {" or ".join(supported_num_types)}'
    assert cfg["num_type"] in supported_num_types, unsupported_num_type_err_msg

  if not cfg.get("dims"):
    cfg["type"] = "scalar"
  else:
    dims_not_list_err_msg = f'dims for {name} should be a list'
    assert isinstance(cfg["dims"], list), dims_not_list_err_msg
    cfg["type"] = "tensor"

  if not is_par: return
  if not cfg.get("alias"):
    cfg["alias"] = { l: name for l in supported_languages }
  else:
    alias = cfg["alias"]
    for l in supported_languages:
      if not alias.get(l):
        alias[l] = name
      else:
        alias_not_string_err_msg = f'{l} alias for {name} should be a string'
        assert isinstance(alias[l], str), alias_not_string_err_msg


def add_output_year_dim(cfg):
  if not cfg.get("dims"):
    cfg["dims"] = ["output_years"]
  else:
    cfg["dims"].append("output_years")


def generate_hpp(template_name, *args, **kwargs):
  template_path = f'cpp/{template_name}.j2'
  dest_path = relative_file_path("..", "..", "leapfrogr", "inst", "include", "generated", f"{template_name}.hpp")
  generate(template_path, dest_path, *args, **kwargs)


def generate_delphi(template_name, *args, **kwargs):
  template_path = f'delphi/{template_name}.j2'
  dest_path = relative_file_path("..", "..", "delphi", f"{template_name}.pas")
  generate(template_path, dest_path, *args, **kwargs)


dat = load_json("..", "modelSchemas", "FullModel.json")
dat = { k: load_children_model_schemas(v) for k, v in dat.items() }


for config in dat["configs"]:
  config["output_state"] = copy.deepcopy(config["state"])

  for name, cfg in config["pars"].items():
    process_var_config(name, cfg, True)

  for name, cfg in config["intermediate"].items():
    process_var_config(name, cfg)

  for name, cfg in config["state"].items():
    process_var_config(name, cfg)

  for name, cfg in config["output_state"].items():
    add_output_year_dim(cfg)
    process_var_config(name, cfg)


file_loader = FileSystemLoader(relative_file_path("..", "templates"))
env = Environment(
  loader = file_loader,
  trim_blocks = True,
  lstrip_blocks = True,
  autoescape = select_autoescape(),
  keep_trailing_newline = True
)

generate_hpp("model_variants", dat)
generate_hpp("state_space", utils.state_space.extract_state_space_info(dat) | vars(utils.state_space))
generate_hpp("concepts", dat | vars(utils.concepts))
generate_hpp("config", dat | vars(utils.config) | vars(utils.general))
generate_hpp("config_mixer", dat)
generate_hpp("py_interface/py_adapters", dat | vars(utils.config) | vars(utils.general))
generate_hpp("cpp_interface/hdf5_utils", dat | vars(utils.config) | vars(utils.general))
generate_hpp("cpp_interface/cpp_adapters", dat | vars(utils.config) | vars(utils.general))
generate_hpp("r_interface/r_adapters", dat | vars(utils.config) | vars(utils.general))
generate_hpp("c_interface/c_adapters", dat | vars(utils.config) | vars(utils.general))
generate_hpp("c_interface/c_types", dat | vars(utils.config) | vars(utils.delphi) | vars(utils.general))
generate_delphi("LeapfrogInterface", dat | vars(utils.delphi) | vars(utils.general))
