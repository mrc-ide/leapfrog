import json
import os

from jinja2 import Environment, FileSystemLoader, select_autoescape

import utils
import utils.concepts
import utils.config
import utils.state_space


def relative_file_path(path):
  dirname = os.path.dirname(__file__)
  return os.path.join(dirname, path)


def load_json(path):
  with open(relative_file_path(path)) as f:
    return json.load(f)


def load_children_model_schemas(paths):
  if isinstance(paths, str):
    return load_json(f'../modelSchemas/{paths}')
  else:
    return [load_json(f'../modelSchemas/{p}') for p in paths]


def generate_hpp(template_name, *args, **kwargs):
  template = env.get_template(f'{template_name}.j2')
  output = template.render(*args, **kwargs)
  with open(relative_file_path(f'../../inst/include/generated/{template_name}.hpp'), "w") as f:
    f.write(output)


dat = load_json("../modelSchemas/FullModel.json")
dat = { k: load_children_model_schemas(v) for k, v in dat.items() }

file_loader = FileSystemLoader(relative_file_path("../templates"))
env = Environment(
  loader = file_loader,
  trim_blocks = True,
  lstrip_blocks = True,
  autoescape = select_autoescape(),
  keep_trailing_newline = True
)

generate_hpp("model_variants", dat)
generate_hpp("state_space", dat | vars(utils.state_space))
generate_hpp("concepts", dat | vars(utils.concepts))
generate_hpp("state_space_mixer", dat)
generate_hpp("config", dat | vars(utils.config))
generate_hpp("config_mixer", dat)
