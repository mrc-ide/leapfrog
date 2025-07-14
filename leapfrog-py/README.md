To build use:
```
uv build
```

Make sure to have your tooling up to date:
```
python -m pip install --upgrade pip
python -m pip install --upgrade setuptools wheel twine check-wheel-contents
```

UV should build in `dist/` directory so you can install the package via:
```
pip install dist/**.tar.gz
```

Once installed you can use it in your python code:
```python
>>> import leapfrog_py
>>> leapfrog_py.hello()
20
```

We are using `scikit-build-core` with `nanobind` for this project.
