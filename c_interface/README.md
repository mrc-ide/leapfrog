# leapfrog c

This dir contains a interface to compile a C-compatible DLL. We use it to call the model from Delphi. You can build the DLL using the helper script

```
./build_dll.sh
```

Or manually

```
cmake -A Win32 -B build
cmake --build build
```
