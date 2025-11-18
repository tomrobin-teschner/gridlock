![Static Badge](https://img.shields.io/badge/version-0.11.0-red)

> [!WARNING]
> The solver is working for the test cases provided, but it does not yet support all features that are planned for this solver. The readme file will be updated once the solver is in working order and the version will be bumped to 1.x.x once ready for use.

## Installation and running

Ensure that Python and pip are installed, as well as  C++ compiler. Then, install the required python packages with pip:

```bash
pip install -r .\requirements.txt
```

Create a conan profile if it is the first time you are using it. Let Conan detect sensible defaults:

```bash
conan profile detect --force
```

Then, execute Conan to install all required external libraries

```bash
conan install . --build=missing --settings=build_type=Debug
```

Now configure the project

```bash
cmake -B build/ -S . -DCMAKE_BUILD_TYPE=Debug -DCMAKE_TOOLCHAIN_FILE='build/Debug/generators/conan_toolchain.cmake' -G Ninja
```

Compile it

```bash
cmake --build build/ --config Debug
```

Finally, run the code:

```bash
./build/gridlock
```

Or, if you are on Windows, use:

```bash
./build/Release/gridlock.exe
```

### TODO

- update writing to reflect correct paths for Ninja generator
- Currently mixed instructions for Debug/Release
- Add some additional information about solver and use/simple tutorial
- update conan instructions to use C++ 23 (default detects C__14, but TOML needs at least C++17)
