## Installation and running

Ensure that Python and pip are installed, as well as  C++ compiler. Then, isntall the required python packages with pip:

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
./build/AI-CFD.exe
```

