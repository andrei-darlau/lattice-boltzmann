Compile:
```
cmake -B build
cmake --build build
```

Run:
```
cd builddir
./HelloKokkos
```

For nvim:
```
cmake -S . -B build -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
```
