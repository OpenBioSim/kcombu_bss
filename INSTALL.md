# Installation Instructions

## Dependencies

* A working C compiler, e.g. gcc, clang, MSVC etc.
* cmake

## Instructions

Create a build directory, e.g.

```
$ mkdir build
```

Change into that build directory, e.g.

```
$ cd build
```

Run `cmake`, giving the path to the `src` directory, e.g.

```
$ cmake ../src
```

Run `make` to compile all of the programs. You can compile
in parallel using the `-j` option, e.g.

```
$ make -j 8
```

Install via

```
$ make install
```

This will install into the default `cmake` run directory (e.g. `/usr/local/bin`).

You can control the install location using either `ccmake` to configure
manually, or by passing this in as an option to `cmake --install`, e.g.

```
$ cmake --install . --prefix /usr
```

would install into `/usr/bin`.
