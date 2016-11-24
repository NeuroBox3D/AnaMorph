------------
all builds require llvm/clang >= 3.4 for full C++11 standard support
along with a compatible C++ STL implementation. all makefiles are
adjusted to use the GNU STL libstdc++. all builds tested and working
with libstdc++ version >= 4.7.2.5-5
------------

1. Dynamically linked executables, arch depending on host

make clean -f Makefile.dynamic
make -j4 -f Makefile.dynamic

2. Statically linked executables, arch depending on host

make clean -f Makefile.static
make -j4 -f Makefile.static

3. Statically linked x86 cross-compiled executables (requires
cross-compilation/multilib support).

make clean -f Makefile.static.x86
make -j4 -f Makefile.static.x86
