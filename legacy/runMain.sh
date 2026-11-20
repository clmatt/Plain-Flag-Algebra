export PATH=$./mosek/9.3/tools/platform/linux64x86/bin/mosek:$PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:./mosek/9.3/tools/platform/linux64x86/bin
export export MOSEKLM_LICENSE_FILE=./mosek.lic
make install -C ./mosek/9.3/tools/platform/linux64x86/src/fusion_cxx
g++ -pthread -std=c++11 -g -I./mosek/9.3/tools/platform/linux64x86/h -I./mosek/9.3/tools/platform/linux64x86/include -I./nauty27r2 -L./mosek/9.3/tools/platform/linux64x86/bin -Wl,-rpath-link,./mosek/9.3/tools/platform/linux64x86/bin '-Wl,-rpath=$$ORIGIN/mosek/9.3/tools/platform/linux64x86/bin' -o a.out main.cpp -lfusion64 -lmosek64 ./nauty27r2/nauty.a
./a.out
