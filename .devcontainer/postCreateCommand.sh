#!/usr/bin/env bash
# install JuliaFormatter in global environment for VSCode extension
julia --project=/workdir -e 'using Pkg; Pkg.instantiate();'
julia --project=/workdir -e 'using Pkg; Pkg.activate(); Pkg.add("JuliaFormatter");'
# show versions
julia --version
python --version
R --version
