#!/bin/env bash

#--------------------------------------------------------------------------------------------------
# Helper functions
echo_usage() {
    echo "Usage:"
    echo "    $0 <EquationSystem> <Mesh> <2D/3D> <ExampleName> <-n num_MPI> <-b build_dir>"
}

execute() {
    local run_cmd=$1
    local run_dir=$2
    echo "----------------------------------------------------------------------------------------------------"
    echo "Executing [$run_cmd] in [$run_dir]"
    cd "$run_dir" 
    eval "$run_cmd"
    cd -
    echo "----------------------------------------------------------------------------------------------------"
}

generate_run_dir() {
    local eg_dir="$1"
    local run_dir="$2"
    
    if [ -e "$run_dir" ]; then
        read -p "Overwrite existing run directory at $run_dir? (Y/N): " choice && [[ $choice == [yY] || $choice == [yY][eE][sS] ]] || exit 5
        \rm -rf "$run_dir"
    fi
    mkdir -p "$(dirname $run_dir)"
    cp -r "$eg_dir" "$run_dir"
}

set_default_build_dir() {
    build_dir=$(find -L "$REPO_ROOT" -maxdepth 2 -type d -name "spack-build*" -printf "%TY-%Tm-%Td %TT %p\n" | sort -n|tail -1|cut -d " " -f 3)
}

parse_args() {
    if [ $# -lt 2 ]; then
        echo_usage
        exit 1
    fi
    POSITIONAL_ARGS=()
    while [[ $# -gt 0 ]]; do
    case $1 in
        -b|--build-dir)
        if [ -z "$2" ]; then
            echo "Error: Empty build dir supplied"
            exit 2
        fi
        build_dir=$(realpath "$2")
        shift 2
        ;;
        -n|--num_mpi)
        nmpi="$2"
        shift 2
        ;;
        -*|--*)
        echo "Error: Unknown option $1"
        exit 3
        ;;
        *)
        # Save positional args in an array
        POSITIONAL_ARGS+=("$1")
        shift
        ;;
    esac
    done

    # Restore and extract positional args
    set -- "${POSITIONAL_ARGS[@]}"
    eq_sys=$1
    mesh=$2
    dim=$3
    eg_name=$4
}

report_options() {
    echo "--------------------------------------------------"
    echo "Options:"
    echo "      EquationSystem : $eq_sys"
    echo "      Mesh : $mesh"
    echo "      Dimension : $dim"
    echo "      Example : $eg_name"
    echo "      n MPI : $nmpi"
    echo "--------------------------------------------------"
}

set_run_cmd() {
    local run_dir=$1
    local solver_exec="$2"
    local nmpi="$3"
    run_cmd_file="$run_dir/run_cmd_template.txt"
    if [ -f "$run_cmd_file" ]; then
        run_cmd=$(sed -e 's|<SOLVER_EXEC>|'"$solver_exec"'|g' -e 's|<NMPI>|'"$nmpi"'|g'< "$run_cmd_file")
        \rm "$run_cmd_file"
    else
        echo "Error: Can't read template run command from $run_cmd_file."
        exit 4
    fi
}

validate_paths() {
    local build_dir=$1
    local solver_exec=$2

    if [ -z "$build_dir" ]; then
        # default build dir wasn't found and none was supplied
        echo "Error: No build directory at ./build-*/spack-build-*"
        echo "Set solver location with '-b <build_dir>'"
        echo "--------------------------------------------------"
        exit 5
    fi

    if [ ! -f "$solver_exec" ]; then
        echo "Error: No solver found at $solver_exec"
        echo "--------------------------------------------------"
        exit 6
    fi

}

validate_eg() {
    local eq_sys=$1
    local mesh=$2
    local dim=$3
    local eg_name=$4

    if [ ! -d "$REPO_ROOT/examples/$eq_sys" ]; then
        echo "Error: No examples found for $eq_sys"
        echo
        echo "Valid equation systems are:"
        echo "$(find $REPO_ROOT/examples -maxdepth 1 -mindepth 1 -type d -exec basename {} \;)"
        echo "--------------------------------------------------"
        exit 7
    fi
    if [ ! -d "$REPO_ROOT/examples/$eq_sys/$mesh" ]; then
        echo "Error: No examples found for $eq_sys/$mesh"
        echo
        echo "Valid meshes for $eq_sys are:"
        echo "$(find $REPO_ROOT/examples/$eq_sys -maxdepth 1 -mindepth 1 -type d -exec basename {} \;)"
        echo "--------------------------------------------------"
        exit 8
    fi
    if [ ! -d "$REPO_ROOT/examples/$eq_sys/$mesh/$dim" ]; then
        echo "Error: No examples found for $eq_sys/$mesh/$dim"
        echo
        echo "Valid dimensions for $eq_sys/$mesh are:"
        echo "$(find $REPO_ROOT/examples/$eq_sys/$mesh -maxdepth 1 -mindepth 1 -type d -exec basename {} \;)"
        echo "--------------------------------------------------"
        exit 9
    fi
    if [ ! -d "$REPO_ROOT/examples/$eq_sys/$mesh/$dim/$eg_name" ]; then
        echo "Error: No examples found for $eq_sys/$mesh/$dim/$eg_name"
        echo
        echo "Valid example cases for $eq_sys/$mesh/$dim are:"
        echo "$(find $REPO_ROOT/examples/$eq_sys/$mesh/$dim -maxdepth 1 -mindepth 1 -type d -exec basename {} \;)"
        echo "--------------------------------------------------"
        exit 10
    fi
}
#--------------------------------------------------------------------------------------------------

REPO_ROOT=$( cd -- "$(realpath $( dirname -- "${BASH_SOURCE[0]}" )/..)" &> /dev/null && pwd )

# Default options
solver_name="penknife"
eq_sys='Not set'
mesh='Not set'
dim='Not set'
eg_name='Not set'
nmpi='4'
build_dir='Not set'
set_default_build_dir

# Parse command line args and report resulting options
parse_args $*
report_options

# Set paths to the solver executable and example directory
solver_exec="$build_dir/$solver_name"
# Validate exec, examples paths
validate_paths "$build_dir" "$solver_exec"
validate_eg "$eq_sys" "$mesh" "$dim" "$eg_name"

# Set up run directory, confirming overwrite if it already exists
eg_dir="$REPO_ROOT/examples/$eq_sys/$mesh/$dim/$eg_name"
run_dir="$REPO_ROOT/runs/$eq_sys/$mesh/$dim/$eg_name"
generate_run_dir "$eg_dir" "$run_dir"

# Read run command template and populate it
run_cmd="Not set"
set_run_cmd "$run_dir" "$solver_exec" "$nmpi"

# Execute run_cmd in run_dir
execute "$run_cmd" "$run_dir"
