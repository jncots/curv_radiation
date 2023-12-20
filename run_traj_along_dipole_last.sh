#!/usr/bin/bash
(set -o igncr) 2>/dev/null && set -o igncr; # this comment is needed for Windows system

# Main script
lib_mount="/home/antonpr/work/pulsar_rad/retro_code/mount_library.sh"
# Launching in current directory
cur_dir=$(pwd)
source_dir="source_dir=$cur_dir/src/fortran"
out_dir="out_dir=$cur_dir/bin"

# Launching in other directories
#dir_scr='dir_scr=/home/anton/work/project/halo/chp_prop1/'
#dir_out='dir_out=/home/anton/work/project/halo/chp_prop1/'
# Source files and Output executable file
source_files="source_files=traj_along_dipole_last_gam.f90  util_test.f90"
out_file="out_file=test_traj_along_dipole"

# File with parameters
param="param=file_with_param.dat"

# Launching in separate window:
#x-terminal-emulator -e $lib_mount "$source_dir" "$source_files" "$out_dir" "$out_file" "$param"

#echo $lib_mount "$source_dir" "$source_files" "$out_dir" "$out_file" "$param"
# Launching in current window:
$lib_mount "$source_dir" "$source_files" "$out_dir" "$out_file" "$param"
