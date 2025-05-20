#!/usr/bin/env bash

set -e

usage() {
  echo "Build the DLL for Win32 using CMake, optionally run the configure step"
  echo "Usage: $0 [OPTIONS]"
  echo "Options:"
  echo " -h, --help            Show help"
  echo " --configure           If set, CMake configuration step will be run, required if this is your first time"
}

CONFIGURE=0

handle_options() {
  while [ $# -gt 0 ]; do
    case $1 in
      -h | --help)
        usage
        exit 0
      ;;
      --configure)
        CONFIGURE=1
        ;;
      *)
        echo "Invalid option: $1" >&2
        usage
        exit 1
        ;;
    esac
    shift
  done
}

handle_options "$@"

current="$(pwd)"
echo $current
trap "cd $current" EXIT

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ ! -d "$here/build" ]; then
  if [ "$CONFIGURE" -eq 0 ]; then
    echo "'build' directory not found, CMake has not been configured"
    read -r -p "Would you like to run the CMake configure step now? [y/N] " response
    case "$response" in
      [yY][eE][sS]|[yY])
        CONFIGURE=1
        ;;
      *)
        echo "Aborting build. Run again with --configure or answer 'yes' to the prompt."
        exit 1
        ;;
    esac
  fi
fi

if [ "$CONFIGURE" -eq 1 ]; then
  echo "Running CMake configure step..."
  cd "$here" && cmake -A Win32 -B build
fi

cd "$here" && cmake --build build
