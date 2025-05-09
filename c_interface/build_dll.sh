#!/usr/bin/env bash

set -e

current="$(pwd)"
echo $current
trap "cd $current" EXIT

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

cd $here && cmake --build build
