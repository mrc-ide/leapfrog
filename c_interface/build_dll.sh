set -e

current="$(pwd)"
echo $current
trap "cd $current" EXIT

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd $here

cmake --build build

# TODO: Copying into Rob's specific directory structure and won't work for everyone
# let's think about how to manage this more generally
echo "Copying DLL into place"
cp build/Debug/leapfrog.dll $here/../../Spec5/
