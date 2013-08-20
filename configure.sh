#!/bin/bash 
echo "Configuring Hom3..."

./tools/cleanup.sh

echo "Running cmake..."
while getopts "drhx" option
do
    case "${option}" in
        h)
            echo "HOM3 build script. Following build modes are available:"
            echo "-d    Debug"
            echo "-r    Release"
            exit 1
            ;;
        d)
            cmake -DCMAKE_C_COMPILER=/usr/local/bin/clang \
		-DCMAKE_CXX_COMPILER=/usr/local/bin/clang++ \
                -DLIBCXX_INCLUDE_DIR=$LIBCXX/include \
                -DLIBCXX_LIBRARY=$LIBCXX/lib \
                -DCMAKE_BUILD_TYPE=Debug
            echo "... hom3-debug configuration finished!"
            exit 1
            ;;
        r)
            cmake -DCMAKE_C_COMPILER=/usr/local/bin/clang \
		-DCMAKE_CXX_COMPILER=/usr/local/bin/clang++ \
                -DLIBCXX_INCLUDE_DIR=$LIBCXX/include \
                -DLIBCXX_LIBRARY=$LIBCXX/lib \
                -DCMAKE_BUILD_TYPE=Release
            echo "... hom3-release configuration finished!"
            exit 1
            ;;
    esac
done
