#!/bin/bash 
echo "Configuring Hom3..."

./tools/cleanup.sh

echo "Running cmake..."
while getopts "dramthx" option
do
    case "${option}" in
        h)
            echo "HOM3 build script. Following build modes are available:"
            echo "-d    Debug"
            echo "-r    Release"
            echo "-a    Asan"
            echo "-m    Msan"
            echo "-t    Tsan"
            exit 1
            ;;
        d)
            cmake -DCMAKE_C_COMPILER=clang \
		-DCMAKE_CXX_COMPILER=clang++ \
                -DLIBCXX_INCLUDE_DIR="${LIBCXX_INCLUDE}" \
                -DLIBCXX_LIBRARY="${LIBCXX_LIB}" \
                -DCMAKE_BUILD_TYPE=Debug
            echo "... hom3-debug configuration finished!"
            exit 1
            ;;
        r)
            cmake -DCMAKE_C_COMPILER=/usr/local/bin/clang \
		-DCMAKE_CXX_COMPILER=/usr/local/bin/clang++ \
                -DLIBCXX_INCLUDE_DIR=$LIBCXX_INCLUDE \
                -DLIBCXX_LIBRARY=$LIBCXX_LIB \
                -DCMAKE_BUILD_TYPE=Release
            echo "... hom3-release configuration finished!"
            exit 1
            ;;
        a)
            cmake -DCMAKE_C_COMPILER=/usr/local/bin/clang \
		-DCMAKE_CXX_COMPILER=/usr/local/bin/clang++ \
                -DLIBCXX_INCLUDE_DIR=$LIBCXX_INCLUDE \
                -DLIBCXX_LIBRARY=$LIBCXX_LIB \
                -DCMAKE_BUILD_TYPE=Asan
            echo "... hom3-release configuration finished!"
            exit 1
            ;;
        m)
            cmake -DCMAKE_C_COMPILER=/usr/local/bin/clang \
		-DCMAKE_CXX_COMPILER=/usr/local/bin/clang++ \
                -DLIBCXX_INCLUDE_DIR=$LIBCXX_INCLUDE \
                -DLIBCXX_LIBRARY=$LIBCXX_LIB \
                -DCMAKE_BUILD_TYPE=Msan
            echo "... hom3-release configuration finished!"
            exit 1
            ;;
        t)
            cmake -DCMAKE_C_COMPILER=/usr/local/bin/clang \
		-DCMAKE_CXX_COMPILER=/usr/local/bin/clang++ \
                -DLIBCXX_INCLUDE_DIR=$LIBCXX_INCLUDE \
                -DLIBCXX_LIBRARY=$LIBCXX_LIB \
                -DCMAKE_BUILD_TYPE=Tsan
            echo "... hom3-release configuration finished!"
            exit 1
            ;;
    esac
done

#                -DLIBCXX_INCLUDE_DIR=$LIBCXX_INCLUDE \
#                -DLIBCXX_LIBRARY=$LIBCXX_LIB \
