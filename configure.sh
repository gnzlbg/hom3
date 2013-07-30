rm CMakeCache.txt
rm -rf CMakeFiles

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
            cmake -DCMAKE_CXX_COMPILER=/usr/local/bin/clang++ -DLIBCXX_INCLUDE_DIR=$LIBCXX/include -DLIBCXX_LIBRARY=$LIBCXX/lib -DCMAKE_BUILD_TYPE=Debug
            exit 1
            ;;
        r)
            cmake -DCMAKE_CXX_COMPILER=/usr/local/bin/clang++ -DLIBCXX_INCLUDE_DIR=$LIBCXX/include -DLIBCXX_LIBRARY=$LIBCXX/lib -DCMAKE_BUILD_TYPE=Release
            exit 1
            ;;
        x)
            cmake -G Xcode -DCMAKE_CXX_COMPILER=/usr/local/bin/clang++ -DLIBCXX_INCLUDE_DIR=$LIBCXX/include -DLIBCXX_LIBRARY=$LIBCXX/lib -DCMAKE_BUILD_TYPE=Release
            exit 1
            ;;
    esac
done
