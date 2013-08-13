# Removes all cmake-generated files
# Note: it doesn't remove any binaries

TO_REMOVE="Doxyfile CMakeCache.txt CMakeFiles CTestTestfile.cmake Makefile cmake_install.cmake"

for FILE in $TO_REMOVE
do
    echo "removing ${FILE} ..."
    find . -name $FILE -print0 | xargs -0 rm -rf
done

# Clean GoogleTest stuff
echo "removing ext/gtest/src..."
rm -rf ext/gtest/src
echo "removing ext/gtest/tmp..."
rm -rf ext/gtest/tmp
