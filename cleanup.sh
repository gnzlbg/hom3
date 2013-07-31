# Removes all cmake-generated files
# Note: it doesn't remove any binaries

TO_REMOVE="Doxyfile CMakeCache.txt CMakeFiles CTestTestfile.cmake Makefile cmake_install.cmake"

for FILE in $TO_REMOVE
do
    echo "removing ${FILE} ..."
    find . -name $FILE -print0 | xargs -0 rm -rf
done
