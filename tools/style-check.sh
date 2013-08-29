#!/bin/bash 
EXTENSIONS="h cpp hpp"

echo "Checking coding-style..." # using cpplint

find -E ./ -regex '.*\.(hpp|h|cpp)$' -print0 | xargs -0 tools/cpplint.py \
    --counting=detailed \
    --filter=-,\
+readability/namespace,+build/include_what_you_use,+runtime/explicit,\
+runtime/sizeof,+build/include_order,+whitespace/blank_line,+whitespace/labels,\
+whitespace/parens,+whitespace/braces,+whitespace/comma,+whitespace/indent,\
+whitespace/operators,+whitespace/comments,+whitespace/semicolon,\
+whitespace/line_length,+whitespace/indent,+whitespace/ending_newline,\
+whitespace/forcolon,+whitespace/end_of_line,+whitespace/empty_loop_body,\
+whitespace/newline,+whitespace/tab,+whitespace/todo

echo "...checking coding-style done!"
