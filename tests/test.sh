#!/bin/sh

for test in `ls test_*.py`
do
  echo $test
  python3 $test || exit 1
done
