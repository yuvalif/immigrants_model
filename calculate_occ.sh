#!/bin/bash

#prog=estimation_test_full_random
prog=estimation_test_full

echo running "$prog"...
filename="tmp.log"
./"$prog" input.txt output.txt > "$filename"
line_number=$(grep -n counter "$filename" | cut -d':' -f 1)
sed "$line_number,$ !d" "$filename"
sed -i "$line_number,$ d" "$filename"

for i in {5..16}
do
  total=$(awk "{if (\$$i >= 0) print 1}" "$filename" | wc -l)
  ue=$(awk "{if (\$$i >=0 && \$$i <= 6) print 1}" "$filename" | wc -l)
  blue=$(awk "{if (\$$i > 6 && \$$i <= 13) print 1}" "$filename" | wc -l)
  white=$(awk "{if (\$$i > 13) print 1}" "$filename" | wc -l)

  echo "$total" "$ue" "$white" "$blue"
done
