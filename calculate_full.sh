#!/bin/bash

prog=estimation_test_fullwage

echo running "$prog"...

filename=tmp.log
./"$prog" input.txt output.txt > "$filename"
line_number=$(grep -n counter "$filename" | cut -d':' -f 1)
sed "$line_number,$ !d" "$filename"
sed -i "$line_number,$ d" "$filename"


printf "\n\nOccupation Distribution Count\n"
printf "Time\tTotal\tUE\tWhite\tFull\tPart\n"

for i in {5..24..2}
do
  total=$(awk "BEGIN{count=0}{if (\$$i >= 0) count=count+1}END{print count}" "$filename")
  ue=$(awk "BEGIN{ue=0}{if (\$$i >=0 && \$$i <= 6) ue=ue+1}END{print ue}" "$filename")
  full=$(awk "BEGIN{full=0}{if (\$$i > 6 && \$$i <= 13) full=full+1}END{print full}" "$filename")
  part=$(awk "BEGIN{part=0}{if (\$$i > 13 && \$$i <= 20) part=part+1}END{print part}" "$filename")
  white=$(awk "BEGIN{white=0}{if (\$$i > 20) white=white+1}END{print white}" "$filename")

  printf "%d\t%d\t%d\t%d\t%d\t%d\n" $((i-5)) "$total" "$ue" "$white" "$full" "$part"
done
printf "\n\nOccupation Distribution Percent\n"
printf "Time\tTotal\tUE\tWhite\tFull\tPart\n"

for i in {5..24..2}
do
  total=$(awk "BEGIN{count=0}{if (\$$i >= 0) count=count+1}END{print count}" "$filename")
  ue=$(awk "BEGIN{count=0;ue=0}{if (\$$i >= 0) count=count+1; if (\$$i >=0 && \$$i <= 6) ue=ue+1}END{printf(\"%.4f\n\", ue/count)}" "$filename")
  full=$(awk "BEGIN{count=0;full=0}{if (\$$i >= 0) count=count+1; if (\$$i > 6 && \$$i <= 13) full=full+1}END{printf(\"%.4f\n\", full/count)}" "$filename")
  part=$(awk "BEGIN{count=0;part=0}{if (\$$i >= 0) count=count+1; if (\$$i > 13 && \$$i <= 20) part=part+1}END{printf(\"%.4f\n\", part/count)}" "$filename")
  white=$(awk "BEGIN{count=0;white=0}{if (\$$i >= 0) count=count+1; if (\$$i > 20) white=white+1}END{printf(\"%.4f\n\", white/count)}" "$filename")

  printf "%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\n" $((i-5)) "$total" "$ue" "$white" "$full" "$part"
done

printf "\n\nWhite Wage in Last Period\n"
printf "Type\t1\t3\t4\t5\t267\tAverage\n"

for type in {0..2}
do
  average=$(awk -v type="$type" 'BEGIN{sum=0; count=0}{if ($3 == type && $(NF-1) > 20) {count=count+1; sum=sum+$NF}}END{print sum/count}' "$filename")

  rg1=$(awk -v type="$type" 'BEGIN{sum=0; count=0}{if ($3 == type && $(NF-1) > 20 && $(NF-1)%7 == 0) {count=count+1; sum=sum+$NF}}END{print sum/count}' "$filename")
  rg3=$(awk -v type="$type" 'BEGIN{sum=0; count=0}{if ($3 == type && $(NF-1) > 20 && $(NF-1)%7 == 2) {count=count+1; sum=sum+$NF}}END{print sum/count}' "$filename")
  rg4=$(awk -v type="$type" 'BEGIN{sum=0; count=0}{if ($3 == type && $(NF-1) > 20 && $(NF-1)%7 == 3) {count=count+1; sum=sum+$NF}}END{print sum/count}' "$filename")
  rg5=$(awk -v type="$type" 'BEGIN{sum=0; count=0}{if ($3 == type && $(NF-1) > 20 && $(NF-1)%7 == 4) {count=count+1; sum=sum+$NF}}END{print sum/count}' "$filename")
  rg267=$(awk -v type="$type" 'BEGIN{sum=0; count=0}{if ($3 == type && $(NF-1) > 20 && ($(NF-1)%7 == 1 || $(NF-1)%7 == 5 || $(NF-1)%7 == 6)) {count=count+1; sum=sum+$NF}}END{print sum/count}' "$filename")
  printf "%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n" "$type" "$rg1" "$rg3" "$rg4" "$rg5" "$rg267" "$average"
done

printf "\n\nBlue Full Wage in Last Period\n"
printf "Type\t1\t2\t3\t4\t5\t6\t7\tAverage\n"

for type in {0..2}
do
  average=$(awk -v type="$type" 'BEGIN{sum=0; count=0}{if ($3 == type && $(NF-1) > 6 && $(NF-1) <= 13) {count=count+1; sum=sum+$NF}}END{print sum/count}' "$filename")

  printf "%d\t" "$type"
  for rg in {0..6}
  do
    wage=$(awk -v type="$type" -v rg="$rg" 'BEGIN{sum=0; count=0}{if ($3 == type && $(NF-1) > 6 && $(NF-1) <= 13 && $(NF-1)%7 == rg) {count=count+1; sum=sum+$NF}}END{print sum/count}' "$filename")
    printf "%.2f\t" "$wage"
  done
  printf "%.2f\n" "$average"
done

printf "\n\nBlue Part Wage in Last Period\n"
printf "Type\t1\t2\t3\t4\t5\t6\t7\tAverage\n"

for type in {0..2}
do
  average=$(awk -v type="$type" 'BEGIN{sum=0; count=0}{if ($3 == type && $(NF-1) > 13 && $(NF-1) <= 20) {count=count+1; sum=sum+$NF}}END{print sum/count}' "$filename")

  printf "%d\t" "$type"
  for rg in {0..6}
  do
    wage=$(awk -v type="$type" -v rg="$rg" 'BEGIN{sum=0; count=0}{if ($3 == type && $(NF-1) > 13 && $(NF-1) <= 20 && $(NF-1)%7 == rg) {count=count+1; sum=sum+$NF}}END{print sum/count}' "$filename")
    printf "%.2f\t" "$wage"
  done
  printf "%.2f\n" "$average"
done

