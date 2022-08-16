#!/usr/bin/ksh

for dir in $(<$1); do
  print ${dir}
  cp -R ${dir} gastdocked_complex_folders/ &
done
