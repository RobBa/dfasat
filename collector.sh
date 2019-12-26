HEADER=source/evaluators.h
rm $HEADER
echo "#ifndef __ALL_HEADERS__" > $HEADER
echo "#define __ALL_HEADERS__" >> $HEADER
for file in source/evaluation/*.h
do
    echo "#include \"$file\"" >> $HEADER
done
echo "#endif" >> $HEADER
