

ls "$1"/results/depth/*.depth.txt | sed 's,.txt$,,' > "$1"/results/depth/temp_fofn

while read i 
do

	csvtk plot line \
		-t \
		-x position \
		-y depth \
		--point-size 0.01 \
		-o "$i".pdf \
		"$i".txt

done < "$1"/results/depth/temp_fofn

rm "$1"/results/depth/temp_fofn