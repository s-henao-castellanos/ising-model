for file in *.pdf; do
    magick "$file" "${file%.pdf}.png";
done