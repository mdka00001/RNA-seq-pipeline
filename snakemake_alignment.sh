
#Merge feature counts
TARGET_DIR="$1"
ls -1 "$TARGET_DIR"featureCounts_*.txt | parallel --env TARGET_DIR 'cat {} | sed "1d" | cut -f7 > "$TARGET_DIR"{/.}_filtered.txt'


ls -1  $1featureCounts_*.txt | head -1 | xargs cut -f1 > $1genes.txt
paste $1genes.txt $1*_filtered.txt > $1output.txt
