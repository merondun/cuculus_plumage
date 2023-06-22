blastn -query GCA_017976375.1_bCucCan1.pri_mtDNA.fa -db blastdbs/Cuculus_canorus -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand' -num_threads 5 | \
	awk '{OFS="\t"; if ($13 == "minus") print $2,$10,$9,$3,$13; else print $2,$9,$10,$3,$13}' | \   #fix strand coordinates
    sed 's/minus/-/g' | sed 's/plus/+/g' | \ 
    bedtools sort -i - | \
    bedtools merge -i - -d 500 -c 4,5 -o mean,distinct | \ #average % identity / strand across hits within 500bp
    awk '{OFS="\t"}{print $0, $3-$2+1}'  #calculate length and convert from 0-based
