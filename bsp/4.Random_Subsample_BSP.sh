for j in $(cat ../Groups.list); do 
cat ${j}.pop | shuf -n 5 > sub.tmp
seqtk subseq ../chr_W.BSP.AllSites.50bp.fa sub.tmp | seqret -sequence stdin -outseq nex/${j}.50.nex -osformat nexus
seqtk subseq ../chr_W.BSP.AllSites.10bp.fa sub.tmp | seqret -sequence stdin -outseq nex/${j}.10.nex -osformat nexus
done
