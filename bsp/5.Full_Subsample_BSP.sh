for j in $(cat ../Groups.list); do 
seqtk subseq ../chr_W.BSP.AllSites.50bp.fa ${j}.pop | seqret -sequence stdin -outseq nex/${j}.50.nex -osformat nexus
seqtk subseq ../chr_W.BSP.AllSites.10bp.fa ${j}.pop | seqret -sequence stdin -outseq nex/${j}.10.nex -osformat nexus
done
