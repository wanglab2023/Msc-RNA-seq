STAR --runThreadN N  --runMode alignReads --genomeDir $starIndex \
--readFilesIn ${sample}_2.fq.extracted  \
--outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 41143265264 \
--outFileNamePrefix $sample.

featureCounts -t exon -g gene_id -a $gtf -o gene_assigned_gene -R BAM $sample.Aligned.sortedByCoord.out.bam -T N
mv $sample.Aligned.sortedByCoord.out.bam.featureCounts.bam $sample.Aligned.sortedByCoord.featureCountsUniq.bam

samtools sort -@ N $sample.Aligned.sortedByCoord.featureCountsUniq.bam -o $sample.Aligned.sortedByCoord.featureCountsUniq.sorted.bam

samtools index $sample.Aligned.sortedByCoord.featureCountsUniq.sorted.bam

umi_tools count --per-gene --gene-tag=XT --per-cell --wide-format-cell-counts -I $sample.Aligned.sortedByCoord.featureCountsUniq.sorted.bam -S $sample.counts.tsv

python selectResult.py $sample.counts.tsv $thres $sample.selected

