bowtie2-build index/viral.fasta index/viral

mkdir virus_bam
find clean/*_R1_paired*|sed 's/clean\///'|sed 's/_R1.*//'|parallel --gnu "bowtie2 -p 5 -x index/viral -1 clean/{}_R1_paired.fastq.gz -2 clean/{}_R2_paired.fastq.gz|samtools view -Sh -q 30 -F 4 -|grep -v 'XS:'|samtools view -Shub|samtools sort - -o virus_bam/{}.bam"
find virus_bam/*.bam|parallel --gnu "samtools index {} {}.bai"
find virus_bam/*.bam|parallel --gnu "samtools view -c -F 4 {} > {}.cnt"
find virus_bam/*.bam|parallel --gnu "samtools view -hub -F 16 {}|samtools depth -d 1000000 - > {}.cov"
head virus_bam/*cnt|sed 's/==> virus_bam\///'|sed -r 's/.bam.*<==//' > virus_bam/mapping_stat

mkdir virus_track
python script/cal_virus_track.py virus_bam/EV-1_a-H3.bam.cov virus_bam/EV-1_Input.bam.cov virus_track/H3_vs_input_EV.1.track
python script/cal_virus_track.py virus_bam/EV-1_a-H3K27me3.bam.cov virus_bam/EV-1_Input.bam.cov virus_track/H3K27me3_vs_input_EV.1.track
python script/cal_virus_track.py virus_bam/EV-1_a-H3K27me3.bam.cov virus_bam/EV-1_a-H3.bam.cov  virus_track/H3K27me3_vs_H3_EV.1.track
python script/cal_virus_track.py virus_bam/TYLCV-1_a-H3.bam.cov virus_bam/TYLCV-1_Input.bam.cov virus_track/H3_vs_input_TYLCV.1.track
python script/cal_virus_track.py virus_bam/TYLCV-1_a-H3K27me3.bam.cov virus_bam/TYLCV-1_Input.bam.cov virus_track/H3K27me3_vs_input_TYLCV.1.track
python script/cal_virus_track.py virus_bam/TYLCV-1_a-H3K27me3.bam.cov virus_bam/TYLCV-1_a-H3.bam.cov virus_track/H3K27me3_vs_H3_TYLCV.1.track
python script/cal_virus_track.py virus_bam/TYLCV-C3-mut-1_a-H3.bam.cov virus_bam/TYLCV-C3-mut-1_Input.bam.cov virus_track/H3_vs_input_TYLCV-C3mut.1.track
python script/cal_virus_track.py virus_bam/TYLCV-C3-mut-1_a-H3K27me3.bam.cov virus_bam/TYLCV-C3-mut-1_Input.bam.cov virus_track/H3K27me3_vs_input_TYLCV-C3mut.1.track
python script/cal_virus_track.py virus_bam/TYLCV-C3-mut-1_a-H3K27me3.bam.cov virus_bam/TYLCV-C3-mut-1_a-H3.bam.cov virus_track/H3K27me3_vs_H3_TYLCV-C3mut.1.track


perl cgview/cgview_xml_builder/cgview_xml_builder.pl -sequence index/sequence.gbk -title 'Tomato yellow leaf curl virus-[Almeria]' -output virus_track/EV.1.xml -tick_density 0.7 -feature_labels T -gc_content F -gc_skew F -analysis virus_track/H3_vs_input_EV.1.track virus_track/H3K27me3_vs_input_EV.1.track virus_track/H3K27me3_vs_H3_EV.1.track
perl cgview/cgview_xml_builder/cgview_xml_builder.pl -sequence index/sequence.gbk -title 'Tomato yellow leaf curl virus-[Almeria]' -output virus_track/TYLCV.1.xml -tick_density 0.7 -feature_labels T -gc_content F -gc_skew F -analysis virus_track/H3_vs_input_TYLCV.1.track virus_track/H3K27me3_vs_input_TYLCV.1.track virus_track/H3K27me3_vs_H3_TYLCV.1.track
perl cgview/cgview_xml_builder/cgview_xml_builder.pl -sequence index/sequence.gbk -title 'Tomato yellow leaf curl virus-[Almeria]' -output virus_track/TYLCV-C3mut.1.xml -tick_density 0.7 -feature_labels T -gc_content F -gc_skew F -analysis virus_track/H3_vs_input_TYLCV-C3mut.1.track virus_track/H3K27me3_vs_input_TYLCV-C3mut.1.track virus_track/H3K27me3_vs_H3_TYLCV-C3mut.1.track

java -jar -Xmx1500m cgview/cgview.jar -i virus_track/EV.1.xml -o figure/EV.1.png -f png
java -jar -Xmx1500m cgview/cgview.jar -i virus_track/TYLCV.1.xml -o figure/TYLCV.1.png -f png
java -jar -Xmx1500m cgview/cgview.jar -i virus_track/TYLCV-C3mut.1.xml -o figure/TYLCV-C3mut.1.png -f png


#mkdir virus_peak
#macs14 -t virus_bam/EV-1_a-H3.bam       -c virus_bam/EV-1_Input.bam -g 2969810994 -p 1e-2 -n virus_peak/H3_vs_input_EV.1 --nomodel --keep-dup=all
#macs14 -t virus_bam/EV-1_a-H3K27me3.bam -c virus_bam/EV-1_Input.bam -g 2969810994 -p 1e-2 -n virus_peak/H3K27me3_vs_input_EV.1 --nomodel --keep-dup=all
#macs14 -t virus_bam/EV-1_a-H3K27me3.bam -c virus_bam/EV-1_a-H3.bam  -g 2969810994 -p 1e-2 -n virus_peak/H3K27me3_vs_H3_EV.1 --nomodel --keep-dup=all
#macs14 -t virus_bam/TYLCV-1_a-H3.bam       -c virus_bam/TYLCV-1_Input.bam -g 2969810994 -p 1e-2 -n virus_peak/H3_vs_input_TYLCV.1 --nomodel --keep-dup=all 
#macs14 -t virus_bam/TYLCV-1_a-H3K27me3.bam -c virus_bam/TYLCV-1_Input.bam -g 2969810994 -p 1e-2 -n virus_peak/H3K27me3_vs_input_TYLCV.1 --nomodel --keep-dup=all
#macs14 -t virus_bam/TYLCV-1_a-H3K27me3.bam -c virus_bam/TYLCV-1_a-H3.bam  -g 2969810994 -p 1e-2 -n virus_peak/H3K27me3_vs_H3_TYLCV.1 --nomodel --keep-dup=all
#macs14 -t virus_bam/TYLCV-C3-mut-1_a-H3.bam       -c virus_bam/TYLCV-C3-mut-1_Input.bam -g 2969810994 -p 1e-2 -n virus_peak/H3_vs_input_TYLCV-C3mut.1 --nomodel --nolambda --keep-dup=all
#macs14 -t virus_bam/TYLCV-C3-mut-1_a-H3K27me3.bam -c virus_bam/TYLCV-C3-mut-1_Input.bam -g 2969810994 -p 1e-2 -n virus_peak/H3K27me3_vs_input_TYLCV-C3mut.1 --nomodel --nolambda --keep-dup=all
#macs14 -t virus_bam/TYLCV-C3-mut-1_a-H3K27me3.bam -c virus_bam/TYLCV-C3-mut-1_a-H3.bam  -g 2969810994 -p 1e-2 -n virus_peak/H3K27me3_vs_H3_TYLCV-C3mut.1 --nomodel --nolambda --keep-dup=all
