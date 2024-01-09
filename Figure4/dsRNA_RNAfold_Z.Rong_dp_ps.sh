#!/bash/bin


# After extracting individual fasta files (using bedtools getfasta) based on chromosome coordinates of target sequences
# save all fasta files in one folder and change directory to that fasta folder, and run:
# the outputs are RNA secondary structures saved in ss.ps files and dp.ps files are used to calculate for MFEs later
for f in *.fa
	do RNAfold -p $f
done

# MFE values are calcualted using mountain plot script from ViennaRNA (software)
for k in *_dp.ps

	do ./ViennaRNA-2.6.1/src/Utils/mountain.pl $k >> ../Ctrl_outputs/$k.bed

done

# Extracting different scores from the mountain plot as follows: 
for m in *_dp.ps.bed

	do 	
		#echo $m
		num=`echo $(cat $m|wc -l)`
	   	#echo $num
		num_1=$((num-2))
		num_2=$((num_1/3))
		num_left=$((num_2+2))
		num_right_ori=$(expr $num-$num_2)
		num_right=$((num_right_ori-1))
		cat $m|head -"$num_2"|awk '{print $2}' >> $m'_1.bed'
		cat $m|sed -n $num_left,$num_right'p'|awk '{print $2}'>> $m'_2.bed'
		cat $m|tail -"$num_2"|awk '{print $2}' >> $m'_3.bed'
done
# Extracting the maximum values for the three different scores from the mountain plot as follows:
for n in *_1.bed; do cat $n|sort -n|tail -1 >> max_file1.bed; done
for n in *_2.bed; do cat $n|sort -n|tail -1 >> max_file2.bed; done
for n in *_3.bed; do cat $n|sort -n|tail -1 >> max_file3.bed; done