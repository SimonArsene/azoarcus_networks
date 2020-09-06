folder=/workbench/my_seq/
prefix=my_seq_chk_
nb_chunks=6
position=151
nb_indexes=3
nb_reads=319434126
python concatenate_reads.py my_seq_R1.fastq my_seq_R2.fastq my_seq.fastq --sample_list=my_seq_sample_list.xlsx -N $nb_reads -distance=hamming
split -d -a 1 -l 212956084 my_seq.fastq my_seq_chk_ --additional-suffix=.fastq
python dpBC_extraction.py $folder $prefix $nb_chunks $position -distance=hamming -N $nb_reads
python align_with_bowtie2.py $folder $prefix $nb_chunks $nb_indexes
python combine_sam.py $folder $prefix $nb_chunks $nb_indexes -N $nb_reads
cat my_seq_chk_*_dpBC.csv > my_seq_dpBC_hamming.csv
python extract_meta_rcd.py $folder $prefix $nb_chunks $nb_indexes $position -distance=hamming -N $nb_reads
cat my_seq_chk_*_meta.csv > my_seq_meta_hamming.csv
