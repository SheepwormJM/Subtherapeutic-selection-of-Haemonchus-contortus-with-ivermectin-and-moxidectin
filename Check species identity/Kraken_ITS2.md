Obtained script and custom made Kraken ITS2 database from Benedict Karani, University of Glasgow, August 2025.  

```
#!/bin/sh 
#PBS -l walltime=50:00:00 
#PBS -l cput=20:00:00 
#PBS -l nodes=1:centos7
#PBS -d /export/home2/jmi45g/lainglab/HCON/Selection_lines_poolseq/Fastq_files/renamed_fastq_files

export PATH="/export/home2/jmi45g/anaconda3/envs/Kraken2/bin:$PATH"

for i in $(ls *_R1.fastq.gz | sed "s/_R1.fastq.gz//"); do
  R1_fastq="${i}_R1.fastq.gz"
  R2_fastq="${i}_R2.fastq.gz"

  if [[ -f "$R2_fastq" ]]; then
    echo "Processing $i - $R1_fastq and $R2_fastq" | tee -a kraken.log
    kraken2 --db /export/home2/jmi45g/lainglab/Benedicts_kraken2_helminths_ITS2_db \
      --output "${i}.kraken2" \
      --report "${i}.ITS2.k2report" \
      --threads 16 \
      --paired \
      "$R1_fastq" "$R2_fastq" >> kraken.log 
  else                                          
    echo "Warning: R2 file not found for $i" | tee -a kraken.log
  fi                                                            
done
```
