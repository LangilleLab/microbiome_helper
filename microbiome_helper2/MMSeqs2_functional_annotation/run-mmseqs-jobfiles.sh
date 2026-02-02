DB_PATH="/home/dhwani/databases/mmseqsRefSeqCompleteDB"
echo "Creating Index .... mmseqs createindex $DB_PATH tmp"
#mmseqs createindex $DB_PATH tmp
echo "loading DB index into memory ... mmseqs touchdb targetDB"
mmseqs touchdb $DB_PATH

for jobfile in `cat $PWD/jobfiles.list`
do
    #bowtie2-build $SET-7500.fa $SET-7500
    cat $jobfile
    bash $jobfile

done
