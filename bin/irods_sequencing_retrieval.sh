BASE=/lustre/scratch117/sciops/team87/MiSeq
DIR="$BASE"/irods_holding_pen/
cd $DIR
HOUR=`date +%H`
YEAR=`date +%Y` #Testing
DATE=`date +$(($YEAR - 1))-%m-%dT$(($HOUR - 1)):%M:%S`

jq -n '{avus: [{attribute: "study", value: "iPS MiSEQ Genotyping BGE", o: "="}, {attribute: "target", value: "1", o: "="}], timestamps: [{"created": "'"$DATE"'", "operator": "n>="}]}' | /software/npg/20170228/bin/baton-metaquery --zone seq --obj | jq -r '.[]' | /software/npg/20170228/bin/baton-get --save --verbose

echo "Files retrieved from irods"

REGEX="\/([0-9]*)_1#([0-9]*)."
for FILE in "$DIR"/*
do
    if [[ $FILE =~ $REGEX ]]
    then
        RUN="${BASH_REMATCH[1]}"
        FOLDER="${BASE}/run_${RUN}"

        if [ ! -d "$FOLDER" ]
        then
            cd $BASE
            mkdir "run_${RUN}"
            echo "Created"
        fi
        INDEX="${BASH_REMATCH[2]}"
        echo "Converting $INDEX"

        /software/solexa/pkg/biobambam/current/bin/bamtofastq inputformat=cram exclude=SECONDARY,SUPPLEMENTARY,QCFAIL I=$FILE F="$FOLDER/${INDEX}_S${INDEX}_L001_R1.fastq" F2="$FOLDER/${INDEX}_S${INDEX}_L001_R2.fastq"
    fi
done

cd $DIR
rm *
