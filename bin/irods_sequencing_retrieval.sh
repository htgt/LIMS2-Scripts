#BASE=/lustre/scratch117/sciops/team87/MiSeq
BASE=/lustre/scratch117/sciops/team87/Mark
DIR="$BASE"/irods_holding_pen/
cd $DIR
HOUR=`date +%H`
DAY=`date +%d`
DATE=`date +%Y-%m-$(($DAY - 1))T$(($HOUR - 1)):%M:%S`

#jq -n '{avus: [{attribute: "study", value: "hiPSC Endo diff v1", o: "="}, {attribute: "target", value: "1", o: "="}], timestamps: [{"created": "2017-06-20T00:00:00", "operator": "n>="}]}' | /software/npg/20170228/bin/baton-metaquery --zone seq --obj | jq -r '.[]' | /software/npg/20170228/bin/baton-get --save --verbose
jq -n '{avus: [{attribute: "id_run", value: "23524", o: "="}, {attribute: "target", value: "1", o: "="}], timestamps: [{"created": "2017-01-01T00:00:00", "operator": "n>="}]}' | /software/npg/20170228/bin/baton-metaquery --zone seq --obj | jq -r '.[]' | /software/npg/20170228/bin/baton-get --save --verbose


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

        /software/solexa/pkg/biobambam/current/bin/bamtofastq inputformat=cram I=$FILE F="$FOLDER/${INDEX}_S${INDEX}_L001_R1.fastq" F2="$FOLDER/${INDEX}_S${INDEX}_L001_R2.fastq"
    fi
done

cd $DIR
rm *
