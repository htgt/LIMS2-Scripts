DATE=`date +%Y-%m-%dT%H:%M:%S`
echo $DATE
#jq -n '{avus: [{attribute: "study", value: "iPS MiSEQ Genotyping", o: "="}, {attribute: "target", value: "1", o: "="}], timestamps: [{"created": "2014-01-01T00:00:00", "operator": "n>="}]}' | /software/npg/20170228/bin/baton-metaquery --zone seq --obj  | jq -r '.[]' | /software/npg/20170228/bin/baton-get --save --verbose
