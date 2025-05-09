# Run this in your terminal
cd /home/atoffano/PFP_baselines/
for db in 2023_01 2022_01 2021_01 2020_01 2019_01 2018_01 2017_01 2016_01 2015_01 2014_01 2013_01 2012_01 2011_01 2010_01 2009_03 2008_01 2007_03 2006_02 2005_01 2003_12
do
    for aspect in bp mf cc
    do
        old="${db}/swissprot_${db}_${aspect}_annotations.tsv"
        if [ -f "$old" ]; then
            case $aspect in
                bp) new="${db}/swissprot_${db}_BPO_annotations.tsv" ;;
                mf) new="${db}/swissprot_${db}_MFO_annotations.tsv" ;;
                cc) new="${db}/swissprot_${db}_CCO_annotations.tsv" ;;
            esac
            mv "$old" "$new"
        fi
    done
done