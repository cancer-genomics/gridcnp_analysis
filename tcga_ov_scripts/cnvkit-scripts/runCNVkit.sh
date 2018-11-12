#!/bin/bash

#$ -cwd
#$ -j y
#$ -l mem_free=4G
#$ -l h_vmem=4G
#$ -l h_rt=96:00:00
#$ -l h_fsize=100G
#$ -l gwas

module load python/2.7.9

#-------
# Input
#-----------------------------------------------------------------------------------------------------------------------------------------
targetBED=/dcl01/scharpf1/data/dbruhm/elsa/melanoma/CNV/dataFiles/hg19SureSelect50AllExons_Regions.noalts.tab.bed
normalBams=/dcl01/scharpf1/data/dbruhm/elsa/melanoma/CNV/bamLocations/normalBams.txt
tumorBams=/dcl01/scharpf1/data/dbruhm/elsa/melanoma/CNV/bamLocations/tumorBams.txt
outDir=/dcl01/scharpf1/data/dbruhm/elsa/melanoma/CNV/cnvkit_output
refFlat=/dcl01/scharpf1/data/dbruhm/elsa/melanoma/CNV/dataFiles/refFlat.txt
refGenome=/dcl01/scharpf/data/reference/no_alternates/hg19/hg19.fa
access=/dcl01/scharpf1/data/dbruhm/elsa/melanoma/CNV/dataFiles/access-5k-mappable.hg19.bed
cnvkitDir=/users/dbruhm/Tools/cnvkit
#-----------------------------------------------------------------------------------------------------------------------------------------


#---------------------------------
# Setting up directory structure
#-----------------------------------------------
mkdir -p $outDir/target/
mkdir -p $outDir/coverage/
mkdir -p $outDir/coverage/normal/
mkdir -p $outDir/coverage/tumor/
mkdir -p $outDir/reference/
mkdir -p $outDir/tumor-cn/
mkdir -p $outDir/summary/heatmap-tumors
mkdir -p $outDir/status
mkdir -p $outDir/status/single
mkdir -p $outDir/status/tumors
mkdir -p $outDir/status/normals
mkdir -p $outDir/leave1out
#-----------------------------------------------


#------------------------------------------
# Generate target and antitarget bin files
#-----------------------------------------------------------------
python $cnvkitDir/cnvkit.py target $targetBED \
    --split \
    --annotate $refFlat \
    --short-names \
    -o $outDir/target/target.bed

python $cnvkitDir/cnvkit.py antitarget $targetBED \
    --access $access \
    -o $outDir/target/antitarget.bed
#-----------------------------------------------------------------


#----------------
# Count normals
#-----------------------------------------------------------------------
sn=countNormals.sh

echo -e '#!/bin/bash\n' >> $sn
echo '#$ -cwd' >> $sn
echo '#$ -j y' >> $sn
echo '#$ -l mem_free=8G' >> $sn
echo '#$ -l h_vmem=8G' >> $sn
echo '#$ -l h_rt=24:00:00' >> $sn
echo '#$ -l h_fsize=100G' >> $sn
echo -n '#$ -t 1-' >> $sn
echo $(cat $normalBams | wc -l) >> $sn
echo >> $sn
echo 'module load python/2.7.9' >> $sn
echo >> $sn
echo -n 'bamFile=$(cat ' >> $sn
echo -n $normalBams >> $sn
echo ' | head -n $SGE_TASK_ID | tail -n 1)' >> $sn
echo 'sampleID=$(basename $bamFile | cut -d "." -f1)' >> $sn

cat <<EOF >> $sn

if [ ! -f $outDir/coverage/normal/\$sampleID.targetcoverage.cnn ]; then
  python $cnvkitDir/cnvkit.py coverage \\
    \$bamFile \\
    $outDir/target/target.bed \\
    -o $outDir/coverage/normal/\$sampleID.targetcoverage.cnn
fi

if [ ! -f $outDir/coverage/normal/\$sampleID.antitargetcoverage.cnn ]; then
  python $cnvkitDir/cnvkit.py coverage \\
    \$bamFile \\
    $outDir/target/antitarget.bed \\
    -o $outDir/coverage/normal/\$sampleID.antitargetcoverage.cnn
fi

touch $outDir/status/normals/countNormals.\$sampleID.ck
EOF

qsub countNormals.sh
#-----------------------------------------------------------------------

#----------------
# Count tumors
#-----------------------------------------------------------------------
sn=countTumors.sh

echo -e '#!/bin/bash\n' >> $sn
echo '#$ -cwd' >> $sn
echo '#$ -j y' >> $sn
echo '#$ -l mem_free=8G' >> $sn
echo '#$ -l h_vmem=8G' >> $sn
echo '#$ -l h_rt=24:00:00' >> $sn
echo '#$ -l h_fsize=100G' >> $sn
echo -n '#$ -t 1-' >> $sn
echo $(cat $tumorBams | wc -l) >> $sn
echo >> $sn
echo 'module load python/2.7.9' >> $sn
echo >> $sn
echo -n 'bamFile=$(cat ' >> $sn
echo -n $tumorBams >> $sn
echo ' | head -n $SGE_TASK_ID | tail -n 1)' >> $sn
echo 'sampleID=$(basename $bamFile | cut -d "." -f1)' >> $sn

cat <<EOF >> $sn

if [ ! -f $outDir/coverage/tumor/\$sampleID.targetcoverage.cnn ]; then
  python $cnvkitDir/cnvkit.py coverage \\
    \$bamFile \\
    $outDir/target/target.bed \\
    -o $outDir/coverage/tumor/\$sampleID.targetcoverage.cnn
fi

if [ ! -f $outDir/coverage/tumor/\$sampleID.antitargetcoverage.cnn ]; then
  python $cnvkitDir/cnvkit.py coverage \\
    \$bamFile \\
    $outDir/target/antitarget.bed \\
    -o $outDir/coverage/tumor/\$sampleID.antitargetcoverage.cnn
fi

touch $outDir/status/tumors/countTumors.\$sampleID.ck
EOF

qsub countTumors.sh
#-----------------------------------------------------------------------


while [ $(ls -1v $outDir/status/normals | wc -l) -lt $(cat $normalBams | wc -l) ]; do
  sleep 10
  echo -n '.'
done

rm $outDir/status/normals/countNormals*

#-----------------------------------------------
# Create a reference from the panel of normals
#-------------------------------------------------------------------------------
python $cnvkitDir/cnvkit.py reference \
    --output $outDir/reference/reference.cnn \
    --fasta $refGenome \
    $outDir/coverage/normal/*.targetcoverage.cnn \
    $outDir/coverage/normal/*.antitargetcoverage.cnn 2> \
    $outDir/reference/reference.cnn.log
#-------------------------------------------------------------------------------


while [ $(ls -1v $outDir/status/tumors | wc -l) -lt $(cat $tumorBams | wc -l) ]; do
  sleep 10
  echo -n '.'
done

rm $outDir/status/tumors/countTumors*


#-------------------------------------------------------------
# Normalize each tumor sample against the panel of normals
#-------------------------------------------------------------------
while read path
do
sampleID=$(basename $path | cut -d '.' -f1)
mkdir $outDir/tumor-cn/$sampleID
python $cnvkitDir/cnvkit.py fix \
     $outDir/coverage/tumor/$sampleID.targetcoverage.cnn \
     $outDir/coverage/tumor/$sampleID.antitargetcoverage.cnn \
     $outDir/reference/reference.cnn \
    --output  $outDir/tumor-cn/$sampleID/$sampleID.cnr 2> \
    $outDir/tumor-cn/$sampleID/$sampleID.cnr.log
done < $tumorBams
#--------------------------------------------------------------------


#-------------------------------------------------------------
# Segment the bin log ratios for each tumor sample
#--------------------------------------------------------------------
while read path
do
sampleID=$(basename $path | cut -d '.' -f1)
python $cnvkitDir/cnvkit.py segment \
    $outDir/tumor-cn/$sampleID/$sampleID.cnr \
    --output $outDir/tumor-cn/$sampleID/$sampleID.cns 2> \
    $outDir/tumor-cn/$sampleID/$sampleID.cns.log
done < $tumorBams
#--------------------------------------------------------------------


#-------------------------------------------------------
# Generate genome-wide segmented scatter for tumors
#--------------------------------------------------------------------
while read path
do
sampleID=$(basename $path | cut -d '.' -f1)
python $cnvkitDir/cnvkit.py scatter \
    $outDir/tumor-cn/$sampleID/$sampleID.cnr \
    -s $outDir/tumor-cn/$sampleID/$sampleID.cns \
    -i $sampleID \
    --output $outDir/tumor-cn/$sampleID/$sampleID.scatter.pdf
done < $tumorBams
#--------------------------------------------------------------------


#--------------------------------------------------------------
# Make cohort-level heatmap of CN segments from all tumors
#-------------------------------------------------------------------------------------
while read path
do
sampleID=$(basename $path | cut -d '.' -f1)
ln -s $outDir/tumor-cn/$sampleID/$sampleID.cns $outDir/summary/heatmap-tumors/
done < $tumorBams

python /users/dbruhm/Tools/cnvkit/cnvkit.py heatmap \
    $outDir/summary/heatmap-tumors/*.cns \
    --output $outDir/summary/heatmap-tumors/tumors.heatmap.pdf
#-------------------------------------------------------------------------------------


#--------------
# Leave 1 out
#---------------------------------------------------------------------------------------------------------------
while read path
do
sampleID=$(basename $path | cut -d '.' -f1)
mkdir -p $outDir/leave1out/$sampleID
mkdir -p $outDir/leave1out/$sampleID/coverage
mkdir -p $outDir/leave1out/$sampleID/coverage/normal
mkdir -p $outDir/leave1out/$sampleID/coverage/tumor
mkdir -p $outDir/leave1out/$sampleID/reference
mkdir -p $outDir/leave1out/$sampleID/target
mkdir -p $outDir/leave1out/$sampleID/normal-cn/
mkdir -p $outDir/leave1out/summary/heatmap-normals

ln -s $outDir/target/* $outDir/leave1out/$sampleID/target
 
ln -s $outDir/coverage/normal/* $outDir/leave1out/$sampleID/coverage/normal
unlink $outDir/leave1out/$sampleID/coverage/normal/$sampleID.targetcoverage.cnn
unlink $outDir/leave1out/$sampleID/coverage/normal/$sampleID.antitargetcoverage.cnn
ln -s $outDir/coverage/normal/$sampleID.targetcoverage.cnn $outDir/leave1out/$sampleID/coverage/tumor
ln -s $outDir/coverage/normal/$sampleID.antitargetcoverage.cnn $outDir/leave1out/$sampleID/coverage/tumor

python $cnvkitDir/cnvkit.py reference \
    --output $outDir/leave1out/$sampleID/reference/reference.cnn \
    --fasta $refGenome \
    $outDir/leave1out/$sampleID/coverage/normal/*.targetcoverage.cnn \
    $outDir/leave1out/$sampleID/coverage/normal/*.antitargetcoverage.cnn 2> \
    $outDir/leave1out/$sampleID/reference/reference.cnn.log

python $cnvkitDir/cnvkit.py fix \
    $outDir/leave1out/$sampleID/coverage/tumor/$sampleID.targetcoverage.cnn \
    $outDir/leave1out/$sampleID/coverage/tumor/$sampleID.antitargetcoverage.cnn \
    $outDir/leave1out/$sampleID/reference/reference.cnn \
    --output  $outDir/leave1out/$sampleID/normal-cn/$sampleID.cnr 2> \
    $outDir/leave1out/$sampleID/normal-cn/$sampleID.cnr.log

python $cnvkitDir/cnvkit.py segment \
    $outDir/leave1out/$sampleID/normal-cn/$sampleID.cnr \
    --output $outDir/leave1out/$sampleID/normal-cn/$sampleID.cns 2> \
    $outDir/leave1out/$sampleID/normal-cn/$sampleID.cns.log

python $cnvkitDir/cnvkit.py scatter \
    $outDir/leave1out/$sampleID/normal-cn/$sampleID.cnr \
    -s $outDir/leave1out/$sampleID/normal-cn/$sampleID.cns \
    -i $sampleID \
    --output $outDir/leave1out/$sampleID/normal-cn/$sampleID.scatter.pdf

ln -s $outDir/leave1out/$sampleID/normal-cn/$sampleID.cns $outDir/leave1out/summary/heatmap-normals
done < $normalBams

python /users/dbruhm/Tools/cnvkit/cnvkit.py heatmap \
    $outDir/leave1out/summary/heatmap-normals/*.cns \
    --output $outDir/leave1out/summary/heatmap-normals/normals.heatmap.pdf

echo "Done."
#---------------------------------------------------------------------------------------------------------------


