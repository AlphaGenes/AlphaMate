
rm -Rf AlphaMate
pullED /exports/cmvm/eddie/eb/groups/hickey_group/ggorjanc/AlphaSuite/AlphaMate.tgz
tar -xzvf AlphaMate.tgz

O="Linux"
V="0.3.0-devel"

# rm -Rf AlphaMate_$V
mkdir -p AlphaMate_$V/{bin,doc,example,util}

cp -f AlphaMate/{NEWS.txt,README.txt} AlphaMate_${V}
cp -fr AlphaMate/bin/AlphaMate AlphaMate_${V}/bin/AlphaMate_${O}_${V}

cp -fr AlphaMate/doc/AlphaMate.pdf AlphaMate_${V}/doc/AlphaMate_${V}.pdf

cp -fr AlphaMate/example/{Animal,PlantCross,PlantSelf,GenomeEditing,CoreDiversitySet} AlphaMate_${V}/example/.
rm -Rf AlphaMate_${V}/example/*/{GenerateData.R,*Summary.txt,Contributors*,MatingPlan*,Optimisation*,SeedUsed.txt,Targets.txt}

cp -fr AlphaMate/util/PlotFrontier.R AlphaMate_${V}/util/.

rm -Rf AlphaMate_${V}/*/.git
rm -Rf AlphaMate_${V}/*/*/.git
rm -Rf AlphaMate_${V}/*/*/*/.git

rm -f AlphaMate_${V}/*/.gitignore
rm -f AlphaMate_${V}/*/*/.gitignore
rm -f AlphaMate_${V}/*/*/*/.gitignore

rm -f AlphaMate_${V}/.DS_Store
rm -f AlphaMate_${V}/*/.DS_Store
rm -f AlphaMate_${V}/*/*/.DS_Store

zip -r9v AlphaMate_${V}.zip AlphaMate_${V}