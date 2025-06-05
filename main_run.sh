#!/bin/bash
# Usage: bash this_script.sh target.sdf database.sdf flexible/rigid metric(tanimoto/carbo) chargemethod(gasteiger/mmff/espaloma) sortmethod(ESP*Shape/ESP/Shape) extract(yes/no) numbertoextract
if [ "$3" = "flexible" ]; then
obabel $2 -O query.mol2 -m
obabel $1 -O target.mol2
for j in query*.mol2; do ./LSalign $j target.mol2  -rf 1 -md 1 -acc 1 -H 1 -o output.txt;grep QUE output.txt > `head -n 2 $j | tail -n 1`.pdb; obabel -ipdb `head -n 2 $j| tail -n 1`.pdb -O out.sdf 2>&1; rm $j `head -n 2 $j | tail -n 1`.pdb ;cat out.sdf >> aligned.sdf; done;	
cat $1 aligned.sdf > database.sdf;
python no_align.py --target database.sdf --metric $4 --chargemethod $5 --sort $6 --sanitize False;
fi

if [ "$3" = "rigid" ]; then
obabel $2 -O query.mol2 -m
obabel $1 -O target.mol2
for j in query*.mol2; do ./LSalign $j target.mol2  -rf 0 -H 1 -o output.txt;grep QUE output.txt > `head -n 2 $j | tail -n 1`.pdb; obabel -ipdb `head -n 2 $j| tail -n 1`.pdb -O out.sdf 2>&1; rm $j `head -n 2 $j | tail -n 1`.pdb ;cat out.sdf >> aligned.sdf; done
cat $1 aligned.sdf > database.sdf;
python no_align.py --target database.sdf --metric $4 --chargemethod $5 --sort $6 --sanitize False;
fi

if [ "$3" = "rdkit" ];then
	cat $1 $2 > database.sdf;
python option_5.py --target database.sdf --metric $4 --chargemethod $5 --sort $6 --sanitize False;
fi


if [ "$7" = "yes" ]; then
	a=$8;b=$(($a + 1));for j in `head -n $b sorted.csv | sed '1d' | sort -u |awk -F , '{print $1}'` ; do echo -n $j" " ;done|  rev | cut -c2- | rev >> list.txt; for j in `cat list.txt`; do sed -ne "/$j\\b/,/\$\$\$\$/{/\$$$$/!p;/\$\$\$\$/q}" $2 >> Extracted.sdf
done;fi
