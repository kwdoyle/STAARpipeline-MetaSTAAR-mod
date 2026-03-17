#!/bin/bash

basedir=~/STAARpipeline-MetaSTAAR-mod/
ts -S 5

for i in {1..22}; do
	ts ${basedir}/aux_scripts/anno_variants.sh $i
done

