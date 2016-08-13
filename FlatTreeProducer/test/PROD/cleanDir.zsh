#!/bin/env zsh

ntuplev="MantaRay-patch7-v20150822"

cleanf=1 # clean all ntuples ?
opath="/dpm/in2p3.fr/home/cms/phedex/store/user/kskovpen/FlatTree/${ntuplev}/"

if [[ ${cleanf} == 1 ]]; then
/usr/bin/rfrm -rf ${opath}
/usr/bin/rfmkdir ${opath}
fi

rm -rf crab_*
