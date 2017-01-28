# FlatTree

Code to produce Flat Tree from MINIAOD

```c++
git cms-init

git clone https://github.com/kskovpen/FlatTree

git cms-merge-topic shervin86:Moriond2017_JEC_energyScales
cd EgammaAnalysis/ElectronTools/data; git clone git@github.com:ECALELFS/ScalesSmearings.git; cd -
git cms-merge-topic ikrav:egm_id_80X_v2

git cms-merge-topic -u cms-met:fromCMSSW_8_0_20_postICHEPfilter

scram b -j10
```
