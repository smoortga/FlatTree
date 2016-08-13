#!/bin/env zsh

# use blob in the link
wpath="https://github.com/CERN-PH-CMG/cmg-cmssw/blob/CMGTools-from-CMSSW_7_4_12/CMGTools/TTHAnalysis/data/leptonMVA/tth"

wget "${wpath}/el_eta_cb_BDTG.weights.xml?raw=true"
mv el_eta_cb_BDTG.weights.xml\?raw\=true el_eta_cb_BDTG.weights.xml

wget "${wpath}/el_eta_ec_BDTG.weights.xml?raw=true"
mv el_eta_ec_BDTG.weights.xml\?raw\=true el_eta_ec_BDTG.weights.xml

wget "${wpath}/el_eta_fb_BDTG.weights.xml?raw=true"
mv el_eta_fb_BDTG.weights.xml\?raw\=true el_eta_fb_BDTG.weights.xml

wget "${wpath}/mu_eta_b_BDTG.weights.xml?raw=true"
mv mu_eta_b_BDTG.weights.xml\?raw\=true mu_eta_b_BDTG.weights.xml

wget "${wpath}/mu_eta_e_BDTG.weights.xml?raw=true"
mv mu_eta_e_BDTG.weights.xml\?raw\=true mu_eta_e_BDTG.weights.xml
