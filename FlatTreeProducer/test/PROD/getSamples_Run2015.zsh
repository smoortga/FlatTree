#!/bin/env zsh

tag=Run2016B-PromptReco-v2

wget --no-check-certificate \
--output-document=samples_Run2015.txt \
"https://cmsweb.cern.ch/das/request?view=plain&instance=prod%2Fglobal&input=dataset%3D%2F*%2F*${tag}*%2FMINIAOD+|+sort+dataset.name"
