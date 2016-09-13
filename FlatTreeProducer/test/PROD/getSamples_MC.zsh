#!/bin/env zsh

das_client.py --limit=0 --query="dataset=/*/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2*/MINIAODSIM" | sort \
> samples_MC.txt
