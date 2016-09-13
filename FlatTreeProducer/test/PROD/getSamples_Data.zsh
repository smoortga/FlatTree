#!/bin/env zsh

das_client.py --limit=0 --query="dataset=/*/Run2016*PromptReco-v2/MINIAOD" | sort \
> samples_Data.txt
