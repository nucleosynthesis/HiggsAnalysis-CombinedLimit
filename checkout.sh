#!/bin/sh

set -x
set -e

#git clone -b combine_tutorial_SWAN https://github.com/nucleosynthesis/HiggsAnalysis-CombinedLimit HiggsAnalysis/CombinedLimit
mkdir -p HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
git init 
git remote add origin https://github.com/nucleosynthesis/HiggsAnalysis-CombinedLimit
git pull origin combine_tutorial_SWAN
