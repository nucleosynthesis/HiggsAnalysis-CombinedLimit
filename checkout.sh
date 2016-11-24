#!/bin/sh

set -x
set -e

git clone -b combine_tutorial_SWAN https://github.com/nucleosynthesis/HiggsAnalysis-CombinedLimit HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
source env_standalone_make.sh
