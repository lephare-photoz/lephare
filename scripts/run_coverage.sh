#!/usr/bin/env bash
set -euo pipefail

MODE="${1:-local}"   # local | ci
ROOT_DIR="$(pwd)"

echo "===== Coverage mode: ${MODE} ====="

echo "===== Cleaning old artifacts ====="
rm -f coverage*.info coverage.txt coverage.lcov
rm -rf cov_html
rm -rf build

echo "===== Building with C++ coverage enabled ====="
export INCLUDE_COVERAGE=1
LEPHARE_ENABLE_AVX2=0 && \
    LEPHARE_NATIVE_OPT=0 &&\
    pip install -qq .[dev]

echo "===== Running Python tests (Python coverage) ====="
pytest -s tests \
       --cov=lephare \
       --cov-report=lcov:coverage.lcov

#find build -name '_bindings.cc.gcda' -delete

echo "===== Running geninfo ====="
geninfo --filter range \
	--ignore-errors gcov,gcov \
	--ignore-errors mismatch \
	--debug build \
	--keep-going

echo "===== Capturing C++ coverage ====="
lcov --filter range \
     --ignore-errors gcov,gcov \
     --output-file coverage.cpp \
     --capture \
     --directory build  \
     --exclude '*_bindings.cc'

lcov --output-file coverage.cpp \
     --extract coverage.cpp $PWD/src/"*"


echo "===== Merging Python + C++ coverage ====="
cat coverage.lcov coverage.cpp > coverage.txt

# lcov --add-tracefile coverage.lcov \
#      --add-tracefile coverage.cpp \
#      --output-file coverage.txt

echo "===== Coverage summary ====="
lcov --summary coverage.txt

if [[ "${MODE}" == "local" ]]; then
  echo "===== Generating HTML report ====="
  genhtml coverage.txt -o cov_html
  echo "HTML report in ./cov_html/index.html"
else
  echo "===== CI mode: skipping HTML generation ====="
fi

echo "===== Coverage complete ====="
