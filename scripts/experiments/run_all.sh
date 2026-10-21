#!/bin/bash
for f in experiments/**/*.yaml; do echo $f; python -m continuumx.cli simulate "$f" || true; done
