#!/bin/bash

$PYTHON setup.py clean
$PYTHON setup.py build_ext
$PYTHON setup.py install --single-version-externally-managed --record=record.txt