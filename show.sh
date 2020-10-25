#!/bin/sh

dot -Tpdf "$1" > temp.pdf && evince temp.pdf
