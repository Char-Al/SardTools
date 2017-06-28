# SardTools

Tools for SARDINe

--------------------------------------------------------------------------------

## make sources files

### NAME

makeSourcesFiles.pl

### VERSION

version 0.0.1b

### SYNOPSIS

makeSourcesFiles.pl -b input.bed -c comment \[-o output.vcf\]

Get a bed of amplicons and transform it on sources files for SARDINe

### DESCRIPTION

This script transform a bed of sources files for SARDINe

### ARGUMENTS

#### General

    -h,--help       Print this help

    -m,--man        Open man page

#### Mandatory arguments

    -b,--bed        input.bed     Specify the bed on input

    -c,--comment    comment       String for the comment column (define version)

#### Optional arguments

    -o,--output     /path/output/repertory    You can specify an output repertory (default: ./)

    -e,--example    launch example (equivalent to: 'perl makeSourcesFiles.pl -v test/test.bed -c test/')

### AUTHORS

- Charles VAN GOETHEM

### TODO

- Add option to change genome assembly

- improved speed (algorithm really not optimized)

- clean and comment code

--------------------------------------------------------------------------------

## make sources files

### NAME

coverageByGene.R

### Usage

Rscript --vanilla coverageByGene.R sample_DPamplicons.csv /path/to/output/

--------------------------------------------------------------------------------

### LICENCE

MIT License

Copyright (c) 2017

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
