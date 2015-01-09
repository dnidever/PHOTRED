#!/bin/csh
set input=${1}
daogrow << DONE
photo.opt

${input}.inf
${input}.ext
3
0.9,0.0
0.2
DONE
