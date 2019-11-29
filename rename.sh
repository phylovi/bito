sed -i s/pcss/pcsp/g src/*
git mv doc/pcss.svg doc/pcsp.svg
sed -i s/pcss/pcsp/ doc/pcsp.svg
sed -i s/PCSS/PCSP/g src/*
sed -i s/PCSS/PCSP/g doc/*
sed -i s/PCSS/PCSP/g doc/tex/main.tex
sed -i s/PCSS/PCSP/g test/*py
sed -i s/PCSS/PCSP/g doc/concepts.rst
