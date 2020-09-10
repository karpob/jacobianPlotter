# jacobianPlotter
you need to make sure to checkout graphics submodule, do this by:
git clone --recursve ...
or after clone
git submodule update --init --recursive
To run it for stuff that comes directly out of kmatrix:
./crisFsrSwirCrtm.py  --path $PWD/data --outpath $PWD/plots
To run it for stuff that goes into the GSI solver:
./crisFsrSwirGsi.py  --path $PWD/data --outpath $PWD/plots


Basic plotter to generate plots of jacobians from ascii files spit out of the GSI. both the jacobians straight out of the CRTM, and ones modified by the GSI to put in terms of control variables, along with some kind of inflation/deflation (especially water vapor above the tropopause). 
