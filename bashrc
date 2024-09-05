export AIPS_ROOT=/share/data/askap/bqlao/EAVN/software/AIPS
export LD_LIBRARY_PATH=/share/data/askap/bqlao/EAVN/software/lib/lib:$LD_LIBRARY_PATH
export CPATH=/share/data/askap/bqlao/EAVN/software/lib/include:$CPATH

source $AIPS_ROOT/LOGIN.SH

export PATH=/share/data/askap/bqlao/EAVN/software/pgplot5.2:$PATH
export PGPLOT_DIR=/share/data/askap/bqlao/EAVN/software/pgplot5.2
export LD_LIBRARY_PATH=/share/data/askap/bqlao/EAVN/software/pgplot5.2:$LD_LIBRARY_PATH
export PGPLOT_FONT=/share/data/askap/bqlao/EAVN/software/pgplot5.2/grfont.dat
export PGPLOT_DEV=/xwine
export PGPLOT_LIB="-L /usr/X11R6/lib -lX11 -L /share/data/askap/bqlao/EAVN/software/pgplot5.2/ -lpgplot"

export PATH=/share/data/askap/bqlao/EAVN/software/uvf_difmap_2.5q:$PATH

export PATH=/share/data/askap/bqlao/EAVN/software/anaconda3/bin:$PATH
