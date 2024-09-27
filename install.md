# [AIPS](http://www.aips.nrao.edu/dec24.shtml)
- Download AIPS wizard file 
```
wget http://www.aips.nrao.edu/31DEC24/install.pl
```
- Assign Execution Permissions
```
chmod +x install.pl
```
- install AIPS 
```
perl install.pl -n
```
- AIPS setting
```
   AIPS_ROOT (screen 3): /share/data/askap/bqlao/EAVN/software/AIPS
       Group (screen 4): polarization
 Group Write (screen 4): YES
Architecture (screen 5): LNX64
   Site name (screen 5): BQLAO
  AIPS hosts (screen 6): LOCALHOST
  Data areas (screen 7): /share/data/askap/bqlao/EAVN/software/AIPS/DATA/LOCALHOST_1
    Printers (screen 8):   Paper type (screen 8): A4
 Tape drives (screen 9): 
  Tape hosts (screen 9): 127.0.0.1
```
- add follows content to /etc/services
```
sssin           5000/tcp        SSSIN      # AIPS TV server
ssslock         5002/tcp        SSSLOCK    # AIPS TV Lock
msgserv         5008/tcp        MSGSERV    # AIPS Message Server
tekserv         5009/tcp        TEKSERV    # AIPS TekServer
aipsmt0         5010/tcp        AIPSMT0    # AIPS remote FITS disk access
aipsmt1         5011/tcp        AIPSMT1    # AIPS remote tape 1
aipsmt2         5012/tcp        AIPSMT2    # AIPS remote tape 2
aipsmt3         5013/tcp        AIPSMT3
aipsmt4         5014/tcp        AIPSMT4
aipsmt5         5015/tcp        AIPSMT5
aipsmt6         5016/tcp        AIPSMT6
aipsmt7         5017/tcp        AIPSMT7

```
- disable the original item 5002 
```
#rfe 5002/udp
#rfe 5002/tcp
```
- add AIPS environment and test
```
>source /share/data/askap/bqlao/EAVN/software/AIPS/LOGIN.SH 
>aips
```
# PGPLOT
- download
```
wget ftp://ftp.astro.caltech.edu/pub/pgplot/pgplot5.2.tar.gz
```
- install
```
../pgplot/makemake ../pgplot linux g77_gcc shared
```

# difmap
- download
```
wget ftp://ftp.astro.caltech.edu/pub/difmap/difmap2.5q.tar.gz
```
- pre-install
```
export PGPLOT_LIB="-L/pgplot5.2_lib_dir -L/X11_lib_dir \

    -Xlinker -R/pgplot5.2_lib_dir:X11_lib_dir -lpgplot -lX11"
```
- install
```
./configure linux-i486-gcc
./makeall
```

# Anaconda (python >=3.8,<3.9.0a0)
- download and install
```
wget https://mirrors.tuna.tsinghua.edu.cn/anaconda/archive/Anaconda3-2021.05-Linux-x86_64.sh
```
or  
```
wget --user-agent="Mozilla" https://mirrors.tuna.tsinghua.edu.cn/anaconda/archive/Anaconda3-2021.05-Linux-x86_64.sh
```
- Assign Execution Permissions  
```
chmod +x Anaconda3-2021.05-Linux-x86_64.sh
```
- install  
```
./Anaconda3-2021.05-Linux-x86_64.sh
```

# [ParselTongue](https://www.jive.eu/jivewiki/doku.php?id=parseltongue:parseltongue) (3.0)
```
conda install -c kettenis parseltongue 
```

# Set environment
## example
```
#aips
export AIPS_ROOT=/media/hero/Intel6/EAVN/software/AIPS
source $AIPS_ROOT/LOGIN.SH

#pgplot
export PATH=/media/hero/Intel6/mwa_software/pgplot:$PATH
export PGPLOT_DIR=/media/hero/Intel6/mwa_software/pgplot
export LD_LIBRARY_PATH=/media/hero/Intel6/mwa_software/pgplot:$LD_LIBRARY_PATH
export PGPLOT_FONT=/media/hero/Intel6/mwa_software/pgplot/grfont.dat
export PGPLOT_DEV=/xwine
export PGPLOT_LIB="-L /usr/lib/x86_64-linux-gnu -lX11 -L /media/hero/Intel6/mwa_software/pgplot -lpgplot"

export PATH=/media/hero/Intel6/EAVN/software/uvf_difmap_2.5q:$PATH
export LD_LIBRARY_PATH=/media/hero/Intel6/EAVN/software/uvf_difmap_2.5q/lib:$LD_LIBRARY_PATH
export CPATH=/media/hero/Intel6/EAVN/software/uvf_difmap_2.5q/inlcude:$CPATH

#anaconda3
export PATH=/media/hero/Intel6/EAVN/software/anaconda3/bin:$PATH

```
