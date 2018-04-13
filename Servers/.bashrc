# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

# Uncomment the following line if you don't like systemctl's auto-paging feature:
# export SYSTEMD_PAGER=

# User specific aliases and functions

#export PYTHONPATH=$PYTHONPATH:/opt/anaconda2/lib/python2.7/site-packages/

export LD_LIBRARY_PATH=/opt/anaconda2/lib:$LD_LIBRARY_PATH

source /usr/local/bin/thisroot.sh

export PATH=/opt/anaconda2/bin/:$PATH

#alias python=/opt/anaconda2/bin/python
#alias f2py=/opt/anaconda2/bin/f2py
