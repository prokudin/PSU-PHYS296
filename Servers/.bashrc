# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

# Uncomment the following line if you don't like systemctl's auto-paging feature:
# export SYSTEMD_PAGER=

# User specific aliases and functions

export PYTHONPATH=$PYTHONPATH:/opt/anaconda2/lib

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/anaconda2/lib

export PATH=$PATH:/opt/anaconda2/bin

alias python=/opt/anaconda2/bin/python
alias f2py=/opt/anaconda2/bin/f2py
