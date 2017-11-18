#!/usr/bin/env python
import os,sys
import subprocess
import time

def gen_xml(conf):

  xml=[]
  xml.append('<Request>')
  xml.append('<Project name="solid"/>')
  xml.append('<Email email="nsato@jlab.org" request="true" job="false"/>')
  xml.append('<Track name="analysis"/>')
  #xml.append('<Track name="debug"/>')
  xml.append('<OS name="centos7"/>')
  xml.append('<CPU core="%d"/>'%conf['cpu'])
  xml.append('<Name name="%s"/>'%conf['name'])
  xml.append('<TimeLimit time ="%d" unit="hours"/>'%conf['time'])
  xml.append('<Job>')
  xml.append('<Command>')
  xml.append('<![CDATA[')
  xml.append('cd %s'%conf['path'])
  xml.append('source setup.csh')
  xml.append('cd fitlab')
  xml.append('%s'%conf['cmd'])
  xml.append(']]>')
  xml.append('</Command>')
  xml.append('<Memory space="6050" unit="MB"/>')
  xml.append('<Stdout dest="%s/out/%s.out"/>'%(conf['path'],conf['name']))
  xml.append('<Stderr dest="%s/err/%s.err"/>'%(conf['path'],conf['name']))
  xml.append('</Job>')
  xml.append('</Request>')
  xml=[l+'\n' for l in xml]
  return xml

def send_job(conf):
  xml=gen_xml(conf)
  F=open('job.xml','w')
  F.writelines(xml)
  F.close()

  while 1:
    cmd=['jsub','-xml','job.xml']
    out=subprocess.check_output(cmd)
    if ('error' in out)==False:
      print 'cmd %s has been submitted'%conf['name']
      break
    else:
      delay=5
      print 'error occured, waiting for %d seconds for resubmission'%delay
      time.sleep(delay)

if __name__=='__main__':

  # itrans == 0: Nab, pions only
  # itrans == 1: Nab, pions+kaons
  itrans=0

  conf={}
  conf['path']='/lustre/volatile/JAM/nsato/tmdpheno/fitpack'
  conf['cpu']  = 1
  conf['time'] = 24

  for irun in range(10):

    conf['name'] = 'trans%d-%d'%(itrans,irun)
    conf['cmd']  = './resman.py -t 3 -i %d inputs/transversity/transversity-%d.py'%(irun,itrans)
    send_job(conf)

  #for irun in range(10):

  #  conf['name'] = 'transl%d-%d'%(itrans,irun)
  #  conf['cmd']  = './resman.py -t 3 -i %d inputs/transversity/transversity-lattice-%d.py'%(irun,itrans)
  #  send_job(conf)

