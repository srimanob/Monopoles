#!/bin/env python
'''python script to plot ecalMaps with monopole and ecay overlays'''


from ROOT import *

#from pyrootlogon import *



data_file=""
sub_directory="demo"
event_number = int(1)
output="null"


def ecalMapPlotter(fileName,wevent=1,subPath='demo',output="null"):
  '''ecalMapPlotter'''

  if fileName == None:
    print('Must give fileName')

  print 'Welcome to ecalMapPlotter'

  file = TFile(fileName)
  tDir = file.Get(subPath)

  tree = tDir.Get('tree')

  count = 0

  run = 0
  lumi = 0
  eventN = 0

  m1eta = 1e3
  m1phi = 1e3
  m2eta = 1e3
  m2phi = 1e3

  m1deta = []
  m1dphi = []
  m2deta = []
  m2dphi = []

  print 'looping over tree to find monopole information'


  for event in tree:
    count += 1
    
    if count == wevent:
      run = event.run
      lumi = event.lumi
      eventN = event.event

      m1eta = event.monoExp_eta
      m1phi = event.monoExp_phi
      m2eta = event.antiExp_eta
      m2phi = event.antiExp_phi

      m1deta = event.md_eta
      m1dphi = event.md_phi
      m2deta = event.ad_eta
      m2dphi = event.ad_phi

      break


  


  eMapName = 'eMap_%d_%d_%d' % (run,lumi,eventN)
  print 'retrieving ECal Map',eMapName
  eMapTemp = tDir.Get(eMapName)

  eMap = TH2D('eMap',eMapName,1,0,10,1,0,10)
  eMapTemp.Copy(eMap)

  

  can = TCanvas('c1','',1200,800)
  gStyle.SetOptStat(0) 
  eMap.Draw('colz')
  can.SetLogz()


  l = [TLatex()]
  l.append(l[0].DrawLatex(m1eta,m1phi,'M'))
  l.append(l[0].DrawLatex(m2eta,m2phi,'A'))

  for i in xrange(0,len(m1deta)):
    l.append(l[0].DrawLatex( m1deta[i],m1dphi[i],'mp' ) )


  for i in xrange(0,len(m2deta)):
    l.append(l[0].DrawLatex( m2deta[i],m2dphi[i],'ap' ) )


  if output == "null":
    raw_input('press Enter to continue...')
  else:
    can.Print(output)

  return [can,eMap,l]





#
# usage
def usage():
  usage_str = '''\necalMapPlotter [options] data_file.root\n
    \tdata_file.root\tfile created by MonoRecAnalysis EDAnalyzer\n
    \t-h, --help \tProduces this help message.\n
    \t-n  --event\tEvent number in data file (first,second, etc)\n
    \t-s, --sub  \tSub-directory in root file. Default \"demo\"\n
    \t-o, --output \tOutput to dump map\n'''

  print usage_str





#
# a function to parse command line options
def do_opts():
  import sys,os,getopt

  try:
    opts,args = getopt.getopt(sys.argv[1:],"hn:s:o:",["help","event=","sub=","output="])
  except getopt.GetoptError, err:
    print str(err)
    usage()
    sys.exit(2)

  global event_number
  global sub_directory
  global output

  for o, a in opts:
    if o in ("-h","--help"):
      usage()
      sys.exit()
    elif o in ("-n","--event"):
      event_number = int(a)
    elif o in ("-s","--sub"):
      sub_directory = a
    elif o in ("-o","--output"):
      output = a
    else:
      assert False, "unhandled option"

  if len(args) != 1:
    print "one argument is required"
    usage()
    sys.exit()

  global data_file
  data_file = args[0]






if __name__ == "__main__":
  do_opts() 

  try:
    ecalMapPlotter(data_file,event_number,sub_directory,output)
  except Exception, err:
    print 'An error occured during ecalMapPlotter'
    print err
    sys.exit(1)


