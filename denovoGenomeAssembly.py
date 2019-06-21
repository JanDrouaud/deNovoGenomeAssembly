#!/users/jdrouau/opt/python_virtualenvs/main-3.6.3/bin/python3.6
# -*- coding: utf8 -*-

# # # # # import clauses
import matplotlib
matplotlib.use('Agg')
import \
  os, sys, subprocess, collections, itertools, datetime, shutil, argparse, math, functools, operator, Bio.SeqIO, re, inspect, socket, psutil, \
  multiprocessing, tempfile, pathlib, gzip, random, string, gzip, magic, matplotlib.pyplot, matplotlib.backends.backend_pdf, seaborn, fileinput

# # # # # define functions
######################################################################################################
def multiLineInput():
  print("paste your text, then Enter and Ctrl-D")
  out=[]
  while True:
    try: line=input()
    except EOFError: break
    out.append(line.rstrip('\\'))
  return ''.join(out).split()

######################################################################################################
def envModuleLoad(modulefile):
  module('load',modulefile)

######################################################################################################
def envModuleUnload(modulefile):
  module('unload',modulefile)

######################################################################################################
def _copy(self, target):
  import shutil
  assert self.is_file()
  shutil.copy(str(self), str(target))

pathlib.Path.copy = _copy

######################################################################################################
def _prefix(self):
  return (self.name).split('.')[0] or self.name

pathlib.Path.prefix = _prefix

######################################################################################################
def _commonsuffix(self):
  if all(isinstance(e,str) for e in self): return os.path.commonprefix(list(e[::-1] for e in self))[::-1]
  elif all(isinstance(e,pathlib.Path) for e in self): return os.path.commonprefix(list(str(e)[::-1] for e in self))[::-1]
  else: return None

os.path.commonsuffix= _commonsuffix
pathlib.Path.commonsuffix= _commonsuffix

######################################################################################################
def _gzOpen(filePath,*nargs,**kwargs):
  import gzip, magic
  if 'mode' not in kwargs: mode=nargs[0]
  if mode[0]=='r' and filePath.is_file() and magic.from_file(str(filePath)).split(', ')[0]!='gzip compressed data':
    return filePath.open(*nargs,**kwargs)
  else:
    return gzip.open(str(filePath),*nargs,**kwargs)

pathlib.Path.gzOpen = _gzOpen

######################################################################################################
def _pigzOpen(filePath,*nargs,**kwargs):
  if 'mode' in kwargs: mode=kwargs['mode']
  elif len(nargs)>0: mode=nargs[0]
  else: log('pigzOpen | mode is not provided') ; return None
  if not (mode in ['r','w']): log('pigzOpen | mode is not valid') ; return None
  if 'cpu' in kwargs: cpu=kwargs['cpu']
  elif len(nargs)>1: cpu=nargs[1]
  else: cpu=multiprocessing.cpu_count()-4
  if mode =='r':
    if filePath.is_symlink(): filePath=filePath.resolve()
    if not filePath.is_file(): log('pigzOpen | file does not exist') ; return None
    elif magic.from_file(str(filePath)).split(', ')[0]!='gzip compressed data': return filePath.open('rb')
    else: a=subprocess.Popen('pigz -f -k -c -d -p '+str(cpu),stdin=filePath.open('rb'),stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True) ; return a.stdout
  else: a=subprocess.Popen('pigz -f -k -c -p'+str(cpu),stdin=subprocess.PIPE,stdout=filePath.open('wb'),stderr=subprocess.PIPE,shell=True) ; return a.stdin

pathlib.Path.pigzOpen = _pigzOpen

######################################################################################################
def pigzMerge(inputFpsList=None,outputFp=None,cpu=None,log=True):
  if log: wlog(functionCallInfo(frame=inspect.currentframe()))
  outputFp=pathlib.Path(outputFp) ; inputFpsList=list(pathlib.Path(fp) for fp in inputFpsList)
  if not(outputFp.is_file() and outputFp.stat().st_size>0):
    if len(inputFpsList)==1: inputFpsList[0].rename(outputFp)
    else:
      with outputFp.pigzOpen(mode='w',cpu=cpu) as outputFh:
        for inputFp in inputFpsList: outputFh.writelines(inputFp.pigzOpen(mode='r',cpu=cpu).readlines())
  return outputFp

######################################################################################################
def _pbzOpen(filePath,*nargs,**kwargs):
  if 'mode' in kwargs: mode=kwargs['mode']
  elif len(nargs)>0: mode=nargs[0]
  else: log('pbzOpen | mode is not provided') ; return None
  if not (mode in ['r','w']): log('pbzOpen | mode is not valid') ; return None
  if 'cpu' in kwargs: cpu=kwargs['cpu']
  elif len(nargs)>1: cpu=nargs[1]
  else: cpu=multiprocessing.cpu_count()-4
  if mode =='r':
    if filePath.is_symlink(): filePath=filePath.resolve()
    if not filePath.is_file(): log('pbzip | file does not exist') ; return None
    elif magic.from_file(str(filePath)).split(', ')[0]!='bzip2 compressed data': return filePath.open('rb')
    else: a=subprocess.Popen('pbzip2 -f -k -c -d -v -p'+str(cpu),stdin=filePath.open('rb'),stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True) ; return a.stdout
  else: a=subprocess.Popen('pbzip2 -f -k -c -z -v -p'+str(cpu),stdin=subprocess.PIPE,stdout=filePath.open('wb'),stderr=subprocess.PIPE,shell=True) ; return a.stdin

pathlib.Path.pbzOpen = _pbzOpen

######################################################################################################
def pbzMerge(inputFpsList=None,outputFp=None,cpu=None,log=True):
  if log: wlog(functionCallInfo(frame=inspect.currentframe()))
  if not(outputFp.is_file() and outputFp.stat().st_size>0):
    if len(inputFpsList)==1: inputFpsList[0].rename(outputFp)
    else:
      with outputFp.pbzOpen(mode='w',cpu=cpu) as outputFh:
        for inputFp in inputFpsList: outputFh.writelines(inputFp.pbzOpen(mode='r',cpu=cpu).readlines())
  return outputFp

######################################################################################################
def timeStamp():
  return datetime.datetime.now().strftime('%Y%m%d%H%M%S')

######################################################################################################
def dictFlatten(dd,sep='|',pref=''):
  return \
    {pref+sep+k if pref else k:v
      for kk,vv in dd.items()
      for k,v in dictFlatten(vv,sep,kk).items()} \
    if isinstance(dd, dict) \
    else {pref+sep+k if pref else k:v
      for kk, vv in enumerate(dd)
      for k,v in dictFlatten(vv,sep,str(kk)).items()} \
    if isinstance(dd,list) \
    else {pref:dd}

######################################################################################################
def rndTmpFileName(tmpDirName='/tmp'):
  return tmpDirName+'/'+''.join(random.sample(string.digits,10))

######################################################################################################
def floatRange(start=None,stop=None,step=None):
  stepSign=(step>0)*2-1
  return [start]+floatRange(start=round(start+step,8),stop=stop,step=step) if (start+step)*stepSign<stop*stepSign else [start]

######################################################################################################
def filesCat(inputFpsList=None,outputFp=None):
  if outputFp==None: outputFile=(lambda tF:(tF,tF.name))(tempfile.NamedTemporaryFile(mode='w+b',dir=str(workDirPaths['tmp']),delete=False))
  else: outputFile=(outputFp.open('rb'),outputFp)
  for fp in inputFpsList: outputFile[0].writelines(fp.open('rb').readlines()) ; outputFile[0].flush()
  outputFile[0].close()
  return outputFile[1]

######################################################################################################
def fillParser(parser=None):
  print(timeStamp(),'fillParser')
  parser.add_argument("--projectName", action="store", dest='projectName', type=str, required=True, help="Name of the project.")
  parser.add_argument("--projectsRootPath", action="store", dest='projectsRootPath', type=pathlib.Path, required=True, help="Full path to the projects directory.")
  parser.add_argument("--projectMetaDataFp", action="store", dest='projectMetaDataFp', type=pathlib.Path, required=False, help="Path to the file holding the project metadata.")
  parser.add_argument("--inputFqFps", action="store", dest='inputFqFps', type=pathlib.Path, nargs='*', default=[], required=False, help="Space separated list of the full path(es) to the files to process.")
  parser.add_argument("--fqSampleSize", action="store", dest='fqSampleSize', type=int, default=100000, required=False, help="Number of reads per sample for parameters tuning.")
  parser.add_argument("--minProb", action="store", dest='minProb', type=float, default=0.95, required=False, help="Minimum probability of error free fragments to consider when preprocessing raw data.")
  parser.add_argument("--minLen", action="store", dest='minLen', type=int, default=20, required=False, help="Minimum length of read fragments for k-mers hashing.")
  parser.add_argument("--chunkSize", action="store", dest='chunkSize', type=int, default=100000, required=False, help="Number of reads per chunk.")
  parser.add_argument("--cpu", action='store', dest='cpu', type=int, default=multiprocessing.cpu_count()-4, required=False, help="number of cpu to use for parallelized steps.")
  parser.add_argument("--memory", action='store', dest='mem', type=int, default=psutil.virtual_memory().total*0.9//1024**3, required=False, help="available mem in Gb.")
  parser.add_argument("--assemblers",action='store',dest='assemblers',type=str, nargs='*',default=['idba','spades'],choices=['abyss','idba','spades','sga'], required=False, help="assemblers.")
  parser.add_argument(\
    "--refGenomeFaFps", action="store", dest='refGenomeFaFps', type=pathlib.Path, nargs='*', required=False, \
    default=[pathlib.Path(os.environ['HOME'])/'opt/data/public/Phytozome/PhytozomeV12/Lusitatissimum/assembly/Lusitatissimum_200_BGIv1.0.fa', \
    pathlib.Path(os.environ['HOME'])/'opt/data/public/ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/224/295/GCA_000224295.2_ASM22429v2/GCA_000224295.2_ASM22429v2_genomic.fna'], \
    help="Path to the fasta file(s) holding the reference genome sequence.")
  parser.add_argument("--refGenomeGffFp", action="store", dest='refGenomeGffFp', type=pathlib.Path, required=False, help="Path to the gff file(s) holding the reference genome annotation.")
  parser.add_argument("--genomeSize", action='store', dest='genomeSize', type=int, default=None, required=True, help="Reference genome size (nt)")
  parser.add_argument("--verbose", action='store', dest='verbose', type=bool, default=False, required=False, help="quiet/verbose script output.")
  return parser

######################################################################################################
def wlog(msg=None):
  params['logFile'].open('a').write('*'*14+'\n'+timeStamp()+'\n'+'*'*14+'\n'+str(msg).rstrip()+'\n'*2)

######################################################################################################
def functionCallInfo(frame=None):
  return \
    'called function "'+inspect.getframeinfo(frame).function+'" with arguments: '+\
    str(dict((k,v if len(str(v))<=100000 else 'TOO LONG FOR DISPLAY') for k,v in inspect.getargvalues(frame).locals.items()))

######################################################################################################
def initializeDirectories(params=None):
  # ! cannot log the functionCallInfo right now because the log file does not yet exist
  # directory structure
  #~ params['workDirRootPath']=(lambda fp:fp.resolve() if fp)(params['projectsRootPath']/'analysis'/params['projectName'])
  params['workDirRootPath']=(params['projectsRootPath']/'analysis'/params['projectName']).resolve()
  params['projectMetaDataFp']=(params['projectsRootPath']/'metaData'/'projectMetaData').resolve()
  if params['inputFqFps'] == []:
    params['inputFqFps']=\
      (lambda fp:list(fp.resolve().glob('*.fastq*')) if fp.is_dir() else sys.exit("projectName does not match any data directory"))\
      (params['projectsRootPath']/'data'/params['projectName'])
  params['workDirPaths'] = dict(\
    [('root',params['workDirRootPath'])]+\
    list(\
      (lambda fp:(fp.name,fp))(params['workDirRootPath']/subDir) \
      for subDir in ('log','input','reference','preprocessed','assembly','tmp','junk','external','alignment'))+\
    list((lambda fP:(fP.name,fP))((params['workDirRootPath']/'assembly')/assembler) for assembler in params['assemblers']))
  for workDirPath in params['workDirPaths'].values():
    if not workDirPath.is_dir(): workDirPath.mkdir(parents=True)
    workDirPath=workDirPath.resolve()
  params['logFile'] = params['workDirPaths']['log']/('logFile_'+now)
  wlog(functionCallInfo(frame=inspect.currentframe()))
  return params

######################################################################################################
def buildMetaDataDict(params=None,log=True):
  if log: wlog(functionCallInfo(frame=inspect.currentframe()))
  with params['projectMetaDataFp'].open('rt') as fh:
    metaDataKeys=fh.readline().strip().split('\t')
    metaDataList=list(filter(lambda dd:dd['cultivar']==params['projectName'],(dict(zip(metaDataKeys,line.strip().split('\t'))) for line in fh.readlines())))
  if len(metaDataList)==0: sys.exit("projectName does not match any metaData entry")
  if set(map(lambda fp:str(fp.name),params['inputFqFps'])) != set(map(lambda dd:dd['fileName'],metaDataList)):
    sys.exit("file names in metadata and file names in data directory are not identical")
  params['metaDataDict']=dict(\
    (k1,dict(\
      (k2,list(zip(*g2))[1]) \
      for k2,g2 in itertools.groupby(sorted((lambda fn:(re.search('_R([12]).fastq',fn).groups()[0],fn))(dd2['fileName']) for dd2 in g1),key=lambda x:x[0]))) \
    for k1,g1 in itertools.groupby(metaDataList,key=lambda dd:dd['key']))
  del metaDataKeys, metaDataList
  return params

######################################################################################################
def buildInputFps(params=None,log=True):
  if log: wlog(functionCallInfo(frame=inspect.currentframe()))
  # symlink the inputFqFps
  inputFqFpsList = list(\
    (lambda y,z:(y,y.symlink_to(z) if not(y.is_symlink()) else None)[0])(params['workDirPaths']['input']/fp.name,fp) for fp in params['inputFqFps'])
  # merge together all R1 files from one hand. and all R2 files from the other hand, by library
  params['projectTag'],dd=list(e[0] for e in zip(*params['metaDataDict'].items()))
  params['inputFqFpsList']=list(\
    (lambda fp:pbzMerge(inputFpsList=list(filter(lambda e:e.name in l,sorted(inputFqFpsList))),outputFp=fp,cpu=params['cpu']))\
    (params['workDirPaths']['input']/(params['projectTag']+'_'+r+'.fastq.bz2')) \
    for r,l in dd.items())
  del params['inputFqFps']
  params['refGenomeFaFpsDict']=dict(
    ('v'+str(v),(lambda y,z:(y,y.symlink_to(z) if not(y.is_symlink()) else None)[0])(params['workDirPaths']['reference']/fp.name,fp)) \
    for v,fp in enumerate(params['refGenomeFaFps'],start=1))
  del params['refGenomeFaFps']
  #~ for p in ('refGenomeFaFp','refGenomeGffFp'):
    #~ if params[p]: params[p]=(lambda fp:(lambda y,z:(y,y.symlink_to(z) if not(y.is_symlink()) else None)[0])(params['workDirPaths']['reference']/fp.name,fp))(params[p])
  return params

######################################################################################################
def preprocParamsTest(inputFqFpsList=None,minProbRange=None,minLenRange=None,projectTag=None,cpu=None,log=True):
  if log: wlog(functionCallInfo(frame=inspect.currentframe()))
  sampleInputFqFpsList=list(fqSample(inputFqFp=fp,numReads=100000,outputFqFp=fp.parent/(fp.prefix()+'_test'+'.fq.gz'),cpu=24,log=True) for fp in inputFqFpsList)
  wlog('sampleInputFqFpsList : '+str(sampleInputFqFpsList))
  for fp in sampleInputFqFpsList: wlog(str(fp)+'\t'+str(len(fp.pigzOpen(mode='r',cpu=cpu).readlines()))+'\t'+str(fp.stat().st_size))
  for mp,ml in list(itertools.product(minProbRange,minLenRange)):
    if log: wlog('**********************'+'\n'+'minProb : '+str(mp)+'\t'+'minLen : '+str(ml))
    samplePreprocFqFpsList=preprocess(inputFqFpsList=sampleInputFqFpsList,projectTag=projectTag+'_test',chunkSize=1000,minProb=mp,minLen=ml,cpu=cpu,log=True)
    if log:
      wlog(\
        'samplePreprocFqFpsList : '+str(samplePreprocFqFpsList)+'\n'+\
        str(list(str(fp)+'\t'+str(len(fp.pigzOpen(mode='r',cpu=cpu).readlines()))+'\t'+str(fp.stat().st_size)+'\n' for fp in samplePreprocFqFpsList)))
      for fp in samplePreprocFqFpsList: wlog(str(fp)+'\t'+str(len(fp.pigzOpen(mode='r',cpu=cpu).readlines()))+'\t'+str(fp.stat().st_size))

######################################################################################################
def fragSizesDist(inputFqFp=None,cpu=None,log=True):
  if log: wlog(functionCallInfo(frame=inspect.currentframe()))
  if inputFqFp == None: log("fragSizeDist | The input fq file was not provided") ; return None
  inputFqFp=pathlib.Path(inputFqFp)
  if not inputFqFp.is_file(): log("fragSizeDist | The input fq file does not exist") ; return None
  outputFragSizeDistFp=inputFqFp.parent/pathlib.Path(inputFqFp.name.split('.')[0]+'_fragSizeDist.txt')
  if not outputFragSizeDistFp.is_file():
    with inputFqFp.pigzOpen('r',cpu) as inputFqFh, outputFragSizeDistFp.open('w') as outputFragSizeDistFh:
      seqStrsIterator=iter(list(i[1] for i in g)[1] for k,g in itertools.groupby(enumerate(inputFqFh),lambda x:x[0]//4))
      sizesCounter=collections.Counter(itertools.chain(*(map(len,re.findall('[^N]+',seqStr.decode().strip())) for seqStr in seqStrsIterator)))
      outputFragSizeDistFh.writelines(iter(str(k)+'\t'+str(v)+'\n' for k,v in sorted(sizesCounter.items())))
      outputFragSizeDistFh.close() ; inputFqFh.close()
  else: sizesCounter=collections.Counter(dict(map(lambda x:tuple(map(int,x.strip().split())),outputFragSizeDistFp.open('r').readlines())))
  return sizesCounter

#######################################################################################################
def Nx(sizesDist=None,proportion=None,log=True):
  # given a size distribution, this function computes the size Sp \
  # such that the cumulated size of fragments at least Sp long \
  # represent a proportion at least p of the total cumulated size
  if log: wlog(functionCallInfo(frame=inspect.currentframe()))
  if sizesDist==None: log("Nx | The input sizes distribution was not provided") ; return None
  if proportion==None: log("Nx | the proportion of total cumulated size was not provided") ; return None
  try: sizesCounter=collections.Counter(dict(sizesDist))
  except: log("Nx | The input sizes distribution cannot be read") ; return None
  cumulatedSizesCounter=collections.Counter(dict(itertools.accumulate(\
    sorted(((k,k*v) for k,v in sizesCounter.items()),reverse=True),\
    lambda x,y:(y[0],y[1]+x[1]))))
  totalCumulatedSize=sorted(cumulatedSizesCounter.items())[0][1]
  return next(itertools.dropwhile(lambda x:x[1]<proportion*totalCumulatedSize,sorted(cumulatedSizesCounter.items(),reverse=True)))[0]

######################################################################################################
def fqIterToFaIter(fqIter=None,log=False):
  for n,l in enumerate(fqIter):
    if n%4==0: yield b'>'+l[1:]
    elif n%4==1: yield l

######################################################################################################
def faIterRenum(its=None,paired=False,log=False):
  pref=b'P' if paired else b'S'
  if len(its)==2 and paired:
    for n,r in enumerate(zip(*(zip(*[it]*2) for it in its))): yield b'>'+pref+bytes(str(n),'utf-8')+b'/1\n'+r[0][1]+b'>'+pref+bytes(str(n),'utf-8')+b'/2\n'+r[1][1]
  elif len(its)==1 and not paired:
    for n,l in enumerate(zip(*its*2)): yield b'>'+pref+bytes(str(n),'utf-8')+b'/1\n'+l[1]

######################################################################################################
def faRename(inputFaFp=None,outputFaFp=None,log=True):
  if log: wlog(functionCallInfo(frame=inspect.currentframe()))
  if outputFaFp==None: outputFaFile=(lambda tF:(tF,tF.name))(tempfile.NamedTemporaryFile(mode='w+t',dir=str(workDirPaths['tmp']),delete=False))
  else: outputFaFile=(outputFaFp.open('wt'),outputFaFp)
  with inputFaFp.open('rt') as fh: outputFaFile[0].writelines('>'+str(e)+'\n'+faSeq[1]+'\n' for e,faSeq in enumerate(Bio.SeqIO.FastaIO.SimpleFastaParser(fh)))
  outputFaFile[0].close()
  return outputFaFile[1]

######################################################################################################
def faSizeFilter(inputFaFp=None,inputFaLines=None,minSize=None,outputFaFp=None,log=False):
  if log: wlog(functionCallInfo(frame=inspect.currentframe()))
  if (outputFaFp and not(outputFaFp.is_file() and outputFaFp.stat().st_size>0)) or not outputFaFp:
    if not inputFaFp:
      if not inputFaLines: return None
      else: inputFaFh=(lambda tF:(tF,tF.writelines(inputFaLines),tF.seek(0))[0])(tempfile.NamedTemporaryFile(mode='w+t',dir=str(workDirPaths['tmp'])))
    else: inputFaFh=inputFaFp.open('rt')
    outputIterator=iter('>'+h+'\n'+s+'\n' for h,s in Bio.SeqIO.FastaIO.SimpleFastaParser(inputFaFh) if len(s)>=minSize)
  else: return outputFaFp
  if outputFaFp:
    outputFaFp.open('wt').writelines(outputIterator)
    return outputFaFp if outputFaFp.is_file() and outputFaFp.stat().st_size>0 else None
  else: return list(outputIterator)

######################################################################################################
def faRename2(inputFaFp=None,inputFaLines=None,pattern=None,appendSize=True,outputFaFp=None,log=False):
  if log: wlog(functionCallInfo(frame=inspect.currentframe()))
  if (outputFaFp and not(outputFaFp.is_file() and outputFaFp.stat().st_size>0)) or not outputFaFp:
    if not inputFaFp:
      if not inputFaLines: return None
      else: inputFaFh=(lambda tF:(tF,tF.writelines(inputFaLines),tF.seek(0))[0])(tempfile.NamedTemporaryFile(mode='w+t',dir=str(workDirPaths['tmp'])))
    else: inputFaFh=inputFaFp.open('rt')
    outputIterator=iter(\
      '>'+re.match(pattern,h).group()+('_'+str(len(s)) if appendSize else '')+'\n'+s+'\n' \
      for h,s in Bio.SeqIO.FastaIO.SimpleFastaParser(inputFaFh))
  else: return outputFaFp
  if outputFaFp:
    outputFaFp.open('wt').writelines(outputIterator)
    return outputFaFp if outputFaFp.is_file() and outputFaFp.stat().st_size>0 else None
  else: return list(outputIterator)

######################################################################################################
def fqSample(inputFqFp=None,numReads=None,outputFqFp=None,cpu=None,log=True):
  if log: wlog(functionCallInfo(frame=inspect.currentframe()))
  if inputFqFp == None: print("fqSample | The input fq file was not provided") ; return None
  inputFqFp=pathlib.Path(inputFqFp)
  if not inputFqFp.is_file(): print("fqSample | The input fq file does not exist") ; return None
  if not(any(isinstance(numReads,t) for t in (float,int))): print("fqSample | No numeric numReads was provided") ; return None
  numReads=str(int(numReads))
  if not outputFqFp: outputFqFp=inputFqFp.parent/pathlib.Path(inputFqFp.name.split('.')[0]+'.sample'+numReads+''.join(inputFqFp.suffixes))
  #~ wlog('outputFqFp : '+str(outputFqFp))
  if not(outputFqFp.is_file() and outputFqFp.stat().st_size>0):
    inputZipState=magic.from_file(str(inputFqFp)).split(', ')[0]
    if inputZipState=='bzip2 compressed data': inputFqFh=inputFqFp.pbzOpen(mode='r',cpu=cpu)
    elif inputZipState=='gzip compressed data': inputFqFh=inputFqFp.pigzOpen(mode='r',cpu=cpu)
    elif inputZipState=='ASCII text': inputFqFh=inputFqFp.open(mode='rb') ; 
    else: return None
    if outputFqFp.suffix=='.bz2': outputFqFh=outputFqFp.pbzOpen(mode='w',cpu=cpu)
    elif outputFqFp.suffix=='.gz': outputFqFh=outputFqFp.pigzOpen(mode='w',cpu=cpu)
    else: outputFqFh=outputFqFp.open(mode='wb')
    envModuleLoad('seqtk/1.1')
    a=subprocess.Popen('seqtk sample /dev/stdin '+numReads,stdin=inputFqFh,stdout=outputFqFh,stderr=subprocess.PIPE,shell=True,universal_newlines=True)
    out = a.communicate() ; envModuleUnload('seqtk/1.1') ; outputFqFh.close() ; inputFqFh.close()
  return None if not(outputFqFp.is_file() and outputFqFp.stat().st_size>0) else outputFqFp

######################################################################################################
def fqMerge(inputFqIterators=None,chunkNum=None,log=True):
  if log: wlog(functionCallInfo(frame=inspect.currentframe()))
  inputFilesList=list(\
    (lambda tF:(tF,tF.name,tF.writelines(it),tF.flush())[:2])(tempfile.NamedTemporaryFile(mode='w+b',dir=str(workDirPaths['tmp']))) for it in inputFqIterators)
  outputFilesList=list((lambda tF:(tF,tF.name))(tempfile.NamedTemporaryFile(mode='w+b',dir=str(workDirPaths['tmp']))) for i in (0,1,2))
  envModuleLoad('usearch/8.1.1861')
  usearchCmdLine=\
    'usearch -fastq_mergepairs '+inputFilesList[0][1]+' -reverse '+inputFilesList[1][1]+' -fastq_maxdiffpct 10 -fastq_maxdiffs 20 -fastq_minovlen 10 '+\
    ' -fastqout '+outputFilesList[0][1]+' -fastqout_notmerged_fwd '+outputFilesList[1][1]+' -fastqout_notmerged_rev '+outputFilesList[2][1]
  a=subprocess.Popen(usearchCmdLine,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True,universal_newlines=True)
  out,err=a.communicate() ; a.wait()
  envModuleUnload('usearch/8.1.1861')
  for tF in inputFilesList: tF[0].close()
  return list(tf[0].readlines() for tf in outputFilesList)

######################################################################################################
def fqJoin(inputFqIterators=None,chunkNum=None,log=True):
  inputFilesList=list((lambda tF:(tF,tF.name,tF.writelines(it),tF.flush())[:2])(tempfile.NamedTemporaryFile(mode='w+b',dir=str(workDirPaths['tmp']))) for it in inputFqIterators)
  outputFile=(lambda tF:(tF,tF.name))(tempfile.NamedTemporaryFile(mode='w+b',dir=str(workDirPaths['tmp'])))
  envModuleLoad('usearch/8.1.1861')
  a=subprocess.Popen(\
    'usearch -fastq_join '+inputFilesList[0][1]+' -reverse '+inputFilesList[1][1]+' -fastqout '+outputFile[1],\
    stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True,universal_newlines=True)
  out,err=a.communicate() ; a.wait()
  envModuleUnload('usearch/8.1.1861')
  for tF in inputFilesList: tF[0].close()
  return [outputFile[0].readlines()]

######################################################################################################
def fqRevComp(inputFqIterator=None,log=True):
  inputFile=(lambda tF:(tF,tF.name,tF.writelines(inputFqIterator),tF.flush())[:2])(tempfile.NamedTemporaryFile(mode='w+b',dir=str(workDirPaths['tmp'])))
  outputFile=(lambda tF:(tF,tF.name))(tempfile.NamedTemporaryFile(mode='w+b',dir=str(workDirPaths['tmp'])))
  envModuleLoad('usearch/8.1.1861')
  a=subprocess.Popen(\
    'usearch -fastx_revcomp '+inputFile[1]+' -fastqout '+outputFile[1],\
    stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True,universal_newlines=True)
  out,err=a.communicate() ; a.wait()
  envModuleUnload('usearch/8.1.1861')
  inputFile[0].close()
  return outputFile[0]

######################################################################################################
def cutRange(r=None,bs=None,log=False):
  return tuple((lambda y:(y[0],y[-1]+1))(tuple(map(operator.itemgetter(1),g))) for k,g in itertools.groupby(enumerate(set(r)-set(bs)), lambda x:x[0]-x[1]))

######################################################################################################
def bestSplit(frag=None,logs=None,minLog=None,minLen=None,log=False):
  b1,e1=frag
  try:
    b2,e2=next(itertools.dropwhile(lambda f:sum(logs[f[0]:f[1]])<minLog,iter((b2,b2+l) for l in range (e1-1,minLen,-1) for b2 in range(b1,e1-l+1))))
    return tuple(itertools.chain(*map(\
      lambda f1:(f1,) if sum(logs[f1[0]:f1[1]])>=minLog else bestSplit(frag=f1,logs=logs,minLog=minLog,minLen=minLen),\
      tuple(filter(lambda f2:f2[1]-f2[0]>=minLen,((b1,max(b1,b2-1)),(b2,e2),(min(e1,e2+1),e1)))))))
  except StopIteration: return ()

######################################################################################################
def readSplit(fqRead=None,minProb=None,minLen=None,namePrefix=None,log=False):
  sH,sS,qH,qS=map(lambda x:bytearray.strip(bytearray(x)),fqRead) ; minLog=math.log(minProb) ; logs=list(logGoodDict[q] for q in qS) ; frag=(0,len(qS))
  goodFrags=list(filter(lambda x:x[1]-x[0]>=minLen,cutRange(r=range(*frag),bs=set(filter(lambda z:logs[z]<=minLog or sS[z]==78,range(*frag))))))
  bestReadSplit=list(itertools.chain(*map(lambda f:((f,)) if sum(logs[f[0]:f[1]])>minLog else bestSplit(frag=f,logs=logs,minLog=minLog,minLen=minLen),goodFrags)))
  bestSubReadsList=\
    (list([l] for l in enumerate((bestReadSplit.pop(i) for i in (0,-1)),start=1)) if min(len(bestReadSplit),3)>1 else [[],[]])+\
    [list(enumerate(bestReadSplit,start=3))]
  return list(\
    list(itertools.chain(*([bytes(namePrefix+'/'+str(r)+'\n','utf-8')]+list(map(lambda ba:bytes(ba)+b'\n',(sS[b:e],qH,qS[b:e]))) for r,(b,e) in l))) if l else b'' \
    for l in bestSubReadsList)

######################################################################################################
def fqSplit(inputFqIterators=None,chunkNum=None,libName=None,minProb=None,minLen=None,log=True):
  if len(inputFqIterators)==3:
    inputFqIterator=itertools.chain(inputFqIterators[0],fqJoin(inputFqIterators[1:]))
  elif len(inputFqIterators)==1: inputFqIterator=inputFqIterators[0]
  splitOutputIterators=list(map(\
    lambda t:list(itertools.chain(*filter(None,t))),\
    zip(*(\
      readSplit(fqRead=fqRead,minProb=minProb,minLen=minLen,namePrefix='@'+libName+':'+str(chunkNum)+':'+str(r)) \
      for r,fqRead in enumerate(zip(*[iter(inputFqIterator)]*4))))))
  return [splitOutputIterators[2],splitOutputIterators[0],fqRevComp(inputFqIterator=splitOutputIterators[1]).readlines()]

######################################################################################################
def preprocess(inputFqFpsList=None,projectTag=None,chunkSize=100000,minProb=None,minLen=None,cpu=None,log=True):
  if log: wlog(functionCallInfo(frame=inspect.currentframe()))
  mergedFqFpsList=parallelProcess(\
      inputFpsList=inputFqFpsList,\
      outputFpsList=list(workDirPaths['preprocessed']/(projectTag+'_'+r+'.fq.gz') for r in ('0','1','2')),\
      outputZip='gzip',cpu=cpu,linesPerItem=4,chunkSize=chunkSize,fnName='fqMerge',fnArgs=(),log=True)
    #~ for libName,libFqFpsDict in inputFqFpsDict.items())
  if log: wlog('mergedFqFpsList : '+str(mergedFqFpsList))
  correctedFqFpsList=lighter(\
      inputFqFpsList=mergedFqFpsList,\
      outputDirPath=workDirPaths['preprocessed'],\
      kmerSize=99,genomeSize=genomeSize,samplingProp=1,cpu=cpu,log=True)
  if log: wlog('correctedFqFpsList : '+str(correctedFqFpsList))
  joinedFqFpsList=parallelProcess(\
      inputFpsList=correctedFqFpsList[1:],\
      outputFpsList=[workDirPaths['preprocessed']/(projectTag+'_joined.fq.gz')],\
      outputZip='gzip',cpu=cpu,linesPerItem=4,chunkSize=chunkSize,fnName='fqJoin',fnArgs=(),log=True)[0]
  if log: wlog('joinedFqFpsList : '+str(joinedFqFpsList))
  allFqFp=pigzMerge(\
      inputFpsList=[correctedFqFpsList[0],joinedFqFpsList],\
      outputFp=workDirPaths['preprocessed']/(projectTag+'_all.fq.gz'),\
      cpu=cpu,log=True)
  if log: wlog('allFqFp : '+str(allFqFp))
  splittedFqFpsList=parallelProcess(\
      inputFpsList=[allFqFp],\
      outputFpsList=list(workDirPaths['preprocessed']/(projectTag+'_'+tag+'_2.fq.gz') for tag in ('0','1','2')),\
      outputZip='gzip',cpu=cpu,linesPerItem=4,chunkSize=chunkSize,fnName='fqSplit',fnArgs=(projectTag,minProb,minLen),log=True)
  if log: wlog('splittedFqFpsList : '+str(splittedFqFpsList))
  preprocFqFpsList=sorted(splittedFqFpsList)
  if log: wlog('preprocFqFpsList : '+str(preprocFqFpsList))
  return preprocFqFpsList

######################################################################################################
def parallelProcess(inputFpsList=None,outputFpsList=None,outputZip=None,cpu=None,linesPerItem=None,chunkSize=None,fnName=None,fnArgs=None,log=True):
  if log: wlog(functionCallInfo(frame=inspect.currentframe()))
  if not all(fp and fp.is_file() and fp.stat().st_size>0 for fp in outputFpsList):
    inputZipState=list(set(magic.from_file(str(fp)).split(', ')[0] for fp in inputFpsList))[0]
    if inputZipState=='bzip2 compressed data': inputFhs=list(fp.pbzOpen(mode='r',cpu=cpu) for fp in inputFpsList)
    elif inputZipState=='gzip compressed data': inputFhs=list(fp.pigzOpen(mode='r',cpu=cpu) for fp in inputFpsList)
    elif inputZipState=='ASCII text': inputFhs=list(fp.open(mode='rb') for fp in inputFpsList)
    else: return None
    if outputZip=='gzip': outputFhs=list(fp.pigzOpen(mode='w',cpu=cpu) for fp in outputFpsList)
    elif outputZip=='bzip': outputFhs=list(fp.pbzOpen(mode='w',cpu=cpu) for fp in outputFpsList)
    else: outputFhs=list(fp.open(mode='wb') for fp in outputFpsList)
    chunkIterator=iter(list(list(chunk) for chunk in chunkTuple) for chunkTuple in zip(*(zip(*[fh]*chunkSize*linesPerItem) for fh in inputFhs)))
    def writeOut(fnOut):
      for outputFh,outChunk in zip(outputFhs,fnOut): outputFh.writelines(outChunk)
    with multiprocessing.Pool(processes=cpu) as pool:
        for chunkNum,chunk in enumerate(chunkIterator):
          pool.apply_async(globals()[fnName],(chunk,chunkNum,)+fnArgs+((True,) if log else()),callback=writeOut)
        pool.close() ; pool.join()
    for fh in outputFhs: fh.close()
  if all(fp and fp.is_file() and fp.stat().st_size>0 for fp in outputFpsList): return outputFpsList
  else: log("parallel | output file(s) is(are) not present") ; return None

######################################################################################################
def lighter(inputFqFpsList=None,outputDirPath=None,kmerSize=None,genomeSize=None,samplingProp=1,cpu=None,log=True):
  if log: wlog(functionCallInfo(frame=inspect.currentframe()))
  outputFqFpsList=list(fp.parent/(''.join([fp.prefix(),'.cor']+fp.suffixes)) for fp in inputFqFpsList)
  if not all(fp.is_file() and fp.stat().st_size>0 for fp in outputFqFpsList):
    lighterArgsList=list(\
      ('-r',fp) for fp in inputFqFpsList)+[('-k',' '.join(map(str,[kmerSize,genomeSize,samplingProp]))),('-od',outputDirPath),('-t',cpu),('-zlib',5)]
    lighterCmdLine='lighter '+' '.join(itertools.chain.from_iterable((str(k),str(v)) for k,v in lighterArgsList))
    envModuleLoad('lighter/1.1.1')
    a=subprocess.Popen(lighterCmdLine,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True,universal_newlines=True)
    out,err = a.communicate() ; wlog('lighter out'+'\n'+out) ; wlog('lighter err'+'\n'+err)
    envModuleUnload('lighter/1.1.1')
  if not all(fp.is_file() and fp.stat().st_size>0 for fp in outputFqFpsList): return None
  else: return outputFqFpsList

######################################################################################################
def abyssAssembly(readsFqFpsList=None,longFaFp=None,kmerRange=None,target=None,cpu=None,log=True):
  if log: wlog(functionCallInfo(frame=inspect.currentframe()))
  outputDirPathesDict=collections.OrderedDict((k,workDirPaths['abyss']/('k'+str(k))) for k in kmerRange)
  outputFpsDict=collections.OrderedDict((k,outputDirPath/('k'+str(k)+'-'+target+'.fa')) for k,outputDirPath in outputDirPathesDict.items())
  for k,outputDirPath in outputDirPathesDict.items():
    if not outputDirPath.is_dir(): outputDirPath.mkdir()
    outputFp=outputFpsDict[k]
    if not(outputFp.is_file() and outputFp.stat().st_size>0):
      name=projectName+'_k'+str(k)+'_abyss'
      abyssArgsDict=collections.OrderedDict(\
      [('--directory ',outputDirPath)]+([('se=',readsFqFpsList[0])] if readsFqFpsList else [])+([('long=',longFaFp)] if longFaFp else [])+\
      ([('lib=','\'l1\''),('l1=','\''+' '.join(map(str,readsFqFpsList[1:]))+'\'')] if len(readsFqFpsList)==3 else [])+\
      [('name=','k'+str(k)),('j=',cpu),('np=',cpu),('a=',2),('c=',3),('b=',2*k+2),('e=',2),('k=',k),('l=',max(50,k))]+\
      [('m=',max(50,kmerSize)),('n=',3),('N=',4),('p=',0.99),('q=',0),('Q=',0),('s=',2*k),('S=',2*k),('t=',k)]+\
      [('POPBUBBLES_OPTIONS','"--branches=2 --bubble-length='+str(2*k+2)+' --identity='+str(math.floor((1-1/k)*1000)/1000)+' --coverage=1.5"')]+\
      [('SIMPLEGRAPH_OPTIONS=','"--no-scaffold"'),('OVERLAP_OPTIONS=','"--no-scaffold"'),('MERGEPATH_OPTIONS=','"--greedy"'),('',target)])
      abyssCmdLine='abyss-pe '+' '.join(list(str(k)+str(v) for k,v in abyssArgsDict.items()))
      log('abyssCmdLine : '+abyssCmdLine)
      envModuleLoad('abyss/2.0.2_mpich') ; envModuleLoad('bwa/0.7.17')
      a=subprocess.Popen(abyssCmdLine,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True,universal_newlines=True)
      out,err = a.communicate() ; wlog('abyss out'+'\n'+out) ; wlog('abyss err'+'\n'+err)
      envModuleUnload('bwa/0.7.17') ; envModuleUnload('abyss/2.0.2_mpich')
  return None if not all(outputFp.is_file() and outputFp.stat().st_size>0 for outputFp in outputFpsDict.values()) else outputFp

######################################################################################################
def spadesAssembly(readsFqFpsList=None,kmersRange=None,trustedContigsFaFp=None,untrustedContigsFaFp=None,outputDirPath=None,cont=None,cpu=None,mem=None,log=True):
  if log: wlog(functionCallInfo(frame=inspect.currentframe()))
  outputFpsDict={'contigs':outputDirPath/'contigs.fasta','scaffolds':outputDirPath/'scaffolds.fasta'}
  if cont: spadesCmdLine='spades.py --continue -o '+str(outputDirPath)
  else: 
    spadesArgsDict=collections.OrderedDict(\
      [('-1',readsFqFpsList[1]),('-2',readsFqFpsList[2]),('-s',readsFqFpsList[0]),('--threads',cpu),('--memory',mem)]+\
      [('--tmp-dir',workDirPaths['tmp']),('-o',outputDirPath),('--only-assembler',''),('--cov-cutoff',3),('-k',','.join(map(str,kmersRange)))]+\
      ([('--trusted-contigs',trustedContigsFaFp)] if trustedContigsFaFp else [])+\
      ([('--untrusted-contigs',untrustedContigsFaFp)] if untrustedContigsFaFp else []))
    spadesCmdLine='spades.py '+' '.join(list(str(k)+' '+str(v) for k,v in spadesArgsDict.items()))
  wlog('spadesCmdLine : '+spadesCmdLine)
  if not all(fp.is_file() and fp.stat().st_size>0 for fp in outputFpsDict.values()):
    envModuleLoad('spades/3.11.1')
    a=subprocess.Popen(spadesCmdLine,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True,universal_newlines=True)
    out,err = a.communicate() ; wlog('spades out'+'\n'+out) ; wlog('spades err'+'\n'+err)
    envModuleUnload('spades/3.11.1')
  return outputFpsDict if all(fp.is_file() and fp.stat().st_size>0 for fp in outputFpsDict.values()) else None

######################################################################################################
def idbaAssembly(readsFqFpsList=None,analysisPrefix=None,kmersRange=None,longFaFp=None,outputDirPath=None,cpu=None,log=True):
  if log: wlog(functionCallInfo(frame=inspect.currentframe()))
  outputFpsDict={'contigs':outputDirPath/'contig.fa','scaffolds':outputDirPath/'scaffold.fa'}
  if not all(fp.is_file() and fp.stat().st_size>0 for fp in outputFpsDict.values()):
    mergedFaReadsFp=workDirPaths['preprocessed']/(analysisPrefix+'_mergedReads.fa')
    if not (mergedFaReadsFp.is_file() and mergedFaReadsFp.stat().st_size>0):
      mergedFaReadsFp.open('wb').writelines(\
          (lambda its:itertools.chain(faIterRenum(its=[its[0]],paired=False),faIterRenum(its=its[1:],paired=True)))\
          (list(fqIterToFaIter(fqIter=fp.pigzOpen('r',cpu=cpu)) for fp in readsFqFpsList)))
    mink,maxk,stepk=(lambda r:(r[0],r[-1],r.step))(kmersRange)
    idbaArgsDict=collections.OrderedDict(\
      [('--out',outputDirPath),('--read',mergedFaReadsFp)] +([('--long_read',longFaFp)] if longFaFp else [])+\
      [('--mink',mink),('--maxk',maxk),('--step',stepk)]+[('--inner_mink',mink),('--inner_step',stepk)]+[('--prefix',3),('--min_count',3),('--min_support',3)]+\
      [('--num_threads',cpu),('--seed_kmer',36),('--min_contig',2*maxk),('--similar',str(math.floor((1-1/maxk)*1000)/1000)),('--min_pairs',3)]+\
      [('--no_bubble',''),('--no_correct','')])
    idbaCmdLine='idba_ud  '+' '.join(str(k)+' '+str(v) for k,v in idbaArgsDict.items())
    wlog('idbaCmdLine : '+idbaCmdLine)
    envModuleLoad('idba/1.1.3')
    a=subprocess.Popen(idbaCmdLine,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True,universal_newlines=True)
    out,err = a.communicate() ; wlog('idba out'+'\n'+out) ; wlog('idba err'+'\n'+err)
    envModuleUnload('idba/1.1.3')
  return outputFpsDict if all(fp.is_file() and fp.stat().st_size>0 for fp in outputFpsDict.values()) else None

######################################################################################################
def sgaAssembly(readsFqFpsList=None,minOverlap=50,prefix=None,cpu=None,log=True):
  if log: wlog(functionCallInfo(frame=inspect.currentframe()))
  outPrefix=prefix.lower()+'SgaMerge' ; outputFaFp=inputFaFp.parent/(outPrefix+'-contigs.fa')
  mergedInputFqTmpfile=\
    (lambda tF:(tF,tF.name,pigzMerge(inputFpsList=readsFqFpsList,outputFp=tF.name,cpu=cpu))[:2])\
    (tempfile.NamedTemporaryFile(mode='w+b',dir=str(workDirPaths['tmp'])))
  tmpFifoFile=tmpFifoFile=rndTmpFileName(tmpDirName=str(workDirPaths['tmp']))
  os.system('mkfifo '+tmpFifoFile+' ; pigz -dck '+' '.join(map(str,readsFqFpsList))+' > '+tmpFifoFile+ ' &')
  if not(outputFaFp.is_file() and outputFaFp.stat().st_size>0):
    sgaArgsDict=collections.OrderedDict([\
    ('index',collections.OrderedDict(\
      [('--verbose',''),('--algorithm=','sais'),('--disk=',1000000)]+\
      [('--threads=',cpu),('--gap-array=',32),('',tmpFifoFile)])),\
    ('rmdup',collections.OrderedDict([('--verbose',''),('--error-rate=',0),('--threads=',cpu),('--sample-rate=',8),('',inputFaFp)])),\
    ('overlap',collections.OrderedDict(\
      [('--verbose',''),('--threads=',cpu),('--error-rate=',0),('--min-overlap=',int(minOverlap))]+\
      [('--sample-rate=',8),('',(lambda fp:fp.parent/(fp.stem+'.rmdup.fa'))(inputFaFp))])),\
    ('assemble',collections.OrderedDict(\
      [('--verbose',''),('--min-overlap=',int(minOverlap)),('--max-edges=',512),('--bubble=',0),('--max-divergence=',0)]+\
      [('--max-gap-divergence=',0),('--max-indel=',1000),('--out-prefix=',outPrefix),('--cut-terminal=',10)]+\
      [('--min-branch-length=',50),('',(lambda fp:fp.parent/(fp.stem+'.rmdup.asqq.gz'))(inputFaFp))]))])
    sgaCmdLine=' ; '.join('sga '+sp+' '+' '.join(str(k)+str(v) for k,v in dd.items()) for sp,dd in sgaArgsDict.items())
    log('sgaCmdLine : '+sgaCmdLine)
    envModuleLoad('sga/0.10.15')
    a=subprocess.Popen(sgaCmdLine,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True,universal_newlines=True)
    out,err = a.communicate() ; wlog('sga out'+'\n'+out) ; wlog('sga err'+'\n'+err)
    envModuleUnload('sga/0.10.15')
  return outputFaFp if (outputFaFp.is_file() and outputFaFp.stat().st_size>0) else None

######################################################################################################
def bwa(queryCtgsFaFp=None,subjectCtgsFaFp=None,workDirPath=None,log=True):
  if log: wlog(functionCallInfo(frame=inspect.currentframe()))
  queryCtgLens=dict((ctgId.split(' ')[0],len(ctgSeq)) for ctgId,ctgSeq in Bio.SeqIO.FastaIO.SimpleFastaParser(queryCtgsFaFp.open('rt')))
  subjectCtgLens=dict((ctgId.split(' ')[0],len(ctgSeq)) for ctgId,ctgSeq in Bio.SeqIO.FastaIO.SimpleFastaParser(subjectCtgsFaFp.open('rt')))
  bwaDbPrefixPath=workDirPath/subjectCtgsFaFp.stem ; bwaDbFps=list(pathlib.Path(str(bwaDbPrefixPath)+'.'+suffix) for suffix in ('bwt','pac','ann','amb','sa'))
  outputSamFp=workDirPath/(queryCtgsFaFp.stem+'_vs_'+subjectCtgsFaFp.stem+'.sam')
  if not all(fp.is_file() and fp.stat().st_size>0 for fp in bwaDbFps):
    bwaCmdLine='bwa index -p '+str(workDirPath/subjectCtgsFaFp.stem)+' '+str(subjectCtgsFaFp)
    log('bwaCmdLine : '+bwaCmdLine)
    envModuleLoad('bwa/0.7.17')
    a=subprocess.Popen(bwaCmdLine,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True,universal_newlines=True)
    out,err = a.communicate() ; log('bwa out'+'\n'+out) ; log('bwa err'+'\n'+err)
    envModuleUnload('bwa/0.7.17')
  if not all(fp.is_file() and fp.stat().st_size>0 for fp in bwaDbFps): return None
  if not (outputSamFp.is_file() and outputSamFp.stat().st_size>0):
    tmpOutputSamFile=(lambda tF:(tF,tF.name))(tempfile.NamedTemporaryFile(mode='w+b',dir=str(workDirPaths['tmp'])))
    bwaCmdLine='bwa mem -k 30 -r 1.1 -w 0 -d 50 -T 30 -A 1 -B 8 -O50 -E 50 -c 50 -t '+str(cpu)+' -o '+tmpOutputSamFile[1]+' '+str(bwaDbPrefixPath)
    log('bwaCmdLine : '+bwaCmdLine)
    envModuleLoad('bwa/0.7.17')
    a=subprocess.Popen(bwaCmdLine,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True,universal_newlines=True)
    out,err = a.communicate() ; wlog('bwa out'+'\n'+out) ; wlog('bwa err'+'\n'+err)
    envModuleUnload('bwa/0.7.17')
    if not pathlib.Path(tmpOutputSamFile).is_file() and pathlib.Path(tmpOutputSamFile).stat().st_size>0: return None
    with open(tmpOutputSamFile[0],'rt') as inFh, outputSamFp.open('rt') as outFh:
      for line in inFh.readlines():
        flds=line.strip().split('\t')
        if line[0]=='@' or flds[1]=='4': outFh.write(line) ; continue
        queryCtgId=flds[0] ; queryCtgLen=queryCtgLens[queryCtgId] ; subjCtgId=flds[3] ;subjCtgLen=subjectCtgLens[subjCtgId]
        cigarFlds=re.findall('[0-9]+[MSHDINPX]',flds[5]) ; cigarSums=collections.Counter(itertools.chain(*((op[-1]*int(op[:-1])) for op in cigarFlds)))
        tags=dict((lambda l:(l[0],l[-1]))(tagStr.split(':'))for tagStr in flds[11:])
        if all(k in ('M','S') for k in cigarSums.keys()) and cigarSums['M']+cigarSums['S']==queryCtgLen \
        and tags['NM']==0 and tags['MD']==cigarSums['M'] and tags['AS']==cigarSums['M'] \
        and ((len(cigarFlds)==2 and cigarSums['M']>0 and cigarSums['S']>0) or (len(cigarFlds)==3 and cigarSums['M']==subjCtgLen and queryCtgLen>subjCtgLen)):
          outFh.write(line) ; continue
    if not (outputSamFp.is_file() and outputSamFp.stat().st_size>0): return None
    else: return outputSamFp

######################################################################################################
def cerulean(subjectCtgsFaFp=None,subjectCtgsDotFp=None,queryCtgsFaFp=None,workDirPath=None,cpu=None,log=True):
  if log: wlog(functionCallInfo(frame=inspect.currentframe()))
  outputDirPath=workDirPaths['cerulean']
  outputFpsList=list(outputDirPath.glob('*cerulean*'))
  if len(outputFpsList)==0 or all(fp.stat().st_size==0 for fp in outputFpsList):
    subjectCtgsFaAliasFp=(lambda y,z:(y,y.symlink_to(z) if not(y.is_symlink()) else None)[0])(outputDirPath/(subjectCtgsFaFp.stem+'-contigs.fa'),subjectCtgsFaFp)
    subjectCtgsDotAliasFp=(lambda y,z:(y,y.symlink_to(z) if not(y.is_symlink()) else None)[0])(outputDirPath/(subjectCtgsDotFp.stem+'-contigs.dot'),subjectCtgsDotFp)
    queryCtgsFaAliasFp=(lambda y,z:(y,y.symlink_to(z) if not(y.is_symlink()) else None)[0])(outputDirPath/queryCtgsFaFp.name,queryCtgsFaFp)
    sawriterOutputFp=outputDirPath/(subjectCtgsFaAliasFp.stem+'.sa')
    blasrOutputFp=outputDirPath/(subjectCtgsFaAliasFp.stem+'_pacbio_contigs_mapping.fasta.m4')
    sawriterCmdLine='sawriter '+str(sawriterOutputFp)+' '+str(subjectCtgsFaAliasFp)
    blasrArgsDict=collections.OrderedDict([\
      ('--sa',sawriterOutputFp),('--minMatch',30),('--minPctIdentity',100),('--bestn',30),('--hitPolicy','allbest'),('--randomSeed',0),('-m',4),\
      ('--header',''),('--scoreMatrix',' '.join((['-5']+['100']*5)*4+['100'])),('--nCandidates',30),('--affineOpen',100),('--affineExtend',100),\
      ('-maxScore',-200),('-nproc',cpu),('--noSplitSubreads',''),('-out',blasrOutputFp),('--unaligned',outputDirPath/(queryCtgsFaAliasFp.stem+'_unaligned.fa'))])
    blasrCmdLine='blasr '+str(queryCtgsFaAliasFp)+' '+str(subjectCtgsFaAliasFp)+' '+' '.join(str(k)+' '+str(v) for k,v in blasrArgsDict.items())
    ceruleanArgsDict=collections.OrderedDict([('--nproc',cpu),('--dataname',subjectCtgsFaFp.stem),('--basedir',outputDirPath)])
    ceruleanCmdLine='/usr/bin/python /users/jdrouau/opt/softwares/cerulean/0.1/bin/Cerulean.py '+' '.join(str(k)+' '+str(v) for k,v in ceruleanArgsDict.items())
    wholeCmdLine=' ; '.join([sawriterCmdLine,blasrCmdLine,ceruleanCmdLine])
    wlog('wholeCmdLine : '+wholeCmdLine)
    envModuleLoad('blasr/20180214') ; envModuleLoad('cerulean/0.1')
    a=subprocess.Popen(wholeCmdLine,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True,universal_newlines=True)
    out,err = a.communicate() ; wlog('sawriter/blasr/cerulean out'+'\n'+out) ; wlog('sawriter/blasr/cerulean err'+'\n'+err)
    envModuleUnload('cerulean/0.1') ; envModuleUnload('blasr/20180214')
  outputFpsList=list(outputDirPath.glob('*cerulean*'))
  if len(outputFpsList)==0 or all(fp.stat().st_size==0 for fp in outputFpsList): return None
  else: return outputFpsList

######################################################################################################
def mummerNucmer(referenceFaFp=None,queryFaFp=None,outputDirPath=None,prefix=None,cpu=None,log=True):
  if log: wlog(functionCallInfo(frame=inspect.currentframe()))
  if not prefix: prefix=queryFaFp.prefix()+'_vs_'+referenceFaFp.prefix()
  deltaFp=outputDirPath/(prefix+'.delta')
  if not(deltaFp.is_file() and deltaFp.stat().st_size>0):
    nucmerArgsDict=collections.OrderedDict(\
      [('--maxmatch',''),('--breaklen=',500),('--mincluster=',500),('--diagdiff=',5),('--diagfactor=',0.05),('--maxgap=',500),('--minmatch=',20)]+\
      [('--minalign=',500),('--nosimplify',''),('--prefix=',outputDirPath/prefix),('--threads=',cpu),(referenceFaFp,''),(queryFaFp,'')])
    nucmerCmdLine='nucmer  '+' '.join(str(k)+str(v) for k,v in nucmerArgsDict.items())
    envModuleLoad('mummer/4.0.0beta2')
    a=subprocess.Popen(nucmerCmdLine,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True,universal_newlines=True)
    out,err = a.communicate() ; wlog('nucmer out'+'\n'+out) ; wlog('nucmer err'+'\n'+err)
    envModuleUnload('mummer/4.0.0beta2')
  return None if not(deltaFp.is_file() and deltaFp.stat().st_size>0) else deltaFp

######################################################################################################
def mummerDeltaFilter(inputDeltaFp=None,outputDirPath=None,prefix=None,cpu=None,log=True):
  if log: wlog(functionCallInfo(frame=inspect.currentframe()))
  #~ if not prefix: prefix=inputDeltaFp.prefix()
  outputDeltaFp=outputDirPath/(inputDeltaFp.name+'.filtered')
  if not(outputDeltaFp.is_file() and outputDeltaFp.stat().st_size>0):
    deltaFilterCmdLine='delta-filter -1 -i 95 -l 500 -u 0 -o 100 '+str(inputDeltaFp)+' > '+str(outputDeltaFp)
    envModuleLoad('mummer/4.0.0beta2')
    a=subprocess.Popen(deltaFilterCmdLine,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True,universal_newlines=True)
    out,err = a.communicate() ; wlog('delta-filter out'+'\n'+out) ; wlog('delta-filter err'+'\n'+err)
    envModuleUnload('mummer/4.0.0beta2')
  return outputDeltaFp if (outputDeltaFp.is_file() and outputDeltaFp.stat().st_size>0) else None

######################################################################################################
def mummerShowCoords(deltaFp=None,outputDirPath=None,prefix=None,cpu=None,log=True):
  if log: wlog(functionCallInfo(frame=inspect.currentframe()))
  if not prefix: prefix=deltaFp.prefix()
  coordsFp=outputDirPath/(prefix+'.coords')
  if not(coordsFp.is_file() and coordsFp.stat().st_size>0):
    showCoordsCmdLine='show-coords -c -I 95 -r '+str(deltaFp)+' > '+str(coordsFp)
    envModuleLoad('mummer/4.0.0beta2')
    a=subprocess.Popen(showCoordsCmdLine,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True,universal_newlines=True)
    out,err = a.communicate() ; wlog('show-coords out'+'\n'+out) ; wlog('show-coords err'+'\n'+err)
    envModuleUnload('mummer/4.0.0beta2')
  return coordsFp if (coordsFp.is_file() and coordsFp.stat().st_size>0) else None

######################################################################################################
def mummerShowSnps(deltaFp=None,coordsFp=None,outputDirPath=None,prefix=None,cpu=None,log=True):
  if log: wlog(functionCallInfo(frame=inspect.currentframe()))
  if not prefix: prefix=deltaFp.prefix()
  snpsFp=outputDirPath/(prefix+'.snps')
  if not(snpsFp.is_file() and snpsFp.stat().st_size>0):
    showSnpsCmdLine='show-snps -C -I -l -r -T '+str(deltaFp)+' > '+str(snpsFp)
    envModuleLoad('mummer/4.0.0beta2')
    a=subprocess.Popen(showSnpsCmdLine,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True,universal_newlines=True)
    out,err = a.communicate() ; wlog('show-snps out'+'\n'+out) ; wlog('show-snps err'+'\n'+err)
    envModuleUnload('mummer/4.0.0beta2')
  return snpsFp if (snpsFp.is_file() and snpsFp.stat().st_size>0) else None

######################################################################################################
def mummerShowAligns(deltaFp=None,refId=None,qryId=None,log=False):
  if log: wlog(functionCallInfo(frame=inspect.currentframe()))
  print(refId,qryId)
  showAlignsCmdLine='show-aligns -w 100 '+str(deltaFp)+' '+refId+' '+qryId
  envModuleLoad('mummer/4.0.0beta2')
  a=subprocess.Popen(showAlignsCmdLine,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True,universal_newlines=True)
  mummerAln,err = a.communicate()
  envModuleUnload('mummer/4.0.0beta2')
  return mummerAln

######################################################################################################
def mummerAlnParse(mummerAln=None,log=False):
  if log: wlog(functionCallInfo(frame=inspect.currentframe()))
  numAlnLines=list(enumerate(mummerAln.split('\n')[3:-3]))
  (refId,qryId)=(lambda f:(f[3],f[5]))(numAlnLines[0][1].split(' '))
  alnBlocks=list(\
    ((b[0][0]+2,b[1][0]-1,)+re.fullmatch('-- BEGIN alignment \[ ([+-])1 (\d+) - (\d+) \| ([+-])1 (\d+) - (\d+) \]',b[0][1]).groups()) \
    for b in zip(*[iter(l for l in numAlnLines[2:] if l[1][:2]=='--')]*2))
  return list(\
    (lambda seqsPair:(refId,refStart,refStop,refStrand,seqsPair[0],qryId,qryStart,qryStop,qryStrand,seqsPair[1]))\
    (tuple(map(lambda lg:''.join(lg),zip(*(tuple(e[1][11:] for e in linesGrp[1:3]) for linesGrp in zip(*[iter(numAlnLines[beginLineNum:endLineNum])]*4)))))) \
    for (beginLineNum,endLineNum,refStrand,refStart,refStop,qryStrand,qryStart,qryStop) in alnBlocks)

######################################################################################################
def mummerAlign(deltaFp=None,ctgPairsFp=None,outputFp=None,log=False):
  if log: wlog(functionCallInfo(frame=inspect.currentframe()))
  if not outputFp: outputFp=deltaFp.parent/(deltaFp.prefix()+'.alns')
  if not(outputFp.is_file() and outputFp.stat().st_size>0):
    with outputFp.open('wt') as ofh:
      #~ ctgPairsList=list(tuple(line.strip().split('\t')) for line in ctgPairsFp.open('rt').readlines())
      ofh.writelines(itertools.chain.from_iterable(\
        map(lambda t:'\t'.join(t)+'\n',mummerAlnParse(mummerAln=mummerShowAligns(deltaFp=deltaFp,refId=refId,qryId=qryId))) \
        for refId,qryId in list(tuple(line.strip().split('\t')) for line in ctgPairsFp.open('rt').readlines())))
  return outputFp if (outputFp.is_file() and outputFp.stat().st_size>0) else None

######################################################################################################
def mummerAlign2(ctgPairsListList=None,chunkNum=None,deltaFp=None,log=True):
  if log: wlog(functionCallInfo(frame=inspect.currentframe()))
  return [list(itertools.chain.from_iterable(\
    map(lambda t:bytes('\t'.join(t)+'\n','UTF-8'),mummerAlnParse(mummerAln=mummerShowAligns(deltaFp=deltaFp,refId=refId,qryId=qryId))) \
    for refId,qryId in list(tuple(line.decode().strip().split('\t')) for line in ctgPairsListList[0])))]

######################################################################################################
def getCtxtSeqs(aln=None,refPos=None,log=False):
  refId,refStart,refStop,refStrand,refSeq,qryId,qryStart,qryStop,qryStrand,qrySeq=list(f(e) for f,e in zip((str,int,int,str,str,str,int,int,str,str),aln))
  refCtxtSeq,qryCtxtSeq=\
    (lambda refAlnPos:(refSeq[refAlnPos-100:refAlnPos+101],qrySeq[refAlnPos-100:refAlnPos+101]))\
    (list(itertools.filterfalse(lambda p:p[1]=='.',enumerate(refSeq)))[refPos-refStart][0])

######################################################################################################
def coordsFileProcess(coordsFp=None,outputFp=None,log=True):
  if log: wlog(functionCallInfo(frame=inspect.currentframe()))
  if not outputFp: outputFp=coordsFp.parent/(coordsFp.prefix()+'.goodCtgPairs')
  if not(outputFp.is_file() and outputFp.stat().st_size>0):
    def isN(ctgId=None,seq=None,start=None,stop=None):
      if (stop-start)<=21: return True
      elif ((start<1 or stop>len(seq)) and stop-start+1<1): return sum(l=='N' for l in seq[start-1:stop])/(stop-start+1)>0.95
      else: return False
    alnTags=('refStart','refEnd','qryStart','qryEnd','refAlnLen','qryAlnLen','idtyPerc','refCov','qryCov','refId','qryId')
    with coordsFp.open('rt') as ifh, outputFp.open('wt') as ofh:
      line=next(ifh) ; [refFaFp,qryFaFp]=list(pathlib.Path(f) for f in line.split())
      refCtgDict=dict((h,s) for h,s in Bio.SeqIO.FastaIO.SimpleFastaParser(refFaFp.open('rt')))
      qryCtgDict=dict((h,s) for h,s in Bio.SeqIO.FastaIO.SimpleFastaParser(qryFaFp.open('rt')))
      lines=list(itertools.islice(ifh,4))
      alnsDict=collections.OrderedDict((tuple(k),dict([('alns',list(collections.OrderedDict(zip(alnTags[:-2],list(map(int,a[0:6]))+list(map(float,a[6:-2])))) for a in g))])) \
        for k,g in itertools.groupby(sorted(\
          iter(re.findall('([\w\d._]+)',line.strip()) for line in ifh.readlines()),\
          key=lambda e:(e[9],e[10],int(e[0]))),key=lambda e:(e[9],e[10])))
      for (refId,qryId),alns in alnsDict.items():
        alnsList=alns['alns'] ; refCtgSeq=refCtgDict[refId] ; refLen=len(refCtgSeq) ; qryCtgSeq=qryCtgDict[qryId] ; qryLen=len(qryCtgSeq)
        if int(alnsList[0]['qryStart'])>int(alnsList[0]['qryEnd']):
          qryCtgSeq=Bio.Seq.reverse_complement(qryCtgSeq)
          for aln in alnsList: aln['qryStart']=len(qryCtgDict[qryId])-int(aln['qryStart'])+1 ; aln['qryEnd']=len(qryCtgDict[qryId])-int(aln['qryEnd'])+1
        firstAln,lastAln=alnsList[0],alnsList[-1]
        for k,v in firstAln.items(): exec(k+'='+str(v))
        leftBitsList=[\
          (firstAln['refStart']==1),(firstAln['qryStart']==1),(firstAln['refStart']>=firstAln['qryStart']),\
          (firstAln['refStart']<=firstAln['qryStart']),\
          (isN(seq=refCtgSeq,start=firstAln['refStart']-firstAln['qryStart'],stop=firstAln['refStart']-1)),(isN(seq=qryCtgSeq,start=1,stop=firstAln['qryStart']-1)),\
          (isN(seq=qryCtgSeq,start=firstAln['qryStart']-firstAln['refStart'],stop=firstAln['qryStart']-1)),(isN(seq=refCtgSeq,start=1,stop=firstAln['refStart']-1))]
        leftByte=sum(v*2**i for i,v in enumerate(leftBitsList))
        for k,v in lastAln.items(): exec(k+'='+str(v))
        rightBitsList=[\
          (lastAln['refEnd']==refLen),(lastAln['qryEnd']==qryLen),((refLen-lastAln['refEnd'])>(qryLen-lastAln['qryEnd'])),\
          ((refLen-lastAln['refEnd'])<(qryLen-lastAln['qryEnd'])),\
          (isN(seq=refCtgDict[refId],start=lastAln['refEnd']+1,stop=lastAln['refEnd']+(qryLen-lastAln['qryEnd']))),(isN(seq=qryCtgSeq,start=lastAln['qryEnd']+1,stop=qryLen)),\
          (isN(seq=qryCtgDict[qryId],start=lastAln['qryEnd']+1,stop=lastAln['qryEnd']+(refLen-lastAln['refEnd']))),(isN(seq=refCtgSeq,start=lastAln['refEnd']+1,stop=refLen))]
        rightByte=sum(v*2**i for i,v in enumerate(rightBitsList))
        if (leftByte!=4 and leftByte!=12 and rightByte!=4 and rightByte!=12): ofh.write(refId+'\t'+qryId+'\n')
  return outputFp if (outputFp.is_file() and outputFp.stat().st_size>0) else None

######################################################################################################
def snpsFileProcess(snpsFp=None,ctgPairsFp=None,outputFp=None,log=True):
  if log: wlog(functionCallInfo(frame=inspect.currentframe()))
  if not outputFp: outputFp=snpsFp.parent/(snpsFp.prefix()+'.goodSnps')
  if not(outputFp.is_file() and outputFp.stat().st_size>0):
    ctgPairsSet=set(tuple(line.strip().split()[-2:]) for line in ctgPairsFp.open('rt').readlines())
    with snpsFp.open('rt') as ifh, outputFp.open('wt') as ofh:
      lines=list(itertools.islice(ifh,5))
      headerTags=\
        ('refPos','refAllele','qryAllele','qryPos','nearestSnpDist','nearestEndDist','refLen','qryLen','refStrand','qryStrand','refId','qryId')
      #~ ofh.write('\t'.join(['refPos','refAllele','qryAllele','qryPos','nearestSnpDist','nearestEndDist','refStrand','qryStrand','refId','qryId'])+'\n')
      for line in ifh.readlines():
        #~ flds=line.strip().split()
        #~ for i in 0,3,4,5,6,7,8,9: flds[i]=int(flds[i])
        snpDict=dict(zip(headerTags,line.strip().split()))
        #~ if (snpDict['refId'],snpDict['qryId']) in ctgPairsSet and snpDict['refRepeatAlns']==0 and snpDict['qryRepeatAlns']==0:
        if (snpDict['refId'],snpDict['qryId']) in ctgPairsSet:
          ofh.write('\t'.join(\
            str(snpDict[tag]) \
            for tag in ('refId','refStrand','refPos','refAllele','qryId','qryStrand','qryPos','qryAllele','nearestSnpDist'))+'\n')
  return outputFp if (outputFp.is_file() and outputFp.stat().st_size>0) else None

if not hasattr(sys, 'real_prefix'):
  sys.exit("Beware: this script was not launched from within a virtualenv")

# # # # # enable using environment modules
if socket.gethostname()=='DellPrecisionTower7910-b21': exec(open('/usr/share/modules/init/python.py3','r').read())
else:  exec(open('/cm/local/apps/environment-modules/4.0.0/init/python.py','r').read())
#~ else:  exec(open('/users/jdrouau/usr/share/modules/init/python.py','r').read())


# # # # # main
# # initialization : create the arguments parser - parse the arguments - create the directories - define the global variables
now =                 timeStamp()
logGoodDict =         {offScore:(math.log10(1-10**((33-offScore)/10))) for offScore in range(34,256)}
parser =              fillParser(parser=argparse.ArgumentParser())
args =                sys.argv[1:] or multiLineInput()
params =              vars(parser.parse_args(args))
params =              initializeDirectories(params=params)
wlog('Directories initialized, logfile created')
wlog("Script content : "+"\n"+open(sys.argv[0],'r').read())
params =              buildMetaDataDict(params=params,log=True)
params =              buildInputFps(params=params,log=True)
wlog('parameters : '+str(params))
wlog('max number of cpus : '+str(multiprocessing.cpu_count()))
for key in params: exec(key+"=params[\'"+key+"\']")
os.environ['TMPDIR']=str(workDirPaths['tmp'])
tmpDir=str(workDirPaths['tmp'])
# # merge paired ends reads
#~ libNames=list(inputFqFpsDict.keys())

#~ minProbRange=floatRange(0.9,1,0.01) ; minLenRange=range(20,50,5)
#~ preprocParamsTest(inputFqFpsList=inputFqFpsList,minProbRange=minProbRange,minLenRange=minLenRange,projectTag=projectTag,cpu=cpu,log=True)
#~ sys.exit()
preprocFqFpsList=preprocess(inputFqFpsList=inputFqFpsList,projectTag=projectTag,minProb=minProb,minLen=minLen,cpu=cpu,log=True)
assembliesList=list()
assemblyIterationsNumber=3
for i in range(0,assemblyIterationsNumber):
  idbaAssemblyFpsDict=idbaAssembly(\
      readsFqFpsList=preprocFqFpsList,analysisPrefix=projectTag,kmersRange=range(52,125,8),\
      longFaFp=(None if i==0 else assembliesList[i-1]['spades']['contigs']),outputDirPath=workDirPaths['idba']/('round_'+str(i)),cpu=cpu,log=True)
  wlog('idbaAssemblyFpsDict round '+str(i)+' : '+str(idbaAssemblyFpsDict))
  spadesAssemblyFpsDict=spadesAssembly(\
      readsFqFpsList=preprocFqFpsList,kmersRange=range(51,130,4),trustedContigsFaFp=idbaAssemblyFpsDict['contigs'],\
      outputDirPath=workDirPaths['spades']/('round_'+str(i)),cpu=cpu,mem=mem,log=True)
  wlog('spadesAssemblyFpsDict round '+str(i)+' : '+str(spadesAssemblyFpsDict))
  assembliesList.append({'idba':idbaAssemblyFpsDict,'spades':spadesAssemblyFpsDict})
wlog('assembliesList : '+str(assembliesList))
assemblyFaFp=(workDirPaths['assembly']/(projectName+'_1kb.fa'))
assert assemblyFaFp==faRename2(\
  inputFaLines=faSizeFilter(inputFaFp=assembliesList[assemblyIterationsNumber-1]['spades']['scaffolds'],minSize=1000,log=False),\
  pattern='NODE_\d+',appendSize=True,outputFaFp=assemblyFaFp,log=False)
wlog('assemblyFaFp : '+str(assemblyFaFp))
#~ refAssemblyFaFp=workDirPaths['reference']/'Bethune_1kb.fa'
refAssemblyFaFpsDict=dict((ver,workDirPaths['reference']/('Bethune_'+ver+'_1kb.fa')) for ver in ('v1','v2'))
assert refAssemblyFaFpsDict['v1']==faRename2(\
  inputFaLines=faSizeFilter(inputFaFp=refGenomeFaFpsDict['v1'],minSize=1000,log=False),\
  pattern='^.*$',appendSize=True,outputFaFp=refAssemblyFaFpsDict['v1'],log=False)
wlog('refAssemblyFaFpsDict[v1] : '+str(refAssemblyFaFpsDict['v1']))
assert refAssemblyFaFpsDict['v2']==faRename2(\
  inputFaFp=refGenomeFaFpsDict['v2'],\
  pattern='CP0276[\d]{2}\.1',appendSize=False,outputFaFp=refAssemblyFaFpsDict['v2'],log=False)
wlog('refAssemblyFaFpsDict[v2] : '+str(refAssemblyFaFpsDict['v2']))
mummerDeltaFpsDict=dict() ; mummerDeltaFilteredFpsDict=dict() ; mummerCoordsFpsDict=dict() ; mummerSnpsFpsDict=dict()
goodCtgPairsFpsDict=dict() ; goodSnpsFpsDict=dict() ; goodAlnsFpsDict=dict()
for ver,refAssemblyFaFp in refAssemblyFaFpsDict.items():
  outputDirPath=(workDirPaths['alignment'])
  if not outputDirPath.is_dir(): outputDirPath.mkdir()
  mummerDeltaFpsDict[ver]=mummerNucmer(\
    referenceFaFp=refAssemblyFaFp,queryFaFp=assemblyFaFp,outputDirPath=outputDirPath,\
    prefix=projectName+'_vs_Bethune_'+ver,cpu=cpu,log=True)
  wlog('mummerDeltaFpsDict[ver] : '+str(mummerDeltaFpsDict[ver]))
  mummerDeltaFilteredFpsDict[ver]=mummerDeltaFilter(inputDeltaFp=mummerDeltaFpsDict[ver],outputDirPath=outputDirPath,prefix=None,cpu=cpu,log=True)
  wlog('mummerDeltaFilteredFpsDict[ver] : '+str(mummerDeltaFilteredFpsDict[ver]))
  mummerCoordsFpsDict[ver]=mummerShowCoords(deltaFp=mummerDeltaFilteredFpsDict[ver],outputDirPath=outputDirPath,cpu=cpu,log=True)
  wlog('mummerCoordsFpsDict[ver] : '+str(mummerCoordsFpsDict[ver]))
  mummerSnpsFpsDict[ver]=mummerShowSnps(deltaFp=mummerDeltaFilteredFpsDict[ver],outputDirPath=outputDirPath,cpu=cpu,log=True)
  wlog('mummerSnpsFpsDict[ver] : '+str(mummerSnpsFpsDict[ver]))
  goodCtgPairsFpsDict[ver]=coordsFileProcess(coordsFp=mummerCoordsFpsDict[ver],outputFp=None,log=True)
  wlog('goodCtgPairsFpsDict[ver] : '+str(goodCtgPairsFpsDict[ver]))
  goodSnpsFpsDict[ver]=snpsFileProcess(snpsFp=mummerSnpsFpsDict[ver],ctgPairsFp=goodCtgPairsFpsDict[ver],log=True)
  wlog('goodSnpsFpsDict[ver] : '+str(goodSnpsFpsDict[ver]))
  goodAlnsFpsDict[ver]=parallelProcess(\
    inputFpsList=[goodCtgPairsFpsDict[ver]],\
    outputFpsList=[mummerDeltaFilteredFpsDict[ver].parent/(mummerDeltaFilteredFpsDict[ver].prefix()+'.alns')],\
    outputZip=False,cpu=cpu,linesPerItem=1,chunkSize=100,fnName='mummerAlign2',fnArgs=(mummerDeltaFilteredFpsDict[ver],),log=True)
  wlog('goodAlnsFpsDict[ver] : '+str(goodAlnsFpsDict[ver]))
