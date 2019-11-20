#!/bin/env python

#__________initial imports__________
import binConfig  #s
import checkEnvironment  #s
from datetime import datetime  #s
import optparse,os,time,pickle,subprocess,shutil,sys,getpass,re  #s
import logging  #s
import ROOT  #s
from condor_submit import checkAndRenewVomsProxy  #d
from collections import OrderedDict  #d 
from commands import getoutput  #d
log = logging.getLogger( 'remote' )  #s
#__________end initial imports_______

#__________store generators__________
generators = OrderedDict()
generators['madgraph'] = 'MG'
generators['powheg'] = 'PH'
generators['herwig6'] = 'HW'
generators['herwigpp'] = 'HP'
generators['herwig'] = 'HW'
generators['sherpa'] = 'SP'
generators['amcatnlo'] = 'AM'
generators['alpgen'] = 'AG'
generators['calchep'] = 'CA'
generators['comphep'] = 'CO'
generators['lpair'] = 'LP'
generators['pythia6'] = 'P6'
generators['pythia8'] = 'P8'
generators['pythia'] = 'PY'
generators['gg2ww'] = 'GG'
generators['gg2zz'] = 'GG'
generators['gg2vv'] = 'GG'
generators['JHUGen'] = 'JG'
generators['blackmax'] = 'BM'
generators['unknown'] = ''
#__________end store generators__________

# list of generators used for hadronization on top of another generator (will be removed from name)
showers = [ 'pythia8', 'pythia6', 'pythia', 'herwigpp']

# list of tags which will be removed from name (case insensitive)
blacklist = ['13tev',
             'madspin',
             'FXFX',
             'MLM',
             'NNPDF30',
             'TuneCUEP8M1',
             'TuneCUETP8M1',
             'TuneCUETP8M2T4']

#__________establish sample name__________
def parse_name(dataset, options):
    # format of datasetpath: /.../.../...
    # first part contains name + additional tags ( cme, tune, .. )
    # second part has additional information ( campaign, extention sample? ,... )
    # third part contains sample 'format' (AOD, MINIAOD, ...)
    dataset_split = dataset.split('/')  #splitting the name
    ds_pt1 = dataset_split[1]  #first part of the name  
    ds_pt2 = dataset_split[2]  #second part of the name
    ds_pt3 = "MC" if "SIM" in dataset_split[3] else "data"  #third part of the name
    if ds_pt3 =="data":
        ds_pt3+="_"+ds_pt2.split("-")[0][-1]
    for generator in generators.keys():
        # subn() performs sub(), but returns tuple (new_string, number_of_subs_made)
        # using (?i) at the beginning of a regular expression makes it case insensitive
        ( ds_pt1, n ) = re.subn( r'(?i)[_-]' + generator, '', ds_pt1 )
        if n > 0:
            _generator = generator
            for shower in showers:
                ds_pt1 = re.sub( r'(?i)[_-]' + shower, '', ds_pt1 )
            break
        else:
            _generator = 'unknown'
    #for item in blacklist:  #I don't think that we want the blacklist b/c file names are different in the newer productions.
        #ds_pt1 = re.sub( r'(?i)[_-]*' + item, '', ds_pt1 )
    match = re.search('ext\d\d*',ds_pt2)
    #match = re.search('ext\d\d*',ds_pt2),options.cme
    #name = ds_pt1 + "_" + options.cme + "TeV_" + match.group() + "_" + options.postfix + generators[_generator]+"_"+ds_pt3
    if match:
        name = ds_pt1 + "_" + "TeV_" + match.group() + "_" + generators[_generator]+"_"+ds_pt3  
    else:
        name = ds_pt1 + "_" + "TeV_" + generators[_generator]+"_"+ds_pt3  #This may need to change for data.  Let's check.  
    return name
#__________end establish sample name__________

def bins(file_list, bin_size):

    binning_list = []

    maxsize=0
    current_bin=[]
    for i, file in enumerate(file_list):
        current_bin.append(file)
        maxsize+=file_list[file]
        # one file per bin
        if maxsize>bin_size or (maxsize+file_list[file]*0.5)>bin_size:
            #print(maxsize,current_bin)
            binning_list.append(current_bin)
            current_bin=[]
            maxsize=0
        elif i==len(file_list)-1:
            binning_list.append(current_bin)
            
    return binning_list

#__________be able to read the files via xrootd__________
def getDatasetFileList(sample):
 command = 'dasgoclient --query="file dataset=%s | grep file.name, file.size"' % (sample)
 output = getoutput(command)
 fileList = {}
 for line in output.split(os.linesep):
     file,size=line.split()
     fileList['root://cmsxrootd.fnal.gov//'+file]=int(size)
     #fileList['root://cms-xrd-global.cern.ch/'+file]=int(size)
 return fileList
#__________end be able to read the files via xrootd__________

#__________begin list of samples__________
def getFilesfromFile(cfgFile, options):
    sampleList={}
    file = open(cfgFile,'r')
    
    for line in file:
        if line[0]=="#" or len(line.split())==0:
            continue
        sample=line.strip()  #Skip the guys who are commented.
        
        file_lists=bins(getDatasetFileList(sample),6400)  #size in bytes 3GB
        sampleList[parse_name(sample,options)]=file_lists
    return sampleList
#__________end list of samples__________

#__________setting config files and proper flags__________
def makeExe(options,inputfiles,outputfile,sample):
    from string import Template
    exe="""
#!/bin/bash -e
sleep $[ ( $RANDOM % 30 ) ]
date
cd ${_CONDOR_SCRATCH_DIR}
tar -xvzf exe.tar.gz
ls
isData=$ISDATA
echo $isData
if [ ! -z $isData ]
then
    echo "switch data to true"
    sed -r -i -e 's/(isData\s+)(0|false)/isData true/' -e 's/(CalculatePUS[a-z]+\s+)(1|true)/CalculatePUSystematics false/' \
    $CONFIGDIR/Run_info.in
else
    echo "switch data to false"
    sed -r -i -e 's/(isData\s+)(1|true)/isData false/' -e 's/(CalculatePUS[a-z]+\s+)(0|false)/CalculatePUSystematics true/' \
    $CONFIGDIR/Run_info.in
fi
./Analyzer -in $INPUTFILES -out $OUPUTFILE -C $CONFIGDIR $CONTOLLREGION
xrdcp -sf $_CONDOR_SCRATCH_DIR/$OUPUTFILE $OUTPUTFOLDER$SAMPLE/$OUPUTFILE
"""

    for fileDir in binConfig.cpFiles:
        exe+="rm -r "+fileDir+" \n"
    exe+="rm -r $CONFIGDIR \n"
    exe+="rm -r *.root \n"
    exe+="rm -r *.tar.gz \n"
    CR=""
    if options.CR:
        CR="-CR"
    #isdata= "RunII" in inputfiles[0] or "Tune" in inputfiles[0]
    isdata= not "_Run20" in inputfiles[0]
    d = dict(
            CONFIGDIR=options.configdir,
            INPUTFILES=" ".join(inputfiles),
            OUPUTFILE=outputfile,
            OUTPUTFOLDER=options.outputFolder,
            SAMPLE=sample,
            CONTOLLREGION=CR,
            ISDATA="" if isdata else "false",
        )
    exe=Template(exe).safe_substitute(d)
    exeFile=open("run_"+outputfile.replace(".root","")+".sh","w+")
    exeFile.write(exe)
    exeFile.close()
#__________end setting config files and proper flags__________

def main():

    date_time = datetime.now()
    usage = '%prog [options] CONFIG_FILE'
    parser = optparse.OptionParser( usage = usage )
    parser.add_option( '-C', '--configdir', default = "PartDet", metavar = 'DIRECTORY',
                            help = 'Define the config directory. [default = %default]')
    parser.add_option( '-c', '--CR', action = 'store_true', default = False,
                            help = 'Run with the CR flag. [default = %default]')
    parser.add_option( '-o', '--outputFolder', default = "root://cmseos.fnal.gov//store/user/%s/REPLACEBYTAG/"%(getpass.getuser()), metavar = 'DIRECTORY',
                            help = 'Define path for the output files [default = %default]')
    parser.add_option( '-l', '--local',action = 'store_true', default = False,
                            help = 'run localy over the files [default = %default]')
    parser.add_option( '-f', '--force',action = 'store_true', default = False,
                            help = 'Force the output folder to be overwritten. [default = %default]')
    parser.add_option( '--debug', metavar = 'LEVEL', default = 'INFO',
                       help= 'Set the debug level. Allowed values: ERROR, WARNING, INFO, DEBUG. [default = %default]' )
    parser.add_option( '--filesFromACCRE', action = 'store_true', default = False,
                       help= 'Use the files from ACCRE [default = %default]' )
    parser.add_option( '-t', '--Tag', default = "run_%s_%s_%s_%s_%s"%(date_time.year,
                                                                        date_time.month,
                                                                        date_time.day,
                                                                        date_time.hour,
                                                                        date_time.minute,
                                                                        ), metavar = 'DIRECTORY',
                        help = 'Define a Tag for the output directory. [default = %default]' )
        

    ( options, args ) = parser.parse_args()
    if len( args ) != 1:
        parser.error( 'Exactly one CONFIG_FILE required!' )
    options.outputFolder=options.outputFolder.replace("REPLACEBYTAG",options.Tag)
    
    print("You may enter your grid password here. Do not enter anything to use the available proxy.")
    passphrase = getpass.getpass()
    if passphrase=="":
        passphrase = None
    else:
        checkAndRenewVomsProxy(passphrase=passphrase)



    format = '%(levelname)s from %(name)s at %(asctime)s: %(message)s'
    date = '%F %H:%M:%S'
    logging.basicConfig( level = logging._levelNames[ options.debug ], format = format, datefmt = date )
    log.info("Welcome to the wonders of color!")

    try:
       cmssw_version, cmssw_base, scram_arch = checkEnvironment.checkEnvironment()
    except EnvironmentError as err:
        log.error( err )
        log.info( 'Exiting...' )
        sys.exit( err.errno )


    if os.path.exists(options.outputFolder.replace("root://cmseos.fnal.gov/","/eos/uscms/")) and not options.force:
        log.error("The outpath "+options.outputFolder+" already exists pick a new one or use --force")
        sys.exit(3)
    elif options.force and os.path.exists(options.outputFolder.replace("root://cmseos.fnal.gov/","/eos/uscms/")):
        shutil.rmtree(options.outputFolder.replace("root://cmseos.fnal.gov/","/eos/uscms/"))
    os.makedirs(options.outputFolder.replace("root://cmseos.fnal.gov/","/eos/uscms/"))
    
    cfgFile = args[ 0 ]
    sampleList=getFilesfromFile(cfgFile,options)
    

    thisdir=os.getcwd()
    
    
    exepath=os.path.join(options.Tag,"exe")
    if os.path.exists(exepath):
        shutil.rmtree(exepath)
    os.makedirs(exepath)
    anadir=binConfig.PathtoExecutable.replace("/uscms/home/","/uscms_data/d3/").replace("nobackup/","")
    for fileDir in binConfig.cpFiles:
        if os.path.isdir(os.path.join(anadir,fileDir)):
            shutil.copytree(os.path.join(anadir,fileDir),os.path.join(exepath,fileDir))
        else:
            shutil.copy(os.path.join(anadir,fileDir),os.path.join(exepath,fileDir))
    shutil.copytree(os.path.join(anadir,options.configdir),os.path.join(exepath,options.configdir))
    
    os.chdir(exepath)
    command="tar czf exe.tar.gz *"
    subprocess.call(command, shell=True)
    os.chdir(thisdir)
    pathtozip=os.path.join(os.path.abspath(exepath),"exe.tar.gz")

    n_jobs=0
    for sample in sampleList:
        n_jobs+=len(sampleList[sample])
    print(("There will be %d jobs in total"%n_jobs))  #total number of jobs

    sbumittedjobs=0
    for sample in sampleList:
        os.chdir(thisdir)
        if os.path.exists(os.path.join(options.Tag,sample)) and not options.force:
            log.error("The samplepath "+os.path.join(options.Tag,sample)+" already exists use the --force")
            sys.exit(3)
        elif options.force and os.path.exists(os.path.join(options.Tag,sample)):
            shutil.rmtree(os.path.join(options.Tag,sample))
        os.makedirs(os.path.join(options.Tag,sample))
        os.chdir(os.path.join(options.Tag,sample))
        ## I know not the way we want to trigger stuff
        wrapper="""#!/bin/bash -e
export SCRAM_ARCH=slc7_amd64_gcc630
ls -lrth
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 project CMSSW CMSSW_10_1_9`
cd CMSSW_10_1_9
ls -lrth
eval `scramv1 runtime -sh`
cp ../$@ run.sh
chmod u+x run.sh
./run.sh
        """
        condor_jdl="executable      = "+os.path.join(os.getcwd(),"wrapper.sh")+"\n"
        condor_jdl+="""
universe        = vanilla
Error           = err.$(Process)_$(Cluster)
Output          = out.$(Process)_$(Cluster)
Log             = condor_$(Cluster).log
transfer_input_files = %s
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
request_memory  = 0.5 GB
Notification    = NEVER
x509userproxy = $ENV(X509_USER_PROXY)
"""%(", ".join([pathtozip]+[os.path.join(os.getcwd(),"run_%s_%d.sh"%(sample,i)) for i,binned in enumerate(sampleList[sample]) ]) )
        f=open("wrapper.sh","w")
        f.write(wrapper)
        f.close()
        
        for i,binned in enumerate(sampleList[sample]):
            makeExe(options,binned,"%s_%d.root"%(sample,i),sample)
            condor_jdl+="arguments = %s \n"%("run_%s_%d.sh"%(sample,i))
            condor_jdl+="queue\n"
        condor_jdl+="\n"
        f=open("condor.jdl","w")
        f.write(condor_jdl)
        f.close()
        log.info("Submitting sample %s"%sample)
        command="condor_submit condor.jdl"
        log.debug(command)
        subprocess.call(command, shell=True)
    
    os.chdir(thisdir)
    legacy_file=open("submitted_samples.txt","w")
    legacy_file.write("outFolder:%s\n"%options.outputFolder.replace("root://cmseos.fnal.gov/","/eos/uscms/"))
    for sample in sampleList:
        legacy_file.write("%s\n"%(sample))
    legacy_file.close()
    log.info("Thanks for zapping in, bye bye")
    log.info("The out files will be in "+options.outputFolder)
    log.info("Check the status with condor_q %s"%(getpass.getuser()))
    log.info("When finished run ./add_root_files.py")

if __name__ == '__main__':
    main()
