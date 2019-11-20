# Submission
A python version of the grid submission


In order to submit:
* source your environment (as you would for running the analyzer)
* set ANALYSISDIR (export ANALYSISDIR=/path/to/your/analysisdir or setenv ANALYSISDIR /path/to/your/analysisdir)
* choose the files you want to run in SAMPLES_LIST.cfg (you can name it what you want)
* run ./remote.py SAMPLES_LIST.cfg
* add files ./add_root_files.py
* on first run you have to execute ./make_tester.sh once

One comment on the input for the SAMPLES_LIST.cfg:
* It will automatically update the samples that are not in the folder at the first run
* The sample name has to be something that fits:

        /eos/uscms/store/user/ra2tau/July72017/*/*NAME*
        /eos/uscms/store/user/ra2tau/jan2017tuple/*/*NAME*
        /eos/uscms/store/user/ra2tau/jruizalv/*/*NAME*

There is a bug in the xrd version of root https://root-forum.cern.ch/t/xrdcp-doesnt-work-after-changing-root-to-v6-06-08/22320
You must change to an other root version e.g. by setting up a newer CMSSW (e.g. CMSSW_9_2_X)

Before Submission you have to make a proxy:
voms-proxy-init -voms cms -rfc -valid 192:00

There are some options you might want to look at ./remote.py -h or ./add_root_files.py -h

```bash
./remote.py -h
Usage: remote.py [options] CONFIG_FILE

Options:
  -h, --help            show this help message and exit
  -C DIRECTORY, --configdir=DIRECTORY
                        Define the config directory. [default = PartDet]
  -c, --CR              Run with the CR flag. [default = False]
  -o DIRECTORY, --outputFolder=DIRECTORY
                        Define path for the output files [default = root://cms
                        eos.fnal.gov//store/user/YOURUSERNAME/REPLACEBYTAG/]
  -l, --local           run localy over the files [default = False]
  -f, --force           Force the output folder to be overwritten. [default =
                        False]
  --debug=LEVEL         Set the debug level. Allowed values: ERROR, WARNING,
                        INFO, DEBUG. [default = INFO]
  -t DIRECTORY, --Tag=DIRECTORY
                        Define a Tag for the output directory. [default =
                        output_DATE_short]
```

and:

```bash
./add_root_files.py -h
Usage: add_root_files.py [options]

Options:
  -h, --help            show this help message and exit
  -i FOLDER, --inputFolder=FOLDER
                        Merge all subfolders in these folders, which can be a
                        comma-separated list.[default: none]
  --debug=LEVEL         Set the debug level. Allowed values: ERROR, WARNING,
                        INFO, DEBUG. [default = INFO]
  -o OUTFOLDER, --output=OUTFOLDER
                        Set the output dir [default = output_DATE_long]
  -f, --force           If this option is specifed, all root files will be
                        remerged. [default = False]
  -c, --clean           If this option is specifed, the folders will be cleand
                        up. [default = False]
```
