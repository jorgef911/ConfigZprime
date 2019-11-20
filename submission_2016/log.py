#!/bin/env python
import re
import curses
import time
from curses import wrapper
import sys
import datetime
import optparse
import os
import glob
import getpass
import subprocess
import logging
from collections import OrderedDict

id_format = " {:>40.40}"
formating= id_format + "  | {:>11}  | {:>13}  | {:>11}  |"
start_line = 2


class Subscreen():
    subscr = None
    maxx = 0
    maxy = 0
    total_jobs = 0
    prev_updown = 0
    highlight=0
    highlight_p = 0
    py = 0
    colors = {"COMPLETED": 3, "FAILED": 2, "ERROR": 2,  "RUNNING": 1, "COMPLETING": 1,"IDLE": 1, "PENDING": 5, "TIMEOUT": 2, "FINISHING": 1, "STARTING": 1}
    single_scrolling = False
    
    def __init__(self, y, x, total_jobs):
        self.maxx = x
        self.maxy = y
        self.subscr = curses.newpad(y, x)
        self.total_jobs = total_jobs

    def correct_highp(self):
        if self.highlight_p >= self.py + self.maxy:
            self.py = self.highlight_p
        elif self.highlight_p < self.py and self.highlight_p>=0:
            self.py = self.highlight_p
            
    def display_page(self, item_array):
        start=0
        position=0
        high_id=""
        if not self.single_scrolling:
            for i, line in enumerate(item_array):
                m=re.search("\A\s*((\w|-)+)  \|", item_array[i])
                if m is not None:
                    if self.highlight == start:
                        self.highlight_p = i
                        high_id = m.group(1)
                        break
                    start += 1
            self.correct_highp()
        else:
            if self.highlight_p >= len(item_array):
                self.highlight_p = len(item_array)-1
            m=re.search("\A(\s|-|\w)+\|\s+(\d+)  \|", item_array[self.highlight_p])
            if m is not None:
                high_id=m.group(1)
            m=re.search("\A\s*((\w|-)+)  \|", item_array[self.highlight_p])
            if m is not None:
                self.highlight += self.prev_updown

            
        for i in xrange(self.maxy):
            if i+self.py < len(item_array):
                color = curses.color_pair(1)
                m = re.search(r'\b([A-Z]+)\b', item_array[i+self.py])
                if m is not None: color = curses.color_pair(self.colors[m.group(0)])
                self.subscr.addstr(i, 0, item_array[i+self.py], color)
            else:
                self.subscr.addstr(i, 0, " " * 80)
        if len(high_id) <= 41:
            self.subscr.addstr(self.highlight_p-self.py, 41-len(high_id), high_id, curses.A_STANDOUT)
        else:
            self.subscr.addstr(self.highlight_p-self.py, 0, id_format.format(high_id), curses.A_STANDOUT)
        self.subscr.refresh(0,0,start_line, 0, self.maxy-start_line, self.maxx)

    def job_scroll(self, updown):
        self.single_scrolling = False
        if self.highlight+updown < 0:
            return
        elif self.highlight+updown >= self.total_jobs:
            return
        
        self.highlight += updown

    def single_scroll(self, updown):
        self.single_scrolling = True
        self.highlight_p += updown
        self.prev_updown = updown
        self.correct_highp()
        
    def set_job_scroll(self):
        self.single_scrolling = False
        
        
class Job_Holder():
    stdscr = None
    run_folder=None
    samples = OrderedDict()
    all_jobs = OrderedDict() #  [ state, time ]
    jobs_overview = OrderedDict() # [state, id, max_run]
    options=None
    
    def __init__(self, sample_list,options, stdscr = None):
        self.stdscr = stdscr
        self.options=options
        run_folder=options.runFolder
        for sample in sample_list:
            job_logs=glob.glob(os.path.join(run_folder,sample)+"/condor*.log")
            logging.info(os.path.join(run_folder,sample)+"/condor*.log")
            if len(job_logs)<1:
                sample_list.remove(sample)
                logging.info("wrong here")
                continue
            jobid=job_logs[0].split("_")[-1].split(".")[0]
            njobs=len(glob.glob(os.path.join(run_folder,sample)+"/run_*.sh"))
            self.samples[sample]=[jobid,njobs]
            
    def get_jobsize(self):
        return len(self.all_jobs)

    def refresh(self):
        command="condor_q %s"%(getpass.getuser())
        proc = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
        jobstati = proc.stdout.read()
        jobstati=jobstati.split("\n")
        for sample in self.samples:
            #submitted_time,runtime,frontend_status,size_task
            self.all_jobs[sample]=self.get_sample_info(sample,jobstati)
            overall_state = "COMPLETED"
            for task_id in self.all_jobs[sample]:
                tmp_state = self.all_jobs[sample][task_id][2]
                if tmp_state == "REMOVED":
                    overall_state = "ERROR"
                    break
                if tmp_state == "IDLE":
                    overall_state = "IDLE"
                    break
                if tmp_state != "COMPLETED":
                    overall_state = "RUNNING"
            max_time=max(self.all_jobs[sample][i][1] for i in self.all_jobs[sample])
            self.jobs_overview[sample]  = [overall_state, str(self.samples[sample][0]), str(max_time)]
            
    def overview_array(self, use_id, show_complete=True):
        return_array = []
        for filename, tmp_array in self.jobs_overview.iteritems():
            if not show_complete and tmp_array[0] == "COMPLETED":
                    continue
            job_num = filename
            if use_id:
                job_num = tmp_array[1]
            array_num = "----"
            if tmp_array[0] != "COMPLETED":
                array_num = str(len(self.all_jobs[filename]))
            write_str = formating.format(job_num, array_num, tmp_array[0], tmp_array[2])
            return_array.append(write_str)
        return return_array

    def all_array(self, use_id, show_complete):
        return_array = []
        for filename, job_array in self.all_jobs.iteritems():
            if self.jobs_overview[filename][0] == "COMPLETED":
                if not show_complete:
                    continue
                tmp_array = self.jobs_overview[filename]
                job_num = filename
                if use_id:
                    job_num = self.jobs_overview[filename][1]
                write_str = formating.format(job_num, "----", tmp_array[0], tmp_array[2])
                return_array.append(write_str)
                continue

            first = True
            for array_id, tmp_array in enumerate(job_array):
                array_id += 1
                write_str = ""
                if not show_complete and job_array[tmp_array][2] == "COMPLETED":
                    continue
                if first:
                    job_num = filename
                    if use_id:
                        job_num = self.jobs_overview[filename][1]
                    logging.info(tmp_array)
                    write_str = formating.format(job_num, array_id, job_array[tmp_array][2], job_array[tmp_array][1])
                    first = False
                else:
                    write_str = formating.format("", array_id, job_array[tmp_array][2], job_array[tmp_array][1])
                return_array.append(write_str)
        return return_array

    def get_total_jobs(self):
        return len(self.jobs_overview)
        
    def get_sample_info(self,sample,jobstati):
        sample_info={}
        job_id,ntasks=self.samples[sample]
        this_job_stati=filter(lambda x: job_id in x,jobstati)
        staus_translate={
        "H":"HOLD",
        "R":"RUNNING",
        "I":"IDLE",
        "C":"COMPLETED",
        "X":"REMOVED",
        "S":"SUSPENDED",
        ">":"FINISHING",
        "<":"STARTING",
        }
         #ID      OWNER            SUBMITTED     RUN_TIME ST PRI SIZE CMD               
        #270462.18  kpadeken       10/2  09:30   0+02:57:19 R  0   976.6 wrapper.sh /uscms_
        #270470.80  kpadeken       10/2  09:30   0+02:55:39 R  0   976.6 wrapper.sh /uscms_
        #270470.112 kpadeken       10/2  09:30   0+02:55:30 R  0   976.6 wrapper.sh /uscms_
        for stati in this_job_stati:
            tmp=stati.split()
            if len(tmp)>7:
                comb_id,user,date_subm,time_subm,run_time,status,priority,size_task =tmp[:8]
                task_id=comb_id.split(".")[-1]
                submitted_time=time.strptime("%d "%(time.localtime().tm_year)+date_subm+" "+time_subm, "%Y %m/%d %H:%M")
                days=run_time.split("+")[0]
                rest=run_time.split("+")[1]
                hour,minute,sec=rest.split(":")
                runtime=datetime.timedelta(days=int(days),hours=int(hour),minutes=int(minute),seconds=int(sec))
                try:
                    frontend_status=staus_translate[status]
                except:
                    logging.info(status)
                    logging.info(stati)
                sample_info[task_id]=[submitted_time,runtime,frontend_status,size_task]
        for i in range(ntasks):
            if "%d"%i not in sample_info:
                if sample in self.all_jobs and "%d"%i in self.all_jobs[sample]:
                    oldstatus=self.all_jobs[sample]["%d"%i]
                    sample_info["%d"%i]=[oldstatus[0],oldstatus[1],"COMPLETED",oldstatus[3]]
                else:
                    sample_info["%d"%i]=[time.localtime(),datetime.timedelta(0),"COMPLETED","0"]
        return sample_info
        
    def get_completion_info(self,sample,task_id):
        output_status="Done"
        total_runtime=datetime.timedelta(0)
        job_id,njobs=self.samples[sample]
        if(os.path.getsize(os.path.join(self.options.runFolder,sample,"err.%d_%s"%(task_id,job_id)))>0):
            output_status="Error"
        
        start_time=None
        end_time=None
        condor_log=open(os.path.join(self.options.runFolder,sample,"condor_%s.log"%(job_id)))
        
        
        
        
                
    
class Main_Program():
    max_line = 2
    expandall = False
    use_id = False
    show_complete = True
    prev_size = 0
    maxy = 0
    maxx = 0
    pad=None
    start_time= 0
    
    
    def __init__(self, sample_list,options):
        self.jobs = Job_Holder(sample_list,options)


    def header_footer(self,stdscr):
        title = formating.format("Job ID", "Array ID", "State", "Run Time")
        stdscr.addstr(0,0, title, curses.color_pair(4))
        stdscr.addstr(1,0, formating.format("","","",""), curses.A_UNDERLINE)
        options = ["q - quit", "r - refresh", "e - expand", "w - swap ID/Sample", "h - toggle complete", "j/k - move individual", "n/p - move job" ] 
        First = True
        footer_line = []
        tmp_footer=""
        spacer = 3
        for line in options:
            line = " "*spacer + line + " "*spacer
            if len(tmp_footer)+len(line) < self.maxx:
                if not First:
                    tmp_footer += "|"
                else: First = False
                tmp_footer += line
            else:
                footer_line.append(tmp_footer)
                tmp_footer = line
        footer_line.append(tmp_footer)
        for i, line in enumerate(footer_line):
            stdscr.addstr(self.maxy-len(footer_line)+i, 0, line, curses.color_pair(4))

        stdscr.refresh()    
        return len(footer_line)
        
    def display_pad(self):
        #self.jobs.refresh()
        if self.expandall:
            self.pad.display_page(self.jobs.all_array(self.use_id, self.show_complete))
        else:
            self.pad.display_page(self.jobs.overview_array(self.use_id, self.show_complete))
        
    def main(self, stdscr):
        curses.start_color()
        curses.init_pair(1, curses.COLOR_GREEN, curses.COLOR_BLACK)
        curses.init_pair(2, curses.COLOR_RED, curses.COLOR_BLACK)
        curses.init_pair(3, curses.COLOR_BLUE, curses.COLOR_BLACK)
        curses.init_pair(4, curses.COLOR_WHITE, curses.COLOR_BLACK)
        curses.init_pair(5, curses.COLOR_YELLOW, curses.COLOR_BLACK)
        self.maxy, self.maxx = stdscr.getmaxyx()
        stdscr.clear()

        
        lower_line_size = self.header_footer(stdscr)
        self.jobs.refresh()
        self.pad = Subscreen(self.maxy - lower_line_size + 1, self.maxx, self.jobs.get_total_jobs())
        self.display_pad()

        start_time = time.time()
        while True:
            if (time.time() - start_time) >= 5:
                self.jobs.refresh()
                self.display_pad()
                start_time = time.time()
                
            c = stdscr.getch()
            if c == ord('q'):
                break
            
            elif c == ord('n'):
                self.pad.job_scroll(1)
                self.display_pad()
                
            elif c == ord('p'):
                self.pad.job_scroll(-1)
                self.display_pad()
                
            elif c == ord('e'):
                self.expandall = not self.expandall
                if not self.expandall:
                    self.pad.set_job_scroll()
                self.display_pad()
                
            elif c == ord('w'):
                self.use_id = not self.use_id
                self.display_pad()

            elif c == ord('h'):
                self.show_complete = not self.show_complete
                self.display_pad()

            elif c == ord('j'):
                try:
                    self.pad.single_scroll(1)
                    self.display_pad()
                except:
                    self.pad.single_scroll(-1)
                
            elif c == ord('k'):
                try:
                    self.pad.single_scroll(-1)
                    self.display_pad()
                except:
                    self.pad.single_scroll(1)
            elif c == ord('r'):
                self.jobs.refresh()
                self.display_pad()
                    
                
                
def main():
    date_time = datetime.datetime.now()
    usage = '%prog [options]'
    parser = optparse.OptionParser( usage = usage )
    parser.add_option( '-i','--runFolder' , metavar = 'FOLDER', default=None,
                       help = 'Specify the folder where the jobs are [default: %default]' )
    parser.add_option( '--debug', metavar = 'LEVEL', default = 'INFO',
                       help= 'Set the debug level. Allowed values: ERROR, WARNING, INFO, DEBUG. [default = %default]' )
    ( options, args ) = parser.parse_args()
    
    logging.basicConfig(format='%(levelname)s:%(message)s',filename='log.log',level=logging.DEBUG)

    if (options.runFolder==None and not os.path.exists("submitted_samples.txt")):
        log.error( "You must give either a input file or keep the submitted_samples.txt" )
        sys.exit(3)
    if (options.runFolder==None and os.path.exists("submitted_samples.txt")):
        f=open("submitted_samples.txt","r")
        for line in f:
            if "outFolder:" in line:
                options.runFolder=line.replace("outFolder:","").strip()
                break
    if "/eos/uscms/" in options.runFolder:
        options.runFolder=options.runFolder.split(str(getpass.getuser())+"/")[-1]

    sample_list=glob.glob(options.runFolder+"/*")
    sample_list=[i.split("/")[-1] for i in sample_list]
    sample_list=filter(lambda x: ("exe" not in x),sample_list)

    main = Main_Program(sample_list,options)
        
    wrapper(main.main)
    
if __name__ == '__main__':
    main()

