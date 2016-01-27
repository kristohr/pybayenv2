import sys, os
import commands
import time
import threading
import multiprocessing
import random


# thread class to run the Bayenv test phase in multiple threads
class RunInThread(threading.Thread):
    
    def __init__(self, cmd, thread_id, testdata):
        threading.Thread.__init__(self)
        self.cmd = cmd
        self.thread_id = thread_id
        self.testdata = testdata
        #self.queue = queue

    def run(self):
        # execute the command, queue the result
        #(status, output) = commands.getstatusoutput(self.cmd)
        #elf.queue.put((self.cmd, output, status))
        #result_file = "bf_environ." + env_file
        errors = open("errors.out", 'w')
        outfile = "bf_results_t" + str(self.thread_id) + ".txt"
        k = 0
        dataset = open(self.testdata, 'r')
        thread_id = str(self.thread_id)
        cmd = self.cmd
        while dataset:
            k += 1
            line1 = dataset.readline()
            if line1 == '':
                break
            l = line1.split("\t")
            marker_file = l.pop(0) + "-" + str(self.thread_id)
            FILE = open(marker_file, 'w')
            line1 = "\t".join(l)
            line2 = dataset.readline()
            FILE.write(line1 + line2)
            FILE.close()
            print "BAYENV: thread " + thread_id + " is processing " + marker_file + " (" + str(k) + ")...",
            start_test = time.time()
            sys.stdout.flush()
            failure, output = commands.getstatusoutput(cmd % (marker_file, outfile))
            elapsed = (time.time() - start_test)
            print "done. %f sec to complete" % elapsed
            #total = ((elapsed * 33000)/60)/60
            #print "Total: " + str(total)
        
            if failure:
                print output
                errors.write(output)
                #sys.exit(1)
                error =  "Could not test locus: " + marker_file +"\n"\
                      + line1 + line2
                errors.write(error)
                os.remove(marker_file)
                os.remove(marker_file + ".freqs")
                continue
            os.remove(marker_file)
            os.remove(marker_file + ".freqs")
        errors.close()
        dataset.close()
        #new_result_file = "iter" + str(i) + "_" + result_file
        #os.rename(result_file, new_result_file) #Rename result file.
     
# Process class to run the Bayenv test phase in paralell
class RunInProcess(multiprocessing.Process):
    
    def __init__(self, cmd, thread_id, testdata, testsize):
        multiprocessing.Process.__init__(self)
        self.cmd = cmd
        self.thread_id = thread_id
        self.testdata = testdata
        self.testsize = testsize

    def run(self):
        errors = open("errors.out", 'w')
        outfile = "results/bf_results_t" + str(self.thread_id)
        k = 0
        dataset = open(self.testdata, 'r')
        thread_id = str(self.thread_id)
        cmd = self.cmd
        data_size = self.testsize
        while dataset:
            k += 1
            line1 = dataset.readline()
            if line1 == '':
                break
            l = line1.split("\t")
            marker_file = l.pop(0) + "-" + str(self.thread_id)
            FILE = open(marker_file, 'w')
            line1 = "\t".join(l)
            line2 = dataset.readline()
            FILE.write(line1 + line2)
            FILE.close()
            print "BAYENV: process " + thread_id + " is processing " + marker_file + " (" + str(k) + ")...",
            start_test = time.time()
            sys.stdout.flush()
            failure, output = commands.getstatusoutput(cmd % (marker_file, outfile))
            elapsed = (time.time() - start_test)
            remaining = ((data_size-k)*elapsed)/60
            print "done. %f sec to complete. Estimated time remaining: %f minutes" % (elapsed, remaining)
        
            if failure:
                print output
                errors.write(output)
                error =  "Could not test locus: " + marker_file +"\n"\
                      + line1 + line2
                errors.write(error)
                os.remove(marker_file)
                os.remove(marker_file + ".freqs")
                continue
            os.remove(marker_file)
            os.remove(marker_file + ".freqs")
        errors.close()
        dataset.close()

#Running BAYENV sequentially
def test_run_seq(cmd, thread_id, testdata):

    thread_id = str(thread_id)
    errors = open("errors.out", 'w')
    outfile = "bf_results_t" + thread_id + ".txt"
    k = 0
    dataset = open(testdata, 'r')
    while dataset:
        k += 1
        line1 = dataset.readline()
        if line1 == '':
            break
        l = line1.split("\t")
        marker_file = l.pop(0) + "-" + thread_id
        FILE = open(marker_file, 'w')
        line1 = "\t".join(l)
        line2 = dataset.readline()
        FILE.write(line1 + line2)
        FILE.close()
        print "BAYENV: sequence " + thread_id + " is processing " + marker_file + " (" + str(k) + ")...",
        start_test = time.time()
        sys.stdout.flush()
        failure, output = commands.getstatusoutput(cmd % (marker_file, outfile))
        elapsed = (time.time() - start_test)
        print "done. %f sec to complete" % elapsed
            #total = ((elapsed * 33000)/60)/60
            #print "Total: " + str(total)
        
        if failure:
            print output
            errors.write(output)
                #sys.exit(1)
            error =  "Could not test locus: " + marker_file +"\n"\
                + line1 + line2
            errors.write(error)
            os.remove(marker_file)
            os.remove(marker_file + ".freqs")
            continue
        os.remove(marker_file)
        os.remove(marker_file + ".freqs")
    errors.close()
    dataset.close()
        #new_result_file = "iter" + str(i) + "_" + result_file
        #os.rename(result_file, new_result_file) #Rename result file.
  

############# Running BAYENV sequentially ###############
def test_all_snps_seq(testdata, env_file, no_env_var):
        
    cmds = []
    num_runs = 4
    rand_seeds = [int(random.uniform(1,99999)) for i in range(1, num_runs+1)]

    for i, rand_seed in enumerate(rand_seeds):
        cmds.append("bayenv2 -i %s -m " + "covar_iter" + str(i+1) + ".txt" + \
                    " -e " + env_file + \
                    " -p 3 " + \
                    "-k 20000 " + \
                    " -n " + str(no_env_var) + \
                    " -t -c -r " + str(rand_seed) + " -o %s")

    start = time.time()
    for i in range(len(cmds)):
        test_run_seq(cmds[i], i, testdata)
    print "Elapsed time: %s" % (time.time()-start)
    exit(0)



############# Running BAYENV with threads ###############
def test_all_snps_multit(testdata, env_file, no_env_var):
    
    threads = []
    cmds = []
    num_runs = 4
    rand_seeds = [int(random.uniform(1,99999)) for i in range(1, num_runs+1)]

    for i, rand_seed in enumerate(rand_seeds):
        cmds.append("bayenv2 -i %s -m " + "covar_iter" + str(i+1) + ".txt" + \
                    " -e " + env_file + \
                    " -p 3 " + \
                    "-k 10000 " + \
                    " -n " + str(no_env_var) + \
                    " -t -c -r " + str(rand_seed) + " -o %s")
    start = time.time()
    for i in range(len(cmds)):
        thread = RunInThread(cmds[i], i, testdata)
        threads.append(thread)
        thread.start()
    
    for t in threads:
        t.join()
    print "Elapsed time: %s" % (time.time()-start)

    #Call the second difference method
    #compute_results(num_runs, 1, 1.5) 

    #exit(0)

############# Running BAYENV with multiprocessing ###############
def test_all_snps_multip(testdata, cmds, testsize):
    
    procs = []

    start = time.time()
    for i in range(len(cmds)):
        proc = RunInProcess(cmds[i], i, testdata, testsize)
        procs.append(proc)
        proc.start()
    
    for p in procs:
        p.join()
    print "Elapsed time: %s" % (time.time()-start)
    #Call the second difference method
    #compute_results(len(cmds), 1, 0.5)
    #exit(0)

############# Running Bayenv in test modus ##############
def test_all_snps(testdata, env_file, n_iter, no_env_var):

    result_file = "bf_environ." + env_file
    errors = open("errors.out", 'w')

    for i in range(1, int(n_iter)+1):
        
        dataset = open(testdata, 'r')
        covar_file = "covar_iter" + str(i) + ".txt"
        #cmd = "bayenv 1 3 -47123 30000 %s " + env_file + " " + \
        cmd = "bayenv 1 3 -74123 50000 %s " + env_file + " " + \
              str(no_env_var) + " " + covar_file

        cmd = "bayenv2 -i %s -m " + covar_file + " -e " + env_file + " -p 3 " + \
              "-k 50000 -n " + str(no_env_var) + " -t -c -r 44332" 
        
        k = 0
        while dataset:
            k += 1
            line1 = dataset.readline()
            if line1 == '':
                break
            l = line1.split("\t")
            marker_file = l.pop(0)
            FILE = open(marker_file, 'w')
            line1 = "\t".join(l)
            line2 = dataset.readline()
            FILE.write(line1 + line2)
            FILE.close()
            print "BAYENV processing " + marker_file + " (" + str(k) + ")...",
            start_test = time.time()
            sys.stdout.flush()
            failure, output = commands.getstatusoutput(cmd % marker_file)
            elapsed = (time.time() - start_test)
            print "done. %f sec to complete" % elapsed
            #total = ((elapsed * 33000)/60)/60
            #print "Total: " + str(total)

            if failure:
                print output
                errors.write(output)
                #sys.exit(1)
                error =  "Could not test locus: " + marker_file +"\n"\
                      + line1 + line2
                errors.write(error)
                os.remove(marker_file)
                continue
            os.remove(marker_file)
        errors.close()
        dataset.close()
        new_result_file = "iter" + str(i) + "_" + result_file
        os.rename(result_file, new_result_file) #Rename result file.
        #run_r_cmd(new_result_file, i)

        #make_pdf_report(i)

    

#Calculating covariance matrices every 5000 iterations - Bayenv 1.0
def compute_null_hyp(num_pops, rand_seed, iterations, snpfile):
    
    time_usage = open("t_usage.out", 'w')
    cmd = "bayenv 0 " + num_pops +" "+ rand_seed +" "+ iterations \
          +" "+ snpfile + " > covars.txt"

    #print cmd
    
    print "BAYENV calculating the covariance matrices...",
    start_test = time.time()
    sys.stdout.flush()
    failure, output = commands.getstatusoutput(cmd)
    elapsed = (time.time() - start_test)
    usage = "done. %f sec to complete" % elapsed
    print usage
    time_usage.write(usage)
    time_usage.close()

#Calculating covariance matrices every 500 iterations - Bayenv 2.0
def compute_null_model_bayenv2(num_pops, iterations, snpfile):

    print "#############################################"
    
    time_usage = open("t_usage.out", 'w')

    rand_seed = str(int(random.uniform(1,99999)))
    print "Random seed = " + rand_seed

    cmd_str = open("covar-cmd.txt", "wb")

    cmd = "bayenv2 -i " + snpfile + " -p " + str(num_pops) + " -k " + str(iterations) \
        + " -r " + str(rand_seed) + " > covars.txt"

    print cmd
    cmd_str.write(cmd)
    cmd_str.close()
   
    print "BAYENV calculating the covariance matrices...",
    start_test = time.time()
    sys.stdout.flush()
    failure, output = commands.getstatusoutput(cmd)
    elapsed = (time.time() - start_test)
    usage = "done. %f sec to complete" % elapsed
    print usage
    time_usage.write(usage)
    time_usage.close()
