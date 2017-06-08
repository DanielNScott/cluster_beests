from __future__ import division
import os
import sys
import cStringIO

sys.path.insert(0, '/gpfs_home/dscott3/code/beests-hack-versioned/stopsignal/stopsignal_wtf') # location of src

if (len(sys.argv) != 2):
	print("Please provide exactly one argument, the path to the analysis directory.")
	print("This path must contain 'sst_data.csv' and 'analysis.txt', the BEESTS data and config files.")
	exit(1)

path_analysisDir = sys.argv[1]
dataFile = sys.argv[1] + 'sst_data.csv'
path__analysisDescription = sys.argv[1] + 'analysis.txt'

#make dict for input arguments
vars = dict()
#load analysis.txt with input arguments and automatically set wd
with open(path__analysisDescription) as f:
    for line in f:
        eq_index = line.find(',')
        var_name = line[:eq_index].strip()
        var_name = var_name.replace('"', '').strip()
        number = line[eq_index + 1:].strip()
        number = number.replace('"', '').strip()
        vars[var_name] = number

samples	= int(vars["samples"])
burnIn 	= int(vars["burn-in"])
numberOfChains = int(vars["number of chains"])
thinning = int(vars["thinning"])
estimatesForAll = (vars["estimates for subjects or groups"]=="All")
summaryStatistics = (int(vars["summary statistics"]) != 0)
posteriorDistributions = (int(vars["posterior distributions"]) != 0)
mcmcChains = (int(vars["mcmc chains"]) != 0)
deviance = (int(vars["deviance"]) != 0)
posteriorPredictors = (int(vars["posterior predictors"]) != 0)
posteriorPredictorsSamples = int(vars["posterior predictor samples"])
numCores = int(vars["cpu cores"])
int_lower = int(vars["limits of integration lower"])
int_upper = int(vars["limits of integration upper"])
version = (int(vars["model trigger failure"]) == 0)

if  (version):
        import stopsignal
else:
        import stopsignal_wtf as stopsignal

import multiprocessing
from kabuki.analyze import gelman_rubin
from kabuki.utils import get_traces
from kabuki.utils import save_csv
from kabuki.utils import load_csv
import numpy
import math
from kabuki import Hierarchical

#		args << "-u"
#				<< _programDir.absoluteFilePath("stopsignal/run.py")
#				<< _dataFile
#				<< _analysisDir
#				<< _analysisDescriptionPath;

print(dataFile)
print(path_analysisDir)
print(path__analysisDescription)

data=load_csv(dataFile)
ncol=len(data.iloc[1])

actual_cores = multiprocessing.cpu_count()
if numCores < 1:
        print('Warning: Inputting 0 CPU core is silly. BEESTs will use the default number of ' + str(actual_cores-1) + ' cores.')
        numCores = actual_cores -1

if samples < 1:
        #print "The total number of MCMC samples must be greater than zero."
        sys.exit()

if samples <= burnIn:
        #print "The total number of MCMC samples must be higher than the number of burn-in samples."
        sys.exit()

if thinning < 1:
        #print "The thinning factor must be higher than 0."
        sys.exit()

if  ((samples-burnIn)/thinning) < 1:
        #print "No MCMC samples will be retained. Increse the number of retained samples or decrese the thinning factor."
        sys.exit()

if posteriorPredictors == True:
        if (((samples-burnIn)/thinning)*numberOfChains) < posteriorPredictorsSamples:
        #print "The number of posterior predictive samples cannot be higher than the number of retained MCMC samples."
                sys.exit()

if (int_lower >= int_upper):
        sys.exit()

def check_par(pars):
        for par in pars:
                if  (float(vars[par+" upper"]) <= float(vars[par+" lower"])):
                        sys.exit()
                if  ((float(vars[par+" start"]) < float(vars[par+" lower"])) or (float(vars[par+" start"]) > float(vars[par+" upper"]))):
                        sys.exit()
if  (version):
        print('BEESTS will fit the standard model.')
        check_par(pars=("go mu", "go sigma", "go tau", "stop mu", "stop sigma", "stop tau", "go mu sd",
              "go sigma sd", "go tau sd", "stop mu sd", "stop sigma sd", "stop tau sd"))
else:
        print('BEESTS will fit the trigger failure model.')
        #if (ncol==4):
        #        if (float(vars["stop pf lower"]) < 0 or float(vars["stop pf upper"]) > 1):
        #                sys.exit()
        #check_par(pars=("go mu", "go sigma", "go tau", "go shift", "stop mu", "stop sigma", "stop tau", "stop shift", "stop pf", "go mu sd",
        #      "go sigma sd", "go tau sd", "go shift sd", "stop mu sd", "stop sigma sd", "stop tau sd", "stop shift sd", "stop pf sd"))

#guess if rt data is in msec vs. sec
rts = data["rt"]
rts = rts[rts!=-999]
#print(min(rts))
#print(max(rts))
if ((min(rts)<80) | (max(rts)>2500)):
        print('The maximum and/or minimum RT in your dataset is unlikely. Are you sure that your data are in milliseconds?')


#--------------------------------------------------------------Start sampling
os.chdir(path_analysisDir)
models =[]
local_models = []

def run_model(i):
    ss = stopsignal.StopSignal(data)
    save_stdout = sys.stdout
    sys.stdout = cStringIO.StringIO()
    ss.find_starting_values()
    ss.sample(samples,burn=burnIn,thin=thinning,tune_throughout=False, db='pickle', dbname='remote_traces' + str(i+1) + '.db')
    sys.stdout = save_stdout
    return ss

if __name__ == "__main__":
        if actual_cores < numCores:
                if actual_cores==1:
                        print('Your system doesn\'t have ' + str(numCores) + ' cores. BEESTs will use 1 core.')
                        numCores = actual_cores
                else:
                      print('Your system doesn\'t have ' + str(numCores) + ' cores. BEESTs will use the default number of ' + str(actual_cores-1) + ' core(s).')
                      numCores = actual_cores-1
        else:
                print('Your system has ' + str(actual_cores) + ' cores. BEESTs will use ' + str(numCores) + ' core(s)')

        n_runs = math.ceil((numberOfChains/numCores))
        n_runs = numpy.array(n_runs).astype('int')
        num_remote = numberOfChains-n_runs
        num_remote = numpy.array(num_remote).astype('int')
        n_pool = numCores-1
        #print n_runs
        #print num_remote
        #print n_pool

        if num_remote>0:
                remote_model = []
                pool = multiprocessing.Pool(processes=n_pool)
                results = pool.map_async(run_model, range(num_remote))

        for i in range(n_runs):
                run_id = i+1
                beg = run_id + (i * n_pool)
                end = beg + n_pool
                if i == (n_runs-1):
                        end = numberOfChains
                ran_chains = range(beg,end+1)
                print('Running chain(s) ' + str(ran_chains))

                ss = stopsignal.StopSignal(data)
                print('\n Computing start values. It may take a few minutes. Or hours...')
                save_stdout = sys.stdout
                sys.stdout = cStringIO.StringIO()
                ss.find_starting_values()
                sys.stdout = save_stdout
                ss.sample(samples,burn=burnIn,thin=thinning,tune_throughout=False, db='pickle', dbname='local_traces' + str(i+1) + '.db')
                local_models.append(ss)
                if i == (n_runs-1):
                        print('Waiting for the other chains to finish...')

        if num_remote>0:
                models = results.get()

        for i in range(n_runs):
                models.append(local_models[i])

        #print(len(models))
        print "Finished sampling!"

        if numberOfChains > 1:
                Rhat = gelman_rubin(models)
                print('\n Gelman-Rubin Rhat diagnostic:')
                for key in Rhat.keys():
                        print((key, Rhat[key]))

        name_dataFile = dataFile.replace(".csv","")

        for i in range(numberOfChains):
                save_csv(get_traces(models[i]), 'parameters'+str(i+1)+'.csv', sep = ';')

        print "Posterior samples are saved to file."

        if deviance == True:
                for i in range(numberOfChains):
                        dev = models[i].mc.db.trace('deviance')()
                        numpy.savetxt('deviance'+str(i+1)+'.csv', dev, delimiter=";")
                print "Deviance values are saved to file"
#        if ss.is_group_model:
#          print "DIC: %f" % ss.mc.dic
