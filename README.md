# rmcmc
Rossiter-McLaughlin Monte Carlo simulation

### Installation
`cd /path/to`

`git clone https://github.com/emilknudstrup/rmcmc.git`

`cd /path/to/rmcmc`

`python -m pip install .`

`pip install -r requirements.txt`


### Example: Manually create target
```python
import rmcmc

## Instantiate target
target = rmcmc.target.Target(per=18.71205, T0=2458729.681, rp=0.0690, 
			     a=24.3, b=0.14, ecc=0.0, w=90.0, vsini=5.64,
			     name='TOI-1456')

## Instantiate the Simulator with 300 draws, precision of 1 m/s, and 3 CPUs (multiprocessing)
run = rmcmc.fakit.Simulator(target,ndraws=300,precision=1,nproc=3)
## Simulate some data with an exposure time of 180 s
run.dataSimulator(exposure=180)

## Save results or not
save = 0
##
lams = [0.0,10.] # lambda values to test
run.MC(lams,writeBinary=save)
if save:
	results = rmcmc.fakit.readBinary('TOI-1456.pkl')
	rmcmc.makit.plotKDE(results)
else:
	rmcmc.makit.plotKDE(run.results)


```
### Example: Target from a CSV file 
CSV file should have a format like this with column names corresponding to the parameter name in the Target class

| name     | host          | per    | T0           | a    | b    | vsini | rp    |
| -------- | --------------| ------ | ------------ | ---- | ---- | ----- | ----- |
| TOI-150b | TYC 9191-519-1| 5.8574 | 2458326.2773 | 9.92 | 0.33 | 7.96  | 0.083 |
| ...      | ...           | ...    | ...          | ...  | ...  | ...   | ...   |
| K2-232b  | HD 286123     | 11.168 | 2458326.2773 | 17.1 | 0.59 | 5.15  | 0.091 |

```python
import rmcmc

## Path to CSV file
file = '/home/emil/Desktop/PhD/targets/warm.csv'
## Instantiate target
target = rmcmc.target.Target()
## Grab and set parameters for your planet name in CSV file
target.fromCSV(file,'TOI-677b',skiprows=2)

## Instantiate the Simulator with 400 draws, precision of 5 m/s, and 3 CPUs (multiprocessing)
run = rmcmc.fakit.Simulator(target,ndraws=400,precision=5,nproc=3)

## Simulate some data with an exposure time of 360 s
run.dataSimulator(exposure=360)

## Save results or not
save = 0
## Run the MC
lams = [0.,10.]#Lambdas to test
run.MC(lams,writeBinary=save)
if save:
	results = rmcmc.fakit.readBinary('TOI-677b.pkl')
	rmcmc.makit.plotKDE(results)
else:
	rmcmc.makit.plotKDE(run.results,usetex=True,path='./TOI-677b')#path to plot, name './TOI-677b_MC_KDE.png'


```
This should look something like this:
![TOI-677b_MC_KDE](https://user-images.githubusercontent.com/50403597/191103313-e48a7c02-452c-42f2-a3d1-a510d3cb05ca.png)

### Example: Target and exposure time/precision from a CSV file 


| name     | exp   | precision | jitter |
| -------- | ----- | --------  | ------ |
| WASP-24b | 450   | 1.25      | 3.0    |
| ...      | ...   | ...       | ...    |
| WASP-54b | 450   | 0.83      | 2.0    |

```python

import pandas as pd
import rmcmc
import os

## Path to CSV file
file = '/home/emil/Desktop/PhD/targets/1deg.csv'
## Instantiate target
target = rmcmc.target.Target()
## Grab and set parameters for your planet name in CSV file
name = 'WASP-24b'
path = './mc/'
try:
	os.mkdir(path)
except FileExistsError:
	pass

target.fromCSV(file,name,skiprows=2)

## Read file with target, exposure time, precision, and jitter
obsfile = pd.read_csv('targets_exposure.csv',skiprows=1)
df = obsfile[obsfile['name'] == name]
prec, jitt = df['precision'][0], df['jitter'][0]
exp = df['exp'][0]
## Instantiate the Simulator with 400 draws, precision of prec m/s, jitter of jitt m/s, and 3 CPUs (multiprocessing)
run = rmcmc.fakit.Simulator(target,ndraws=400,precision=prec,jitter=jitt,nproc=3)

## Simulate some data with an exposure time of exp s
run.dataSimulator(exposure=exp)

## Save results or not
save = 0
if save:
	## Lambdas to test
	lams = [0.,10.]
	## Run the MC
	run.MC(lams,writeBinary=save,path=path)
	rmcmc.makit.plotKDE(run.results)
else:
	results = rmcmc.fakit.readBinary(path+name+'.pkl')
	rmcmc.makit.plotKDE(results,usetex=True,path=path+name)#path to plot 'pathname_MC_KDE.png'


```