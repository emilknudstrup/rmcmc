# rmcmc
Rossiter-McLaughlin Monte Carlo simulation

### Installation
`cd /path/to`

`git clone https://github.com/emilknudstrup/rmcmc.git`

`cd /path/to/rmcmc`

`python -m pip install .`

`pip install -r requirements.txt`


### Example: Manual create target
```python
import rmcmc

## Instantiate target
target = rmcmc.target.Target(per=18.71205, T0=2458729.681, rp=0.0690 , 
							 a=24.3, b=0.14, ecc=0.0, w=90.0,vsini=5.64,
							 name='TOI-1456')

## Instantiate the Simulator with 330 draws, precision of 1 m/s, and 3 CPUs (multiprocessing)
run = rmcmc.fakit.Simulator(target,ndraws=300,precision=1,nproc=3)
## Simulate some data with an exposure time of 180 s
run.dataSimulator(exposure=180)

## Save results or not
save = 0
## 
run.MC([0.0,10],writeBinary=save)
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
target.fromCSV(file,'K2-232b',skiprows=2)

## Instantiate the Simulator with 330 draws, precision of 1 m/s, and 3 CPUs (multiprocessing)
run = rmcmc.fakit.Simulator(target,ndraws=300,precision=1,nproc=3)
## Simulate some data with an exposure time of 180 s
run.dataSimulator(exposure=180)

## Save results or not
save = 0
## 
run.MC([0.0,10],writeBinary=save)
if save:
	results = rmcmc.fakit.readBinary('K2-232b.pkl')
	rmcmc.makit.plotKDE(results)
else:
	rmcmc.makit.plotKDE(run.results)


```