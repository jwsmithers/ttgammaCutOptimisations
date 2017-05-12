# ttgamma Cut Optimisation code

## Usage

`$> git clone https://github.com/jwsmithers/ttgammaCutOptimisations.git ; cd ttgammaCutOptimisations` 

`$> ./cutOptimisation.py`

If cutflow is required after all cuts have been optimized:

`$> ./offline_cutflow.py`

## How to optimise variables
### Final cutflows
Open `cutOptimizations.py` and scroll right to the bottom. You will see commented out functions that call main(). For the most basic test, 
an `allCut` jobs has been uncommented. This means the optimisation will run once (for ejets and mujets) and assumes your cuts have been chosen
(look in lines ~255 -> ~266). There are two way to run this:

`_ = Parallel(n_jobs=-1, verbose=5, backend="multiprocessing") \
( delayed(main)(cutName="full_cuts",region=i,cut_range=[1],allCuts = True) for i in regs)`

will run ejets and mujets in parallel. But the output is a little difficult to read. You can also just uncomment 

`main(cutName="full_cuts",region="ejets",cut_range=[1], allCuts = True)`

to do it sequentially.

### Individual variables
Ok, but what if you want to optimise a single variable? To do this, comment out all the main() functions at the bottom of `cutOptimisations.py`
and uncomment out a variable and range of cuts you would like to run over. I.e.

`main(cutName="ph_HFT_MVA",region="ejets",cut_range=[0,0.05,0.10,0.15,0.20,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95])`

This will cut on each value specified in `[]` for the `ph_HFT_MVA` variable in the `ejets` region. *But before you do ./cutOptimisations.py*,
you have to set the `cutValue` in the loop. 

Go to line 252 and change `ph_HFT_MVA_cutValue = 0`, to `ph_HFT_MVA_cutValue = cutValue`. Now run the code. Printouts are useful, but also look in
outputs_ directories to gets some root files and already plotted pngs.

## Systematics
Systematics are quite central to this optimization. These are defined in lines ~640 -> ~685. You can always add more, but remember to account for
this in `optimizationFunctions.py`, `_sigmaSys` function.
