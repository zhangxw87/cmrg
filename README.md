# Multivariate Regression with Grossly Corrupted Observations

This is a collection of MATLAB codes that implement our algorithm and comparison methods in paper.

## Folders
+ `codes/`: contains all the source codes implementing the above 11 algorithms
+ `Datasets/`: default path to local files storing datasets
+ `results/`: stores the computed results. This folder also contains a few scripts which can be used to plot figures in paper

## Usage
+ Run _Example1.m_, _Example2.m_, _Example3.m_, _Example4.m_, _Example5.m_ to get results on simulated data. For each example, there is a corresponding parallel version (e.g. _Example1_parallel.m_) which repeats the experiments 100 times. 
+ Run _Syn_hand0.m_, _Syn_hand1.m_, _Syn_hand2.m_ to get results for simulated hand pose estimation experiment. Each _.m_ file corresponds to a different gross error rate. 