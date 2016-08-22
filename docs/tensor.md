# Tensor decompositions of trial-structured calcium imaging datasets

<i><b>Note:</b> The tensor analyses are under active development. I'll try to update the documentation here continually as I update the code.</i>

## Background material

A *tensor* is a higher-order numeric array. The *order* of a tensor is the number of indices/dimensions it holds. In MATLAB, `randn(10,10,10)` creates a tensor of order 3.

There are many tensor decompositions one can try, but the simplest and easiest to interpret is the canonical polyadic decomposition (CPD). We will focus on fitting this decomposition on an order-3 tensor `neurons x intra-trial phase x trials`, which gives us three sets of low-dimensional factors/components. The first set of components provides a low-dimensional representation of the neurons, the second set of components provide a representation of the within-trial dynamics of the neurons, and the final set of components provide a representation of how the within-trial dynamics change across trials and days.

[Kolda & Bader (2008)](#) is a good review that covers the more technical details of tensor decomposition.

## Basic Walkthrough

Start by loading data and turning it into a `MultiDay` object named `md`. Refer to `docs/ds_quickstart.md` and `docs/md_quickstart.md`. Once you have the data loaded (make sure it is called `md`) you should be able to generate all the relevant plots by running:

```matlab
>> tensor_demo
```

The decomposition can take moderately long to fit when taking data from multiple days (particularly the scree plot, explained below) -- it is probably a good idea to start with just one day.
The sections below walks through the demo script and provides brief explanations of each of the generated plots:

#### Converting a `MultiDay` object to a data tensor

The first step is to export data from the `MultiDay` object into a `K x 1` cell array. Each entry in the cell array is a matrix with `N` rows and a variable number of columns (depending on the number of frames captured on each trial). The neurons are matched across all trials (e.g. the first row in all matrices is the trace corresponding to the same neuron). The following command produces the cell array (`X`) as well as two matrices for mapping the neurons and trials back to their position in the `MultiDay` object.

```matlab
[X, neuron_map, trial_map] = export_multiday_trace(md)`
```

***Filtering trials:*** the `export_multiday_trace` function also has the ability to filter out trials based on their attributes. For example, you may want to construct a tensor containing all trials where the mouse started in the east arm and ended in the north arm:

```matlab
export_multiday_trace(md, 'start', 'east', 'end', 'north')
```

You can also filter out correct/incorrect trials:

```matlab
export_multiday_trace(md, 'start', 'correct', 1) % correct trials
export_multiday_trace(md, 'start', 'correct', 0) % incorrect trials
```
#### Warping the trials to a common length

Next we need to convert `X` from the cell array format to a 3-dimensional array. This is achieved by the poorly named timewarp function, `X = timewarp(X)`. At the moment, we just do linear interpolation/stretching to convert all trials to the same length. In the future I want this function to support other options for warping the data:

- [ ] Dynamic Barycentric Averaging (DBA), which uses dynamic time warping to align all trials.
- [ ] Mapping the trial phase to the position of the mouse

These items should not be too difficult and are up for grabs!

#### Warping the trials to a common length

Next we need to convert `X` from the cell array format to a 3-dimensional array.
This is achieved by the poorly named timewarp function, `X = timewarp(X)`.
At the moment, we just do linear interpolation/stretching to convert all trials to the same length.
In the future I want this function to support other options for warping the data:

- [ ] Dynamic Barycentric Averaging (DBA), which uses dynamic time warping to align all trials.
- [ ] Mapping the trial phase to the position of the mouse

These items should not be too difficult if anyone would like to contribute them!

#### Fitting the CPD model and making a scree plot

Tensor decompositions are non-convex optimization problems and are NP-hard in terms of worst case analysis.
Thus, it is prudent to fit the model with multiple random initializations.
We would like to get a sense of how well fit the model is to the data as a function of model complexity (i.e. we need to determine how many factors to fit to the data).
To do this, we use the `fit_cpd` function:

```matlab
[cpd_list, rsq] = fit_cpd(X) % by default, fit models with rank 1-15 from 10 random starts each
```

The results are returned in a 1-dimensional struct array `cpd_list` which holds the outcome of each optimization.
The vector `rsq` holds the coefficient of determination (R<sup>2</sup>) for each optimzation as a measure of the fit.
These results are commonly summarized with a [scree plot](#).
The following command will produce a nicely formatted scree plot, given the struct array of cpd fits:

```matlab
scree_cpd(cpd_list);
```

![CPD scree plot](cpd_scree.png)

#### Visualizing the factors

Typically, we will just pick the model with the highest R<sup>2</sup> to analyze further:

```matlab
[~,best_idx] = max(rsq);
cpd = cpd_list[best_idx];
```

Next, let's visualize the factors. Each factor is a triplet of three vectors, and the number of factors is the *rank* of the model.
Below, I visualized the factors for a rank 15 model using the command:

```matlab
cpd_factor_plots(cpd,md,trial_map)
```

![CPD factors](cpd_factors.png)

The left column of plots shows the 15 *neuron factors*.
Because the ordering of the neurons in the data tensor `X` isn't especially meaningful, I've sorted the neurons from highest to lowest on the first factor (i.e the top plot).
The second column of plots shows the 15 *across trial factors*, with the east and west starts colored in blue and red respectively.
The third colum of plots shows the 15 *within trial factors*.

You can change the coloring of the across trial factors (middle column) with the `trialcolor` parameter:

```matlab
cpd_factor_plots(cpd,md,trial_map,'trialcolor','correct')
```

Which will produce a plot with the error trials colored red and correct trials colored blue.

#### Visualizing the model fit

It is useful to view the model's prediction and the raw data on the same plot, this can be done with the following commands:

```matlab
% get the full reconstructed tensor from the model
Xest = full(cpd.decomp);
Xest = Xest.data;

% plot fit across neurons
visualize_fit(X,Xest,1,md,trial_map);
```

The last command should produce a series of plots that look like this (press the space button to go to the next plot and `Control-C` / `Command-C` to interrupt the program).

![CPD factors](cpd_fit1.png)

This figure plots a single neuron `X(i,:,:)` on 25 random trials.
The raw data is plotted in red/blue traces respectively denoting east/west trial starts.
The model fit is the blace trace in all plots.

It is also possible to examine the model fit across all trials, i.e. examining slices through the third mode of the tensor `X(:,:,i)`.

```matlab
% plot fit across trials
visualize_fit(X,Xest,3,md,trial_map);
```

![CPD factors](cpd_fit3.png)

#### Visualizing the residuals and outlier detection

Another question of interest is whether the model is better fit to some trials or to some neurons more than others.
To look at this we examine the residuals for each neuron, trial phase, and trial (`X(n,t,k) - Xest(n,t,k)`) as well as the squared error `(X(n,t,k) - Xest(n,t,k))^2`.

```matlab
% plot fit across trials
visualize_resids(X,Xest,md,trial_map);
```

Produces a plot like:

![CPD residual analysis](cpd_resids.png)

## Things to improve:

- [ ] Measuring distance in trial space as a metric of learning / strategy-shifting
- [ ] More detailed residual analysis and better documentation
