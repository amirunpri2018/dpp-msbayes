Table of Contents
=================

 -  [Overview](#overview)
 -  [Requirements](#requirements)
 -  [Installation](#installation)
 -  [Citing dpp-msbayes](#citing-dpp-msbayes)
 -  [Documentation](#documentation)
 -  [Literature Cited](#literature-cited)
 -  [License](#license)

Overview
========

This is a modified version of the `msBayes` package. It is modified by [Jamie
Oaks](http://www.phyletica.com) from version 20130611 (revision 146) available
via SVN here: <http://sourceforge.net/projects/msbayes/>.

`dpp-msbayes` implements approximate Bayesian model choice to infer comparative
phylogeographic histories.

In this distribution of the code, the phylogeographic model has been
re-parameterized, and different probability distributions are made available to
express prior uncertainty about many of the models' parameters.

I strongly recommend you use this software via the Python multiprocessing
package `PyMsBayes` (<https://github.com/joaks1/PyMsBayes>). This package comes
bundled with pre-compiled `dpp-msbayes` executables (and the original `msBayes`
executable, as well, which you can use interchangeably). Also, the test suite
of `PyMsBayes` is *much* more extensive. If you do opt to use `PyMsBayes`,
however, I still encourage you to read the `Documentation` section below.

Requirements
============

The primary requirement is the GNU Scientific Library (GSL;
<http://www.gnu.org/software/gsl/>). To use the bundled `R` code for
local-linear regression/multinomial logistic regression adjustment of posterior
estimates, you will also need `R` and add-on packages `locfit`, `VGAM`, and
`KernSmooth`. However, there are alternative posterior adjustment procedures
availabe (e.g., GLM regression implemented by `ABCtoolbox`;
<http://www.cmpg.iee.unibe.ch/content/softwares__services/computer_programs/abctoolbox/index_eng.html>).

Installation
============

This version of `msBayes` now follows the GNU Build System via Autotools.  For
specifics, please see `INSTALL.md`.

Citing dpp-msbayes
==================

If you publish results obtained using this software, please cite:

 > Oaks, J. R. (2014). An Improved Approximate-Bayesian Model-choice Method for
 > Estimating Shared Evolutionary History. arXiv:1402.6303 [q-bio:PE].
 > <http://arxiv.org/abs/1402.6303>

Also, please cite the work describing the version of `msBayes` this code is
based on:

 > Huang, W., N. Takebayashi, Y. Qi, and M. J. Hickerson, 2011.
 > MTML-msBayes: Approximate Bayesian comparative phylogeographic
 > inference from multiple taxa and multiple loci with rate
 > heterogeneity. BMC Bioinformatics 12:1.

Also, more details of the model implemented in the original version
of `msBayes` can be found in:

 > Oaks, J. R., Sukumaran, J., Esselstyn, J. A., Linkem, C. W., Siler, C. D.,
 > Holder, M. T., & Brown, R. M. (2012). Evidence for climate-driven
 > diversification? A caution for interpreting ABC inferences of simultaneous
 > historical events. Evolution, 67(4), 991-1010.
 > doi:10.1111/j.1558-5646.2012.01840.x.

Documentation
=============

Documentation for the original version of `msBayes` from which this software
has been modified can be found in `documents/mtml-msBayes-manual.pdf` and at
<https://docs.google.com/document/d/15heQlz60cGe6GWKcXqf1AIYMBZ6p2sKuGct-EoEHiNU/edit>.
These sources are still helpful for using this version of the software, but
please read below to understand the differences in the model and how it is
used.

This section provides some information needed to use this software. It is
probably best to use this method via the Python multiprocessing package
`PyMsBayes` (<https://github.com/joaks1/PyMsBayes>), and thus this section
does not provide all the details for doing a full analysis from start to
finish "by hand." All that is needed to use this method via `PyMsBayes` is
an understanding of the `msBayes` configuration file. Below is an example
configuration file:

    concentrationShape = 0
    concentrationScale = 0
    thetaShape = 1
    thetaScale = 0.001
    ancestralThetaShape = 1
    ancestralThetaScale = 0.001
    thetaParameters = 001
    tauShape = 1.0
    tauScale = 2.0
    bottleProportionShapeA = 5
    bottleProportionShapeB = 1
    bottleProportionShared = 1
    numTauClasses = 0
    migrationShape = 1
    migrationScale = 1
    recombinationShape = 0
    recombinationScale = 0
    constrain = 0
    subParamConstrain = 111111111
    
    BEGIN SAMPLE_TBL
    pair0	locus0	1	1	10	10	15.6	1000	0.35	0.28	0.13	pair0.locus0.fasta
    pair1	locus0	1	1	10	10	15.6	1000	0.35	0.28	0.13	pair2.locus0.fasta
    pair2	locus0	1	1	10	10	15.6	1000	0.35	0.28	0.13	pair3.locus0.fasta
    pair3	locus0	1	1	10	10	15.6	1000	0.35	0.28	0.13	pair3.locus0.fasta
    END SAMPLE_TBL

A configuration file has two parts, (1) the *sample table* delimited by the
`BEGIN SAMPLE_TBL` and `END SAMPLE_TBL` lines, and (2) the *preamble*
containing keyword arguments prior to the sample table.

There are two differences between this modified implementation of `msBayes` and
the original when it comes to the preamble of configuration files.

 1. The preamble is processed differently in this version. The original
    implementation is case-sensitive regarding the keywords in the preamble.
    I.e., the camel-case type for `subParamConstrain` is strictly enforced.
    More importantly, the old version quietly uses default settings when any
    options are mis-typed.  in other words if you specify `uppertheta = 0.01`
    in a config (note the lower case "t"), it will quietly use the default
    setting for `upperTheta` and *not* report any warning or error.

    In this version of the software, the keywords in the preamble are case
    insensitive, and if any unrecognized keywords are encountered, an error is
    reported and the process exits (i.e., crashes).

 2. In this version, the options in the preamble are different and allow you to
    parameterize the model and specify prior probability distributions on
    parameters, as follows:
     -  `concentrationShape`/`concentrationScale`
        
        Rather than assume the U-shaped prior of the original `msBayes` (see Oaks
        et al., 2012), this implementation allows the user to specify a
        Dirichlet-process prior over all possible divergence models.  Specifically
        when both of these options are positive numbers, they define the shape and
        scale parameters of a gamma hyperprior controlling the concentration
        parameter of a Dirichlet-process prior on divergence models.
    
        If either or both are zero (as above), this specifies that a uniform prior
        over order-independent divergence models (*not* over the number of
        divergence events; I.e., Psi) is to be used.
    
        If either or both are -1.0 or less, then the original U-shaped prior is
        used. It's an odd combination of a discrete uniform distribution on the
        number of divergence events (Psi) and a multinomial distribution on the
        divergence model given the number of events (with the constraint that there
        must be Psi elements in the model; see Oaks et al. 2012 for full details).
     -  `thetaParameters`
        
        This implementation allows the user full control over the parameterization
        of the population size parameters for each population pair, which we will
        refer to theta_a, theta_d1, and theta_d2 for the ancestral, and two
        descendant populations.
    
        `000` is one extreme, where theta_a, theta_d1, and theta_d2 are all
        constrained to be equal for each population pair (population sizes will
        still vary among pairs).
    
        `012` is the other extreme, where theta_a, theta_d1, and theta_d2 are all
        estimated as independent parameters. This is most similar to the original
        `msBayes`, however, the descendant population sizes are constrained to be
        negatively correlated in the original implementation (see Oaks 2012).
    
        Another example is `001`: the descendant populations share the same size
        parameter, but the ancestral population size is free to vary.
    
        For `011` and `010`, one of the descendant population is constrained to the
        same size as the ancestral, and the other is free to vary.
    
        Note, this indicator sequence should always be three integers long, always
        start with `0` and increment by 1 whenever a free parameter is added.
     -  `thetaShape`/`thetaScale`
    
        These settings define the shape and scale parameters of a gamma prior
        distribution on population size parameters (scaled by the per-site mutation
        rate (u); Nu).
     -  `ancestralThetaShape`/`ancestralThetaScale`
        
        If these settings are both provided and both are positive, they define the
        shape and scale parameters of a gamma prior on the sizes of ancestral
        populations.
        
        If they are excluded, or both are zero, the `thetaShape` and `thetaScale`
        settings are used for the gamma prior on ancestral population size
        parameters.
     -  `tauShape`/`tauScale`
        
        These settings define the shape and scale parameters of a gamma prior
        distribution on divergence time parameters. The units are in coalescent
        units, Nc generations, where "Nc" is a constant reference population size
        based on the mean of the theta prior (defined by `thetaShape` and
        `thetaScale`). If we use theta_prior_mean to represent the mean of the
        theta prior, then Nc is theta_prior_mean / u, where "u" is the per-site
        mutation rate. Thus, you can convert these "Nc generations" units to the
        number of generations by assuming a mutation rate and multiplying by
        (theta_mean_prior/u). See Oaks (2012) for more details.
     -  `bottleProportionShapeA`/`bottleProportionShapeB`
        
        If both are positive, settings define the shape parameters alpha and beta,
        respectively, of a beta prior distribution on the magnitude parameters (in
        units of the proportion of the population size) of a post-divergence
        bottleneck in each of the descendant populations.
    
        The bottleneck proportions are in terms of the proportion of the effective
        population size that remains following the bottleneck. Thus a value of 0.95
        would mean that bottleneck reduces the effective population size by 5%.
    
        If either or both are zero or less, there is no post-divergence population
        bottleneck (i.e., these parameters along with the timing of the bottleneck
        are removed from the model).
    
        NOTE, there are also a parameters in the model for the timing of the end of
        the bottleneck (it begins at speciation in forward time). There is one of
        these parameters for each descendant population (i.e., the descendant
        populations of each pair share the same bottleneck-end-time parameter).
        Thus if either or both of the
        `bottleProportionShapeA`/`bottleProportionShapeB` settings are zero or
        less, you are also removing these bottleneck timing parameters from the
        model. This means you are removing 3*Y parameters from the model, where "Y"
        is the number of taxon pairs.
     -  `bottleProportionShared`
    
        If this option is zero, then there are two free bottleneck-magnitude
        parameters for each population pair (one for each descendant population).
        If it is any other number, then there is one bottleneck-magnitude parameter
        for each population pair (i.e., the descendant populations of each pair
        share the same bottleneck parameter).
    
        NOTE, this setting is overridden if either or both of the
        `bottleProportionShapeA` or `bottleProportionShapeB` settings is zero or
        less (because then there is no bottleneck parameters at all).
     -  `migrationShape`/`migrationScale`
    
        These settings define the shape and scale parameters of a gamma prior
        distribution on symmetric migration parameters (in units of the number of
        gene copies per generation).
    
        If either or both are zero or less, there is no migration in the model.
     -  `recombinationShape`/`recombinationScale`
    
        These settings define the shape and scale parameters of a gamma prior
        distribution on the intragenic recombination rate parameters.
    
        If either or both are zero or less, there is no recombination in the model.
     -  `numTauClasses`
        
        If this setting is zero, the number of divergence events is free to vary
        according to the prior on divergence models specified by
        `concentrationShape` and `concentrationScale`.
    
        If greater than zero, then the model is constrained to `numTauClasses`
        divergence events.
     -  `constrain`/`subParamConstrain`
        
        I strongly recommend *not* changing these settings. The software is
        completely untested in how it behaves when constrained models are specified
        with these options.

Literature Cited
================

 > Oaks, J. R., Sukumaran, J., Esselstyn, J. A., Linkem, C. W., Siler, C. D.,
 > Holder, M. T., & Brown, R. M. (2012). Evidence for climate-driven
 > diversification? A caution for interpreting ABC inferences of simultaneous
 > historical events. Evolution, 67(4), 991-1010.
 > doi:10.1111/j.1558-5646.2012.01840.x.

Acknowledgements
================

This software greatly benefited from funding provided to [Jamie
Oaks](http://www.phyletica.com) from the National Science Foundation (DEB
1011423 and DBI 1308885), University of Kansas (KU) Office of Graduate Studies,
Society of Systematic Biologists, Sigma Xi Scientific Research Society, KU
Ecology and Evolutionary Biology Department, and the KU Biodiversity Institute.

License
=======

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program. If not, see <http://www.gnu.org/licenses/>.

See "LICENSE" for full terms and conditions of usage.

