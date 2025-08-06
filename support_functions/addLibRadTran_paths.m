% --- Add Folders to Path specific to this computer -----


if strcmp(whatComputer,'anbu8374')==true
    % ----- LASP MAC ------
    % LibRadTran Data Folders
    addpath('/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/solar_flux/');
    addpath('/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/wc/');
    addpath('/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/albedo/');
    addpath('/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/altitude/');
    addpath('/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/atmmod/');
    addpath('/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/data/correlated_k/');

    % LibRadTran Bin folder
    addpath('/Users/anbu8374/Documents/LibRadTran/libRadtran-2.0.4/bin/');

elseif strcmp(whatComputer,'andrewbuggee')==true
    % ----- MACBOOK ------


    addpath('/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/LibRadTran/libRadtran-2.0.4/data/solar_flux/');
    addpath('/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/LibRadTran/libRadtran-2.0.4/data/wc/');
    addpath('/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/LibRadTran/libRadtran-2.0.4/data/albedo/');
    addpath('/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/LibRadTran/libRadtran-2.0.4/data/altitude/');
    addpath('/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/LibRadTran/libRadtran-2.0.4/data/atmmod/');
    addpath('/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/LibRadTran/libRadtran-2.0.4/data/correlated_k/');

    % LibRadTran Bin folder
    addpath('/Users/andrewbuggee/Documents/CU-Boulder-ATOC/Hyperspectral-Cloud-Droplet-Retrieval/LibRadTran/libRadtran-2.0.4/bin/');


elseif strcmp(whatComputer,'curc')==true
    % ----- CU SUPERCOMPUTER ------

    % Add paths to all functions needed
    addpath('/projects/anbu8374/Matlab-Research/')
    addpath('/projects/anbu8374/software/gsl-2.6/')
    addpath('/projects/anbu8374/software/libRadtran-2.0.5/')


end

