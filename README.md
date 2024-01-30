---
title: Shiny Single Fit
licence: MIT
colorFrom: green
colorTo: gray
sdk: docker
app_file: app.py
pinned: false
emoji: 📈
---

# Lineshape fitting program

This is a universal service for the 
treatment of the recordings of the single spectral lines' profiles. Leaving beyond 
the scope some details of the spectrometers and measurements techniques used by 
the researchers, as well as primary processing of the raw data acquired, it is focused 
on the analysis of the spectrometer signal-versus-frequency dependencies, which are quite 
common in a large variety of thermodynamic conditions and frequency intervals.

This version includes two approaches to the data treatment. The **first** (and the most
common one) is **recording-by-recording** processing: each recording (i.e. signal-vs-frequency obtained in some certain frequency interval and at some certain
conditions such as absorber pressure, perturber pressure and temperature) is
processed separated from other ones, and the result of the processing is line 
shape parameters (central frequency, half-width, etc) corresponding to each 
single recording from a batch and the conditions it was acquired under. These results 
are further used for the analysis of the pressure- and temperature-dependencies 
of the line shape parameters and getting normalized to pressure coefficients (broadening 
coefficient, shifting coefficient etc.) in case of batch processing of the recordings.

The **second** approach is so-called **multispectrum fitting**, or **multifit**. Here all recordings in a batch are analyzed together: parameters of the line shape related to some certain spectral line of a certain molecule have well-known regular analytical dependencies on the molecular densities (i.e. pressures) of the gases in the studied mixture and on the temperature. When we analyze batch of the recordings acquired at various pressures (and sometimes various temperatures) we can use a priory knoweledge of these dependecies and therefore derive normalized collisional parameters (broadening 
coefficient, shifting coefficient etc.) directly, bypassing the stage of obtaining lihe shape parameters for each recording isolated.

The reason for creating this service was following: a lot of spectroscopy researchers
face quite common routine, the analysis of the recording batches of some spectral line (or
several lines). These batches usually correspond to the evolution of the spectral shape under
conditions of step by step changing pressure of the absorbing of foreign gas. This software 
is an attempt to generalize these routines and work out some time-saving approach to the 
analysis of the recorded spectra for the most common situations.

# Interface and working with the application

The interface consists of several blocks having certain functionality.

## Data upload

The first block is for uploading your spectra and some additional data we need for the analysis. This block is common for line-by-line and multispectrum fitting routines.

### Loading spectra

Press "Browse" button under "Load your spectra" label and upload any amount of the files
containing experimental data. Files should contain two columns: the first one is frequency 
(currently, on 2023.12.26, only MHz units are supported), the second one is spectrometer signal
or absorption coefficient in any arbitrary units. Columns starting from number three, 
if existing, will be ignored. In the "Commented lines" input value you can specify 
how the lines in your files are commented (the most common options are `#` or `//`, the latter
is set as default).

### Loading recordings metadata
 
Press "Browse" button under "Load recordings info" and upload a single file with the list of filenames
and experimental conditions the spectra were recorded under. It should contain 7 columns:

1. filename
2. pressure of the absorbing gas (Torr)
3. foreign gas pressure (Torr)
4. temperature (in C or K, recalculated to K inside the service automatically)
5. frequency deviation (see the details below)
6. recording type (see the detail below)
7. cell length

Column 6 currently supports the following values:

- 0 for the direct observation of the line profile (like in cavity spectrometer)
- 1 for the FM mode recordings from the RAD spectrometer
- 2 for the FM mode recordings from the video spectrometer

FM mode is approach when the frequency of the radiation source oscillates between two close
values, which leads to the corresponding oscillations of the spectrometer signal. 
The amplitude of these oscillations is proportional to the absorption difference between two
mentioned frequency values. The difference between latter two options is due to the 
fact that RAD spectrometer (radio-acoustic detector) produces signal proportional to the
power absorbed by the molecules, while video spectrometer produces signal proportional
to the leftover power passed through the gas.

For recording type (col 6) value of 0 the deviation (col 5) should be set to 0 (this value 
is not used anyway). For col 6 values of 1 and 2, the col 5 value should be non-zero and set
to the frequency deviation value used in FM mode for the corresponding recording (in the same
units that frequency is measured in the recorded spectra).

Here is some example of metadata file lines:

```text
#file           Pself   Pforeign T     dev typ  L 
filename_1.dat	24.768	725.992	297.17	0	0	1.
filename_2.dat	0.095	0.523	296.58	1.5	1	10.
filename_3.dat	0.010	0.034	298.1	0.2	2	100.
```

The first line may contain column headers (just not to confuse the columns), 
in this case you need to switch on the control "Column names in the 1st line"
near the file download button. Filename_1 corresponds to the directly observed line profile
from the resonator (cavity) spectrometer, dev=0 and L=1 (as they are ignored anyway).
Filename_2 corresponds to the FM mode recording from RAD spectrometer (or any similar type
producing signal proportional to the absorbed radiation power). Filename_3 corresponds to
the FM mode recording from videospectrometer producing signal proportional to the 
passthrough power. For both latter cases frequency deviation value is specified (to calculate
proper recorded signal profile shape from the spectral line shape), as well as cell length
value L (to consider Bouger's law, which is essential when observing strong lines, where
absorption multiplied to the pass length is big enough).

### Loading partition function information

In all cases this service is designed for, we need to calculate spectral line integrated
intensity properly in according to the thermodynamical conditions of the experiment, 
namely absorber partial pressure and gas temperature. For this, we should take into 
account partition function which impacts the temperature dependence of the integrated 
intensity of the spectral line in our recordings. Pressing the "Browse" button under
the label "Load partition function" you can upload a partition function file of 
two columns (temperature (K) and partition function value). Such data can be 
downloaded from the [HITRAN database](https://hitran.org). In case where partition function
is absent or unnecessary (e.g. temperature is close to the reference value of 296 K
and we do not need much precision), you can just switch the control "Use no partition 
function" to use value of 1. for processing of your data.

## Line-by-line processing section

## Parameters selection and guess values

First block is to define, what parameters of the line shape we are going to adjust,
and which one are fixed to some constant values. On the left side of this block there is
a checklist with line shape parameters enlisted, checked parameters will be the adjusted ones.

On the right side there are placeholders for input values, which are necessary to model the
spectrum. Note that 
- intensity value is given normalized to the absorber concentration in HITRAN units 
- collisional parameters (such as broadening, shifting, etc) are also given as coefficients,
 i.e. are normalized to pressure

These values are used to calculate guess values for the parameters of the model profile
(like integral intensity, HWHM, mixing parameter, etc) which fitting procedure starts with, 
the guess values are calculated in accordance with the absorbing and foreign gases pressure
and temperature.

## Preview your data

When experimental data, metadata and partition function are loaded, and line parameters for initial guess values are specified, you can preview your data to check it, and also see the model profiles with the guess starting values you used in the first block. If some necessary information is absent, the message will appear saying, what you should add to continue data processing.

## Fitting model profile to the recordings

When you are happy with your uploaded data and preview, press the button "Fit model to 
recordings" to start fitting procedure. When it is over, the residuals will be shown below.

## Results: viewing the plots and downloading the treatment data

After fitting the model function to your experimental data, you have a set of 
values of various profile parameters (HWHM, speed-dependent part of the collisional width,
central frequency, mixing parameter, integrated intensity) corresponding to each 
individual recording. There is an option to preview dependencies of these parameters
versus pressure. Below the residuals, you can find switches: 

- "Foreign pressure": OFF if pure absorbing gas spectrum was recorded, 
ON for the spectrum of mixture of the absorbing and foreign gas  
- "Norm to self P": OFF if partial pressure of the absorbing gas is the same in all
recordings, ON if partial pressure of the absorbing gas is different in different
recordings

"Unshifted center" data input is used for better view of the central frequency versus 
pressure. The specified value is subtracted out of the fitted central frequency values
to make a plot.

"Plot vs pressure" button shows the figures for the pressure dependencies.

"Download results" downloads single text file with the fit results (both parameter value
and fit uncertainty).

"Download residuals" downloads .zip-file with the fit residuals of all uploaded experimental
recordings packed.