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

# Lineshape fitting program, rec-by-rec version

This is a universal service for the 
treatment of the recordings of the single spectral lines' profiles. Leaving beyond 
the scope some details of the spectrometers and measurements techniques used by 
the researchers, as well as primary processing the raw data acquired, it is focused 
on the analysis of the spectrometer signal-versus-frequency dependencies, which are quite 
common in a large variety of thermodynamic conditions and frequency intervals.

This is a version for the recording-by-recording processing: each recording is
processed separated from other ones, and the result of the processing is line 
shape parameters (central frequency, half-width, etc) corresponding to each 
single recording in a batch. This results are further used for the analysis of the 
pressure- and temperature-dependencies of the line shape parameters and getting 
normalized to pressure coefficients in case of batch processing of the recordings.

The reason for creating this service was following: a lot of spectroscopy researchers
face quite common routine, the analysis of the recording batches of some spectral line (or
several lines). These batches usually correspond to the evolution of the spectral shape under
conditions of step by step changing pressure of the absorbing of foreign gas. This software 
is an attempt to generalize these routines and work out some time-saving approach to the 
analysis of the recorded spectra.

# Interface and working with the application

The interface consists of several blocks having certain functionality.

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