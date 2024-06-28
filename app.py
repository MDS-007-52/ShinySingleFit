from shiny.types import FileInfo, NavSetArg
from shiny import *  # App, Inputs, Outputs, Session, render, ui
# from shinywidgets import *
# import ipywidgets as widgets
# import ipysheet
import numpy as np
from pathlib import Path
import zipfile
import io
import math
from matplotlib import pyplot as plt
import pandas as pd
import time

# from SDRP import SDRP
from fitmodels import *
from subr_fit import fit_params, fit_uncertainties
from constants import *
from qpart import *

app_ui = ui.page_fluid(
    # ui.head_content(
    #     ui.tags.script(
    #         src="https://polyfill.io/v3/polyfill.min.js?features=es6"
    #     ),
    #     ui.tags.script(
    #         "if (window.MathJax) MathJax.Hub.Queue(['Typeset', MathJax.Hub]);"
    #     ),
    # ),
    ui.panel_title(ui.h1("Spectral line shape fit demo version")),
    ui.row(ui.h2("Load your data here")),

    # this layout is for uploading "markup" for recordings, i.e. pressures, temperature, recording type etc
    ui.layout_sidebar(
                      ui.panel_sidebar(ui.row(
                                              ui.column(6, ui.input_file("aux_data", "Load recordings info",
                                              accept=["*.*"], multiple=False)),
                                              ui.column(6, ui.input_switch('s_header', "Column names in the 1st line"))
                                              ),
                      ui.row(
                             ui.column(6, ui.input_file("partition", "Load partition function",
                                                        accept=["*.*"], multiple=False)),
                             ui.column(6, ui.input_switch('s_nopart', "Use no partition function"))
                                )
                         ),
                      ui.panel_main(ui.output_ui("aux_data_show"))
                      ),

    # an interface for loading recordings
    ui.row(
        ui.column(4, ui.input_file("input_files", "Load your spectra", accept=["*.*"], multiple=True)),
        ui.column(4, ui.input_text("text_comment", "Commented lines", value="//"))),

    # Tab interface with two tabs for LBL and multifit treatment respectively
                    
    ui.navset_tab(
                  # tab for LBL
                  ui.nav_panel("Line-by-line treatment", 
                               # this sidebar layout is for marking which parameters to adjust and initial parameters values
                               ui.layout_sidebar(
                                                ui.panel_sidebar(ui.h4("Mark adjustable"),
                                                                ui.input_checkbox_group("jac_check", "",
                                                                                        single_params_dict, 
                                                                                        selected=single_params_adjustable_init),
                                                                #ui.input_action_button("b_check", "Check"),
                                                                #ui.output_text_verbatim("jac_flag_out", placeholder=True)
                                                                ),
                                                ui.panel_main(ui.h4("Parameters and coefficients for initial values"),
                                                            ui.row(
                                                                ui.column(4, ui.input_numeric("I0", "Intensity, 1e-25 cm/mol", value=33.)),
                                                                ui.column(4, ui.input_numeric("elow", "Lower level energy, 1/cm", value=0.)),
                                                                ui.column(4, ui.input_numeric("molm", "Molecular mass, a.m.u.", value=28.))
                                                            ),
                                                            ui.row(
                                                                ui.column(3, ui.input_numeric("f0", "Line center, MHz", value=115271.202)),
                                                                ui.column(3, ui.input_numeric("ngam", "T-dep power", value=0.75)
                                                                            ),
                                                                ui.row(
                                                                    ui.column(3, ui.input_numeric("d0", "d0 self, MHz/Torr", value=0.000)),
                                                                    ui.column(3, ui.input_numeric("d0f", "foreign, MHz/Torr", value=0.000))
                                                                ),
                                                                ui.row(
                                                                    ui.column(3, ui.input_numeric("d2", "d2 self, MHz/Torr", value=0.00)),
                                                                    ui.column(3, ui.input_numeric("d2f", "foreign, MHz/Torr", value=0.00))
                                                                ),
                                                                ui.row(
                                                                    ui.column(3, ui.input_numeric("g0", "g0 self, MHz/Torr", value=3.375)),
                                                                    ui.column(3, ui.input_numeric("g0f", "foreign, MHz/Torr", value=2.750))
                                                                ),
                                                                ui.row(
                                                                    ui.column(3, ui.input_numeric("g2", "g2 self, MHz/Torr", value=0.33)),
                                                                    ui.column(3, ui.input_numeric("g2f", "foreign, MHz/Torr", value=0.28))
                                                                ),
                                                                
                                                                ui.row(
                                                                    ui.column(3, ui.input_numeric("y0", "y self, 1/Torr", value=7.e-6)),
                                                                    ui.column(3, ui.input_numeric("y0f", "foreign, MHz/Torr", value=7.e-6))
                                                                ),
                                                                ui.row(
                                                                    ui.column(3, ui.input_numeric("nu", "N_vc self, MHz/Torr", value=0.15)),
                                                                    ui.column(3, ui.input_numeric("nuf", "N_vc foreign, MHz/Torr", value=0.15))
                                                                ),
                                                            )
                                                            )
                                                ),
                                    ui.row(ui.input_action_button("b_preview", "Preview line-by-line")),
    
                                    ui.output_plot("preview_spectrum"),
                                    ui.output_text("preview_status"),

                                    # an interface for fitting the model to the spectra and visualisation of the result

                                    ui.input_action_button("b_fit", "Fit model to recordings"),
                                    ui.output_text("fit_status"),
                                    ui.output_plot("preview_fit"),
                                    ui.output_ui("fit_result_table"),
                                    ui.row(ui.column(3, ui.input_switch('s_foreign', "Foreign pressure")),
                                        ui.column(3, ui.input_switch('s_normself', "Norm. to self P")),
                                        ui.column(4, ui.input_numeric('f0_shift', "Unshifted center", value=115271.202))
                                        ),
                                    ui.row(
                                        ui.column(3, ui.input_action_button("b_coefs", "Plot vs pressure")),
                                        ui.column(3, ui.download_button("download_params", "Download results")),
                                        ui.column(3, ui.download_button("download_resid", "Download residuals")),
                                        ui.column(3, ui.input_text("resid_prefix", "Residual file prefix", value="r_"))
                                        ),

                                    ui.output_plot("preview_params", height="1200px")
                               ),
                   # tab for multifit
                   ui.nav_panel("Multifit treatment", 
                                ui.layout_sidebar(ui.panel_sidebar(ui.h4("Mark adjustable"),
                                       ui.input_checkbox_group("jac_check_multi", "",
                                                               multi_params_dict, 
                                                               selected=multi_params_adjustable_init)), 
                                ui.panel_main(ui.h4("Multifit initial values"),
                                                ui.row(ui.column(4, ui.input_numeric("mint", "Intensity, 1e-25 cm/mol", value=33.)),
                                                    ui.column(4, ui.input_numeric("melow", "Lower level energy, 1/cm", value=0.)),
                                                    ui.column(4, ui.input_numeric("mmolm", "Molecular mass, a.m.u.", value=28.))
                                                    ),
                                                ui.row(ui.column(4, ui.input_numeric("mf0", "Central frq", value=115271.202)),
                                                    ui.column(4, ui.input_numeric("mfrab", "Rabi frq", value=0.1))),
                                                ui.row(ui.column(3, ui.input_numeric("mg0s", "Gamma0 self", value=3.375)),
                                                    ui.column(3, ui.input_numeric("mng0s", "n Gamma0 self", value=0.78)),
                                                    ui.column(3, ui.input_numeric("mg2s", "Gamma2 self", value=0.33)),
                                                    ui.column(3, ui.input_numeric("mng2s", "n Gamma2 self", value=0.78))),
                                                ui.row(ui.column(3, ui.input_numeric("mg0f", "Gamma0 foreign", value=3.33)),
                                                    ui.column(3, ui.input_numeric("mng0f", "n Gamma0 foreign", value=0.78)),
                                                    ui.column(3, ui.input_numeric("mg2f", "Gamma2 foreign", value=0.33)),
                                                    ui.column(3, ui.input_numeric("mng2f", "n Gamma2 foreign", value=0.78))),
                                                ui.row(ui.column(3, ui.input_numeric("md0s", "Delta0 self", value=-0.004)),
                                                    ui.column(3, ui.input_numeric("mnd0s", "n Delta0 self", value=0.78)),
                                                    ui.column(3, ui.input_numeric("md2s", "Delta2 self", value=0.00)),
                                                    ui.column(3, ui.input_numeric("mnd2s", "n Delta2 self", value=0.78))),
                                                ui.row(ui.column(3, ui.input_numeric("md0f", "Delta0 foreign", value=-0.004)),
                                                    ui.column(3, ui.input_numeric("mnd0f", "n Delta0 foreign", value=0.78)),
                                                    ui.column(3, ui.input_numeric("md2f", "Delta2 foreign", value=0.00)),
                                                    ui.column(3, ui.input_numeric("mnd2f", "n Delta2 foreign", value=0.78))),
                                                ui.row(ui.column(3, ui.input_numeric("my0s", "Mixing self", value=0.000006)),
                                                    ui.column(3, ui.input_numeric("mny0s", "n mixing self", value=0.78)),
                                                    ui.column(3, ui.input_numeric("my0f", "Mixing foreign", value=0.000006)),
                                                    ui.column(3, ui.input_numeric("mny0f", "n mixing foreign", value=0.78))),
                                                ui.row(ui.column(3, ui.input_numeric("mnuvcs", "Vel.chn. rate self", value=0.15)),
                                                    ui.column(3, ui.input_numeric("mnnuvcs", "n nu_vc self", value=0.78)),
                                                    ui.column(3, ui.input_numeric("mnuvcf", "Vel.chn. rate foreign", value=0.15)),
                                                    ui.column(3, ui.input_numeric("mnnuvcf", "n nu_vc foreign", value=0.78))),
                                                ui.row(ui.column(3, ui.input_numeric("mcs", "Continuum self", value=5.e-18)),
                                                    ui.column(3, ui.input_numeric("mncs", "n C self", value=0.0)),
                                                    ui.column(3, ui.input_numeric("mcf", "Continuum foreign", value=5.e-18)),
                                                    ui.column(3, ui.input_numeric("mncf", "n C foreign", value=0.0))),
                                                ui.row(ui.column(6, ui.input_numeric("mpow", "power factor", value=0.)),
                                                    ui.column(6, ui.input_numeric("mscl", "scale factor", value=1.0)),),
                                                ui.row(ui.column(3, ui.input_numeric("mbl0", "Baseline 0", value=0.0)),
                                                    ui.column(3, ui.input_numeric("mbl1", "Baseline 1", value=0.0)),
                                                    ui.column(3, ui.input_numeric("mbl2", "Baseline 2", value=0.0)),
                                                    ui.column(3, ui.input_numeric("mbl3", "Baseline 3", value=0.0))),
                                                )
                                ),
                ui.input_action_button("b_preview_multifit", "Preview multifit"),
                ui.output_plot("preview_spectrum_multifit"),
                ui.row(ui.column(4, ui.input_switch('s_verbose', "Verbose fit")),
                       ui.column(4, ui.input_action_button("b_fit_multifit", "Run multifit"))),
                ui.row(ui.column(12, ui.output_text("fit_status_multifit"))),
                ui.row(ui.column(4, ui.input_switch("s_f0_offset_multifit", "Subtract central frequency")), 
                       ui.column(4, ui.input_numeric("f0_shift_multifit", "Central frq", value=115271.202))
                       ),
                ui.output_plot("preview_multifit"),
                ui.row(ui.column(3, ui.download_button("download_params_multifit", "Download results")),
                       ui.column(3, ui.download_button("download_resid_multifit", "Download residuals")),
                       ui.column(3, ui.input_text("resid_prefix_multifit", "Residual file prefix", value="rm_"))
                       ),
                ui.output_ui("multifit_result_table")
                               ),
   
                )
)


def server(input, output, session):

    # some global variables which are visible from any part of the server code

    residuals: reactive.Value[list] = reactive.Value([])
    recording: reactive.Value[list] = reactive.Value([])
    names: reactive.Value[list] = reactive.Value([])
    f_aux = reactive.Value(False)  # flag that aux info loaded
    f_part = reactive.Value(False)  # flag that partition function is loaded
    f_preview = reactive.Value(False)  # flag that spectra are loaded and preview is available
    f_fit = reactive.Value(False)  # flag that fit is finished
    f_fit_multi = reactive.Value(False)  # flag that multifit is finished
    part_data_global: reactive.Value[list] = reactive.Value([])
    params_fit: reactive.Value[list] = reactive.Value([])  # params from fit
    uncert_fit: reactive.Value[list] = reactive.Value([])  # uncertainties
    params_fit_aux: reactive.Value[list] = reactive.Value([])  # aux params (see model function)
    params_fit_multi: reactive.Value[list] = reactive.Value([])
    uncert_fit_multi: reactive.Value[list] = reactive.Value([])
    residuals_multi: reactive.Value[list] = reactive.Value([])    
    nrecs = reactive.Value(0)
    aux_df: pd.DataFrame = reactive.Value()
    text_status_multifit = reactive.Value('no fit result currently')

    # Service part just to test things

    @output
    @render.text
    @reactive.event(input.b_check)
    def jac_flag_out():
        p = ui.Progress(min=0, max=10)
        #with ui.Progress(min=0, max=10) as p:
        for i in range(10):
            p.set(i, message='Test!')
            time.sleep(1)
        p.close()

    ### COMMON SECTION
        
    # Processing of the loaded spectra (common for LBL fit and multifit)
    @reactive.Effect
    @reactive.event(input.input_files)
    def _():
        tmp_recs = []
        tmp_names = []
        f: list[FileInfo] = input.input_files()
        # f has fields: name, size, type, datapath
        for ifil in range(len(f)):
            tmp_names.append(f[ifil]["name"])
            cur_data = np.loadtxt(f[ifil]["datapath"], comments=input.text_comment())
            tmp_recs.append(cur_data)
        recording.set(tmp_recs)
        nrecs.set(len(f))
        names.set(tmp_names)
        if f_aux.get() and f_part.get():
            f_preview.set(True)

    
    # Processing of the loaded metadata on the experimental conditions etc (common)
    @reactive.Effect
    @reactive.event(input.aux_data)
    def _():
        if not (input.aux_data() is None):
            f_aux.set(True)
        if f_aux.get() and f_part.get() and not (input.input_files() is None):
            f_preview.set(True)

    
    # Processing of the loaded partition function (common)
    @reactive.Effect
    @reactive.event(input.partition)
    def _():
        if not (input.partition() is None):
            f_part.set(True)
        if input.s_nopart():
            f_part.set(True)
        if f_aux.get() and f_part.get() and not (input.input_files() is None):
            f_preview.set(True)
        tmp_part_data = np.loadtxt(input.partition()[0]["datapath"],
                                   dtype='float', comments='#', delimiter=None,
                                   unpack=False)  # partition data for intensity
        tmp_part_t = tmp_part_data[:, 0]
        tmp_part_q = tmp_part_data[:, 1]
        out_list = []
        out_list.append(tmp_part_t)
        out_list.append(tmp_part_q)
        part_data_global.set(out_list)


    # Showing loaded metadata on the experimental conditions etc (common)
    @output
    @render.ui
    def aux_data_show():
        if input.aux_data() is None:
            return "Upload a file with recordings conditions. It should contain columns: " \
                   "filename, self pressure (Torr), foreign pressure (Torr), temperature (K)," \
                   "frq deviation (MHz; 0 if not applicable), recording type, cell length (cm; 1 if not applicable). " \
                   "Record types codes: 0 - resonator (deviation = 0), 1 - RAD with freq. manipulation " \
                   "(deviation should be non-zero), " \
                   "2 - video with freq. manipulation (deviation should be non-zero), 3 - video with amplitude" \
                   "modulation (deviation = 0)."
        aux_file: list[FileInfo] = input.aux_data()
        if input.s_header():
            aux_header = 0
            aux_names = None
        else:
            aux_header = None
            aux_names = ['#file', 'Pself', 'Pforeign', 'T', 'dev', 'typ', 'L']
        aux_out = pd.read_csv(aux_file[0]["datapath"], header=aux_header, names=aux_names,
                              delim_whitespace=True, index_col=0)
        aux_df.set(aux_out)
        return ui.HTML(aux_out.to_html(classes="table table-striped"))

    ### LINE-BY-LINE FIT SECTION ###

    # Showing the table of the fitted parameters
    @output
    @render.ui
    def fit_result_table():
        if f_fit.get():
            headers = list(single_params_dict.values())
            headers = ['P self', 'P foreign', 'T'] + headers
            params_show = pd.DataFrame(params_fit.get(), columns=headers)
            return ui.HTML(params_show.to_html(classes="table table-striped"))
    
    
    # this shows the preview of the loaded spectra ("Preview" button)
    @output
    @render.plot(alt="Preview of the loaded spectra")
    @reactive.event(input.b_preview)
    def preview_spectrum():
        if input.input_files() is None:
            fig, ax = plt.subplots()
            ax = plt.text(0.4, 0.4, "Please upload your spectra first")
            return fig
        if input.aux_data() is None:
            fig, ax = plt.subplots()
            ax = plt.text(0.4, 0.4, "Please upload data on the recordings pressures, temperatures, etc")
            return fig
        if (input.partition() is None) and not input.s_nopart():
            fig, ax = plt.subplots()
            ax = plt.text(0.4, 0.4, "Please upload your partition data")
            return fig
        # rearrange information from aux data file in accordance with the
        # order of the names of the loaded spectra
        aux_df.set(aux_df.get().reindex(names.get()))
        # read the list of the loaded filenames etc
        f: list[FileInfo] = input.input_files()
        fig, ax = plt.subplots(1, len(f), sharey=True)
        # print('length of files dataset is ', len(f), ' files')
        # f has fields: name, size, type, datapath
        #aux_df = pd.read_csv(input.aux_data()[0]["datapath"], header=0, delim_whitespace=True, index_col=0)
        aux_data = aux_df.get().to_numpy(dtype=float)
        p_self = aux_data[:, 0]  # self pressure, Torr
        p_for = aux_data[:, 1]  # foreign pressure, Torr
        tmpr = [t+273.5 if abs(t) < 100. else t for t in aux_data[:, 2]]  # temperature, K
        dev = aux_data[:, 3]  # frq deviation where applicable, MHz (for RAD and VID)
        rtype = aux_data[:, 4]  # record type (0=CAV, 1=RAD+dev, 2=VID+dev, 3=VID natural)
        clen = aux_data[:, 5]  # cell length where applicable, cm (for RAD and VIC)

        # read data for partition sum calculations for line strength
        if not(input.s_nopart()):
            part_data = np.loadtxt(input.partition()[0]["datapath"],
                                   dtype='float', comments='#', delimiter=None,
                                   unpack=False)  # partition data for intensity
        #else:
        #    part_data = np.empty((2, 2))
        #    part_data = [[1., 1000.], [1., 1.]]

        for ifil in range(len(f)):
            cur_data = np.loadtxt(f[ifil]["datapath"], comments=input.text_comment())
            if rtype[ifil] == 0:  # cur_data[0, 0] < 1000: # CAV recordings is usually in GHz and 1/cm
                # we convert it to MHz and 1/km
                if cur_data[0, 0] < 1000:
                    cur_data[:, 0] = cur_data[:, 0] * 1000.
                cur_data[:, 1] = cur_data[:, 1] * 1.e+5

            # calculate initial parameters based on the values specified above, pressures and temperatures
            params0 = np.empty(npar)
            aux_params0 = np.empty(nauxpar)
            t_dep = (t_ref / tmpr[ifil]) ** input.ngam()  # T-factor for collisional params
            params0[0] = input.f0() + (input.d0() * p_self[ifil] + input.d0f() * p_for[ifil]) * t_dep  # center
            params0[1] = 1.            
            params0[2] = (input.g0() * p_self[ifil] + input.g0f() * p_for[ifil]) * t_dep  # G0
            params0[3] = (input.g2() * p_self[ifil] + input.g2f() * p_for[ifil]) * t_dep  # G2
            params0[4] = (input.d2() * p_self[ifil] + input.d2f() * p_for[ifil]) * t_dep  # D2
            params0[5] = (input.y0() * p_self[ifil] + input.y0f() * p_for[ifil]) * t_dep  # Y
            params0[6] = (input.nu() * p_self[ifil] + input.nuf() * p_for[ifil])  # Dicke parameter
            params0[7:] = 0.  # some hardware-related params (and non-resonant, if applicable

            # auxilary params
            # statistics-related part of the strength
            partition_factor = 1.
            if not input.s_nopart:
                partition_factor = qpart(tmpr[ifil], t_ref, part_data[:, 0], part_data[:, 1])
            strength = partition_factor \
                       * math.exp(-c2 * input.elow() / tmpr[ifil]) \
                       * (1 - math.exp(-c2 * input.f0() / (tmpr[ifil] * mhz_in_cm))) \
                       / math.exp(-c2 * input.elow() / t_ref) \
                       / (1 - math.exp(-c2 * input.f0() / (t_ref * mhz_in_cm)))
            # molecule and thermodynamics-related factors
            strength *= input.I0() * p_self[ifil] * k_st * 1.e-25 / tmpr[ifil]

            if rtype[ifil] == 0:
                strength *= 1.e05  # 1/cm to 1/km for CAV recordings
            if rtype[ifil] in [1, 2]:
                strength *= clen[ifil]  # multuply to length for better perfomance on the RAD and VID recordings

            aux_params0[0] = strength
            aux_params0[1] = (input.f0() / clight) \
                             * math.sqrt(2 * math.log(2.) * tmpr[ifil] / (input.molm() * k_vs_aem))  # dopler width
            aux_params0[2] = dev[ifil]  # deviation
            aux_params0[3] = clen[ifil]  # cell length
            aux_params0[-1] = rtype[ifil]  # recording type

            tmp_absor = mdl(cur_data[:, 0], params0, aux_params=aux_params0)
            params0[1] = max(cur_data[:, 1])/max(tmp_absor)  # I scale

            # model with initial parameters
            model0 = mdl(cur_data[:, 0], params0, aux_params=aux_params0)

            if nrecs.get() == 1:
                ax.plot(cur_data[:, 0] * 0.001, cur_data[:, 1], 'ro')
                ax.plot(cur_data[:, 0] * 0.001, model0, 'b-')
                ax.plot(cur_data[:, 0] * 0.001, (cur_data[:, 1] - model0) * 10., 'k-')
                ax.text(0.2, 0.8,'P_self = '+str(p_self[ifil]), ha='left', va='center', transform=ax.transAxes)
                ax.text(0.2, 0.7, 'P_foreign = '+str(p_for[ifil]), ha='left', va='center', transform=ax.transAxes)
            else:
                ax[ifil].plot(cur_data[:, 0] * 0.001, cur_data[:, 1], 'ro')
                ax[ifil].plot(cur_data[:, 0] * 0.001, model0, 'b-')
                ax[ifil].plot(cur_data[:, 0] * 0.001, (cur_data[:, 1] - model0) * 10., 'k-')
                ax[ifil].text(0.2, 0.8, 'P_self = '+str(p_self[ifil]), ha='left', va='center', transform=ax[ifil].transAxes)
                ax[ifil].text(0.2, 0.7, 'P_foreign = '+str(p_for[ifil]), ha='left', va='center', transform=ax[ifil].transAxes)
            # ax = plt.subplot(ifil, cur_data[:,0]*0.001, cur_data[:,1]/np.max(cur_data[:,1]))
        # print(fnames)
        # spectra_info()
        return fig

    
    # Showing the status of the data preview (guides user if some data are absent)
    @output
    @render.text
    @reactive.event(input.b_preview)
    def preview_status():
        text_status = ''
        if input.input_files() is None:
            text_status += 'Load spectra. '
        if input.aux_data() is None:
            text_status += 'Load pressures etc. '
        if input.partition() is None:
            text_status += 'Load partition. '
        if len(text_status) == 0:
            text_status = 'All data present, see your preview'
        return text_status

    
    # Run line-by-line fit routine
    @reactive.Effect
    @reactive.event(input.b_fit)  # within this func we get the parameters fitted to each loaded profile
    def _():
        f_fit.set(False)
        aux_df.set(aux_df.get().reindex(names.get()))
        # load the list of files with spectra
        f: list[FileInfo] = input.input_files()
        # f has fields: name, size, type, datapath
        # load aux info and form aux params arrays (pressure, temperature etc)
        #aux_df = pd.read_csv(input.aux_data()[0]["datapath"], header=0, delim_whitespace=True, index_col=0)
        aux_data = aux_df.get().to_numpy(dtype=float)
        p_self = aux_data[:, 0]  # self pressure, Torr
        p_for = aux_data[:, 1]  # foreign pressure, Torr
        tmpr = [t+273.5 if abs(t) < 100. else t for t in aux_data[:, 2]]  # temperature, K
        dev = aux_data[:, 3]  # frq deviation where applicable, MHz (for RAD and VID)
        rtype = aux_data[:, 4]  # record type (0=CAV, 1=RAD+dev, 2=VID+dev, 3=VID natural)
        clen = aux_data[:, 5]  # cell length where applicable, cm (for RAD and VIC)

        params_out = np.empty((len(f), 3 + npar))  # array of params for all the recordings
        uncert_out = np.empty((len(f), 3 + npar))
        params_temp = np.empty(3 + npar)  # array of params for a chosen recording
        uncert_temp = np.empty(3 + npar)

        params_out_aux = np.empty((len(f), nauxpar))  # this array will be transferred to reactive value array

        resid_out = []

        # read data for partition sum calculations for line strength
        if not(input.s_nopart()):
            part_data = np.loadtxt(input.partition()[0]["datapath"],
                                   dtype='float', comments='#', delimiter=None,
                                   unpack=False)  # partition data for intensity
        p = ui.Progress(min=0, max=len(f))
        for ifil in range(len(f)):
            p.set(ifil, message='Fit in progress')
            cur_data = np.loadtxt(f[ifil]["datapath"], comments=input.text_comment())  # recording read from file
            if rtype[ifil] == 0:  # cur_data[0, 0] < 1000: # CAV recordings is usually in GHz and 1/cm
                # we convert it to MHz and 1/km
                if cur_data[0, 0] < 1000:
                    cur_data[:, 0] = cur_data[:, 0] * 1000.
                cur_data[:, 1] = cur_data[:, 1] * 1.e+5

            # calculate initial parameters based on the values specified above, pressures and temperatures
            params0 = np.empty(npar)
            params_aux0 = np.empty(nauxpar)
            t_dep = (t_ref / tmpr[ifil]) ** input.ngam()  # T-factor for collisional params
            params0[0] = input.f0() + (input.d0() * p_self[ifil] + input.d0f() * p_for[ifil]) * t_dep  # center
            params0[1] = 1.  # I scale            
            params0[2] = (input.g0() * p_self[ifil] + input.g0f() * p_for[ifil]) * t_dep  # G0
            params0[3] = (input.g2() * p_self[ifil] + input.g2f() * p_for[ifil]) * t_dep  # G2
            params0[4] = (input.d2() * p_self[ifil] + input.d2f() * p_for[ifil]) * t_dep  # D2
            params0[5] = (input.y0() * p_self[ifil] + input.y0f() * p_for[ifil]) * t_dep  # Y
            params0[6] = (input.nu() * p_self[ifil] + input.nuf() * p_for[ifil]) * t_dep  # Dicke parameter
            params0[7:] = 0.  # some hardware-related params (and non-resonant, if applicable

            # auxilary params
            partition_factor = 1.
            if not input.s_nopart:
                partition_factor = qpart(tmpr[ifil], t_ref, part_data[:, 0], part_data[:, 1])
            strength = partition_factor \
                       * math.exp(-c2 * input.elow() / tmpr[ifil]) \
                       * (1 - math.exp(-c2 * input.f0() / (tmpr[ifil] * mhz_in_cm))) \
                       / math.exp(-c2 * input.elow() / t_ref) \
                       / (1 - math.exp(-c2 * input.f0() / (t_ref * mhz_in_cm)))
            strength *= input.I0() * p_self[ifil] * k_st * 1.e-25 / tmpr[ifil]  # absolute line integral intensity

            if rtype[ifil] == 0:
                strength *= 1.e05  # 1/cm to 1/km for CAV recordings
            if rtype[ifil] in [1, 2]:
                strength *= clen[ifil]  # multuply to length for better perfomance on the RAD and VID recordings

            params_aux0[0] = strength
            params_aux0[1] = (input.f0() / clight) \
                             * math.sqrt(2 * math.log(2.) * tmpr[ifil] / (input.molm() * k_vs_aem))  # dopler width
            params_aux0[2] = dev[ifil]  # deviation
            params_aux0[3] = clen[ifil]  # cell length
            params_aux0[-1] = rtype[ifil]  # recording type

            params_out_aux[ifil, :] = params_aux0[:]

            tmp_absor = mdl(cur_data[:, 0], params0, aux_params=params_aux0)
            params0[1] = max(cur_data[:, 1])/max(tmp_absor)  # I scale

            # array shows which parameters are adjusted and which are not
            jac_flag = [int(elem in input.jac_check()) for elem in single_params_dict]

            # calling the fit subroutine
            params1, _ = fit_params(cur_data[:, 0], cur_data[:, 1], params0, mdl, mdljac, jac_flag,
                                 1.e-4, 1.e-12, 1.e4, 10, aux_params=params_aux0)

            uncert_1 = fit_uncertainties(cur_data[:, 0], cur_data[:, 1],
                                         params1, mdl, mdljac, jac_flag,
                                         aux_params=params_aux0)

            model1 = mdl(cur_data[:, 0], params1, aux_params=params_aux0)

            resid1 = cur_data[:, 1] - model1[:]
            frq1 = np.copy(cur_data[:, 0])
            resid1_to_out = np.stack((frq1, resid1), axis=1)

            resid_out.append(resid1_to_out)

            # qual_fit = int(np.max(cur_data[:, 1]) / np.std(cur_data[:, 1] - model1))
            # ax[ifil].text(0.8, 0.8, 'Q = '+str(qual_fit), ha='center', va='center', transform=ax[ifil].transAxes)

            params_temp[0] = p_self[ifil]
            params_temp[1] = p_for[ifil]
            params_temp[2] = tmpr[ifil]
            params_temp[3:3 + npar] = params1[:]
            params_out[ifil, :] = params_temp[:]

            uncert_temp[0] = p_self[ifil]
            uncert_temp[1] = p_for[ifil]
            uncert_temp[2] = tmpr[ifil]
            uncert_temp[3:3 + npar] = uncert_1[:]
            uncert_out[ifil, :] = uncert_temp[:]

        params_fit.set(params_out)
        params_fit_aux.set(params_out_aux)
        uncert_fit.set(uncert_out)
        residuals.set(resid_out)
        f_fit.set(True)
        p.close()

    
    # Download the results of the line-by-line fit (lineshape parameters)
    @render.download(filename='fit_params.txt')
    def download_params():
        if f_fit.get():
            headers = list(single_params_dict.values())
            headers = [elem.replace(' ', '_') for elem in headers]
            headers_err = [elem + '_err' for elem in headers]
            headers_with_err = []
            for i_head in range(len(headers)):
                headers_with_err.append(headers[i_head])
                headers_with_err.append(headers_err[i_head])
            headers = ['#P_self', 'P_foreign', 'T'] + headers_with_err  # + ['\n']
            headers = [elem.ljust(15) for elem in headers]
            yield ' '.join(headers + ['\n'])
            params_array = np.asarray(params_fit.get())
            uncert_array = np.asarray(uncert_fit.get())
            for i in range(params_array.shape[0]):
                params_row = params_array[i, :]
                uncert_row = uncert_array[i, :]
                output_row = np.empty(3 + 2 * npar)
                output_row[0: 3] = params_row[0: 3]
                output_row[3::2] = params_row[3::]
                output_row[4::2] = uncert_row[3::]
                yield ' '.join([('%15.10g' % elem).ljust(15) for elem in output_row] + ['\n'])

        # headers = list(single_params_dict.values())
        # headers = ['# P self', 'P foreign', 'T'] + headers + ['\n']
        # yield '    '.join(headers)
        # for i in range(5):
        #     tmp_data = [25., 125., 296.5, 3250., 318., 0., 7.975e-6]
        #     yield '    '.join(['%.5g' % elem for elem in tmp_data] + ['\n'] )


    # Download the results of the line-by-line fit (residuals)
    @render.download(filename='residuals.zip')
    def download_resid():
        if f_fit.get():
            with io.BytesIO() as buf:
                # this creates in-memory buffer where zip archive is stored
                # without storing file on the disk
                testzip = zipfile.ZipFile(buf, 'w')  # create a ZipFile object, instead of filename use buffer object name
                tmp_names = names.get()
                tmp_resid = residuals.get()
                for i in range(len(tmp_names)):
                    current_resid_name = input.resid_prefix()+tmp_names[i]
                    current_df = pd.DataFrame(tmp_resid[i])
                    # create a file within zip archive
                    testzip.writestr(current_resid_name, current_df.to_string(header=False, index=False))
                testzip.close()
                yield buf.getvalue()  # return the buffer content as file to download


    # Showing the preview of the residuals of the line-by-line fit
    @output
    @render.plot(alt='Fit residuals')
    def preview_fit():
        if f_fit.get():
            f: list[FileInfo] = input.input_files()
            fig, ax = plt.subplots(1, len(f), sharey=True)
            for ifil in range(len(f)):
                if params_fit_aux.get()[ifil, -1] == 0:
                    tmp_xlabel = 'Frequency detuning, GHz'
                    tmp_frq_scale = 1.E-3
                else:
                    tmp_xlabel = 'Frequency detuning, MHz'
                    tmp_frq_scale = 1.

                q_factor = max(recording.get()[ifil][:, 1]) / np.std(residuals.get()[ifil][:, 1])

                if params_fit_aux.get()[ifil, -1] == 0:
                    q_factor *= 1.e5

                if nrecs.get() == 1:
                    ax.text(0.2, 0.8, 'Q = ' + str(round(q_factor)),
                            ha='center', va='center', transform=ax.transAxes)
                    ax.hlines(0.,
                              min((residuals.get()[ifil][:, 0] - params_fit.get()[ifil, 3]) * tmp_frq_scale),
                              max((residuals.get()[ifil][:, 0] - params_fit.get()[ifil, 3]) * tmp_frq_scale),
                              linestyles='dashed', colors='grey')
                    ax.plot((residuals.get()[ifil][:, 0] - params_fit.get()[ifil, 3]) * tmp_frq_scale,
                            residuals.get()[ifil][:, 1], 'b-')
                else:
                    ax[ifil].text(0.2, 0.8, 'Q = ' + str(round(q_factor)),
                                  ha='center', va='center', transform=ax[ifil].transAxes)

                    ax[ifil].hlines(0.,
                                    min((residuals.get()[ifil][:, 0] - params_fit.get()[ifil, 3]) * tmp_frq_scale),
                                    max((residuals.get()[ifil][:, 0] - params_fit.get()[ifil, 3]) * tmp_frq_scale),
                                    linestyles='dashed', colors='grey')

                    ax[ifil].plot((residuals.get()[ifil][:, 0] - params_fit.get()[ifil, 3]) * tmp_frq_scale,
                                  residuals.get()[ifil][:, 1],
                                  'b-')

                    ax[ifil].set_xlabel(tmp_xlabel)
            return fig
        else:
            fig, ax = plt.subplots()
            ax.text(0.5, 0.5, 'No residuals to show')
            return fig
        

    # Showing plots of fitted parameters-vs-pressure for line-by-line fit 
    # (after pressing "Plot vs pressure" button and setting necessary switches on/off)
    @output
    @render.plot(alt='Parameters vs pressure')
    @reactive.event(input.b_coefs)
    def preview_params():
        if not f_fit.get():
            fig, ax = plt.subplots()
            ax.text(0.5, 0.5, "Press FIT button first to get parameters")
            return fig

        fig, axs = plt.subplots(3, 2, sharey='none', figsize=(30, 90), squeeze=False)
        ax = axs.flat

        p_self = params_fit.get()[:, 0]
        p_foreign = params_fit.get()[:, 1]
        tmpr = params_fit.get()[:, 2]

        f0 = params_fit.get()[:, 3]
        i0 = params_fit.get()[:, 4]
        g0 = params_fit.get()[:, 5]
        g2 = params_fit.get()[:, 6]
        d2 = params_fit.get()[:, 7]
        y0 = params_fit.get()[:, 8]
        cc = params_fit.get()[:, 13]

        f0e = uncert_fit.get()[:, 3] + 1.E-10
        i0e = uncert_fit.get()[:, 4] + 1.E-10
        g0e = uncert_fit.get()[:, 5] + 1.E-10
        g2e = uncert_fit.get()[:, 6] + 1.E-15
        d2e = uncert_fit.get()[:, 7] + 1.E-15
        y0e = uncert_fit.get()[:, 8] + 1.E-15
        cce = uncert_fit.get()[:, 13] + 1.E-25

        points_style = 'ro'
        line_style = 'k-'

        if input.s_normself():
            pnorm = p_foreign[:] / p_self[:]

            gam0_coefs, gam0_cov = np.polyfit(pnorm, g0[:] / p_self[:], deg=1, rcond=None, full=False,
                                              w=1. / g0e[:] ** 2, cov=True)
            gam0_self = [gam0_coefs[1], np.sqrt(gam0_cov[1, 1])]
            gam0_for = [gam0_coefs[0], np.sqrt(gam0_cov[0, 0])]

            gam2_coefs, gam2_cov = np.polyfit(pnorm, g2[:] / p_self[:], 1, rcond=None, full=False, w=1. / g2e[:] ** 2,
                                              cov=True)
            gam2_self = [gam2_coefs[1], np.sqrt(gam2_cov[1, 1])]
            gam2_for = [gam2_coefs[0], np.sqrt(gam2_cov[0, 0])]

            del2_coefs, del2_cov = np.polyfit(pnorm, d2[:] / p_self[:], 1, rcond=None, full=False, w=1. / d2e[:] ** 2,
                                              cov=True)
            del2_self = [del2_coefs[1], np.sqrt(del2_cov[1, 1])]
            del2_for = [del2_coefs[0], np.sqrt(del2_cov[0, 0])]

            y0_coefs, y0_cov = np.polyfit(pnorm, y0[:] / p_self[:], 1, rcond=None, full=False, w=1. / y0e[:] ** 2,
                                          cov=True)
            y0_self = [y0_coefs[1], np.sqrt(y0_cov[1, 1])]
            y0_for = [y0_coefs[0], np.sqrt(y0_cov[0, 0])]

            ind_f0 = 0
            ind_y = 1
            ind_g0 = 2
            ind_g2 = 3
            ind_i0 = 4
            ind_c = 5

            ax[ind_f0].errorbar(pnorm, (f0[:] - input.f0_shift()) / p_self[:],
                                xerr=None, yerr=f0e[:] / p_self[:],
                                fmt=points_style)
            ax[ind_f0].set_xlabel('P_foreign/P_self')
            ax[ind_f0].set_ylabel('Center freq. shift/P_self, MHz/Torr')

            ax[ind_g0].errorbar(pnorm, g0[:] / p_self[:],
                                xerr=None, yerr=g0e[:] / p_self[:],
                                fmt=points_style)
            ax[ind_g0].plot(pnorm, gam0_self[0] + pnorm * gam0_for[0], line_style)
            ax[ind_g0].set_xlabel('P_foreign/P_self')
            ax[ind_g0].set_ylabel('Gamma_0/P_self, MHz/Torr')
            ax[ind_g0].text(0.5, 0.9, 'S: %.3f(%0.f) MHz/Torr' % (gam0_self[0], gam0_self[1] * 1.E3),
                            ha='center', va='center', transform=ax[ind_g0].transAxes)
            ax[ind_g0].text(0.5, 0.8, 'F: %.3f(%0.f) MHz/Torr' % (gam0_for[0], gam0_for[1] * 1.E3),
                            ha='center', va='center', transform=ax[ind_g0].transAxes)

            ax[ind_g2].errorbar(pnorm, g2[:] / p_self[:],
                                xerr=None, yerr=g2e[:] / p_self[:],
                                fmt=points_style)
            ax[ind_g2].plot(pnorm, gam2_self[0] + pnorm * gam2_for[0], line_style)
            ax[ind_g2].set_xlabel('P_foreign/P_self')
            ax[ind_g2].set_ylabel('Gamma_2/P_self, MHz/Torr')
            ax[ind_g2].text(0.5, 0.9, 'S: %.3f(%0.f) MHz/Torr' % (gam2_self[0], gam2_self[1] * 1.E3),
                            ha='center', va='center', transform=ax[ind_g2].transAxes)
            ax[ind_g2].text(0.5, 0.8, 'F: %.3f(%0.f) MHz/Torr' % (gam2_for[0], gam2_for[1] * 1.E3),
                            ha='center', va='center', transform=ax[ind_g2].transAxes)

            ax[ind_y].errorbar(pnorm, y0[:] / p_self[:], xerr=None, yerr=y0e[:] / p_self[:], fmt=points_style)
            ax[ind_y].plot(pnorm, y0_self[0] + pnorm * y0_for[0], line_style)
            ax[ind_y].set_xlabel('P_foreign/P_self')
            ax[ind_y].set_ylabel('Mixing parameter/P_self, 1/Torr')

            ax[ind_y].text(0.5, 0.9, 'S: %.3f(%0.f)*1E-6 1/Torr' % (y0_self[0] * 1E6, y0_self[1] * 1.E9),
                           ha='center', va='center', transform=ax[ind_y].transAxes)
            ax[ind_y].text(0.5, 0.8, 'F: %.3f(%0.f)*1E-6 1/Torr' % (y0_for[0] * 1E6, y0_for[1] * 1.E9),
                           ha='center', va='center', transform=ax[ind_y].transAxes)

            ax[ind_i0].errorbar(p_self[:] / (kB * tmpr[:]), i0, xerr=None, yerr=i0e, fmt=points_style)
            ax[ind_i0].set_xlabel('Absorber concentration, 1/m^3')
            ax[ind_i0].set_ylabel('Intensity correction, unitless')

            ax[ind_c].errorbar(pnorm, cc[:] / p_self[:] ** 2, xerr=None, yerr=cce[:] / p_self[:] ** 2, fmt=points_style)
            ax[ind_c].set_xlabel('P_foreign/P_self')
            ax[ind_c].set_ylabel('Continuum parameter')

        if (not input.s_normself()) and input.s_foreign():

            gam0_coefs, gam0_cov = np.polyfit(p_foreign, g0[:], deg=1, rcond=None, full=False,
                                              w=1. / g0e[:] ** 2, cov=True)
            gam0_i = [gam0_coefs[1], np.sqrt(gam0_cov[1, 1])]
            gam0_b = [gam0_coefs[0], np.sqrt(gam0_cov[0, 0])]

            gam2_coefs, gam2_cov = np.polyfit(p_foreign, g2[:], 1, rcond=None, full=False, w=1. / g2e[:] ** 2,
                                              cov=True)
            gam2_i = [gam2_coefs[1], np.sqrt(gam2_cov[1, 1])]
            gam2_b = [gam2_coefs[0], np.sqrt(gam2_cov[0, 0])]

            del2_coefs, del2_cov = np.polyfit(p_foreign, d2[:], 1, rcond=None, full=False, w=1. / d2e[:] ** 2,
                                              cov=True)
            del2_i = [del2_coefs[1], np.sqrt(del2_cov[1, 1])]
            del2_b = [del2_coefs[0], np.sqrt(del2_cov[0, 0])]

            del0_coefs, del0_cov = np.polyfit(p_foreign, f0[:], deg=1, rcond=None, full=False,
                                              w=1. / f0e[:] ** 2, cov=True)
            del0_i = [del0_coefs[1], np.sqrt(del0_cov[1, 1])]
            del0_b = [del0_coefs[0], np.sqrt(del0_cov[0, 0])]

            ind_g0 = 0
            ind_g2 = 1
            ind_f0 = 2
            ind_d2 = 3
            ind_i0 = 4


            ax[ind_f0].errorbar(p_foreign, f0[:] - del0_i[0],
                                xerr=None, yerr=f0e[:],
                                fmt=points_style)
            ax[ind_f0].plot(p_foreign, p_foreign * del0_b[0], line_style)
            ax[ind_f0].set_xlabel('P_foreign')
            ax[ind_f0].set_ylabel('Center freq. shift, MHz')
            ax[ind_f0].text(0.5, 0.8, 'Foreign shift: %.3f(%0.f) Mhz/Torr' % (del0_b[0], del0_b[1]*1.E3),
                            ha='center', va='center', transform=ax[ind_f0].transAxes)

            ax[ind_g0].errorbar(p_foreign, g0[:],
                                xerr=None, yerr=g0e[:],
                                fmt=points_style)
            ax[ind_g0].plot(p_foreign, gam0_i[0] + p_foreign * gam0_b[0], line_style)
            ax[ind_g0].set_xlabel('P$_{foreign}$')
            ax[ind_g0].set_ylabel('Gamma$_0$, MHz/Torr')
            ax[ind_g0].text(0.5, 0.8, 'Foreign g$_0$: %.3f(%0.f) MHz/Torr' % (gam0_b[0], gam0_b[1] * 1.E3),
                            ha='center', va='center', transform=ax[ind_g0].transAxes)

            ax[ind_g2].errorbar(p_foreign, g2[:],
                                xerr=None, yerr=g2e[:],
                                fmt=points_style)
            ax[ind_g2].plot(p_foreign, gam2_i[0] + p_foreign * gam2_b[0], line_style)
            ax[ind_g2].set_xlabel('P_foreign')
            ax[ind_g2].set_ylabel('Gamma$_2$, MHz/Torr')
            ax[ind_g2].text(0.5, 0.8, 'Foreign $g_2$: %.3f(%0.f) MHz/Torr' % (gam2_b[0], gam2_b[1] * 1.E3),
                            ha='center', va='center', transform=ax[ind_g2].transAxes)
            
            ax[ind_d2].errorbar(p_foreign, d2[:],
                                xerr=None, yerr=d2e[:],
                                fmt=points_style)
            ax[ind_d2].plot(p_foreign, del2_i[0] + p_foreign * del2_b[0], line_style)
            ax[ind_d2].set_xlabel('P_foreign')
            ax[ind_d2].set_ylabel('Delta$_2$, MHz/Torr')
            ax[ind_d2].text(0.5, 0.8, 'Foreign d$_2$: %.3f(%0.f) MHz/Torr' % (del2_b[0], del2_b[1] * 1.E3),
                            ha='center', va='center', transform=ax[ind_d2].transAxes)

            ax[ind_i0].errorbar(p_self[:] / (kB * tmpr[:]), i0, xerr=None, yerr=i0e, fmt=points_style)
            ax[ind_i0].set_xlabel('Absorber concentration, 1/m^3')
            ax[ind_i0].set_ylabel('Intensity correction, unitless')

        if (not input.s_normself()) and (not input.s_foreign()):

            gam0_coefs, gam0_cov = np.polyfit(p_self, g0[:], deg=1, rcond=None, full=False,
                                              w=1. / g0e[:] ** 2, cov=True)
            gam0_i = [gam0_coefs[1], np.sqrt(gam0_cov[1, 1])]
            gam0_b = [gam0_coefs[0], np.sqrt(gam0_cov[0, 0])]

            gam2_coefs, gam2_cov = np.polyfit(p_self, g2[:], 1, rcond=None, full=False, w=1. / g2e[:] ** 2,
                                              cov=True)
            gam2_i = [gam2_coefs[1], np.sqrt(gam2_cov[1, 1])]
            gam2_b = [gam2_coefs[0], np.sqrt(gam2_cov[0, 0])]

            del2_coefs, del2_cov = np.polyfit(p_self, d2[:], 1, rcond=None, full=False, w=1. / d2e[:] ** 2,
                                              cov=True)
            del2_i = [del2_coefs[1], np.sqrt(del2_cov[1, 1])]
            del2_b = [del2_coefs[0], np.sqrt(del2_cov[0, 0])]

            del0_coefs, del0_cov = np.polyfit(p_self, f0[:], deg=1, rcond=None, full=False,
                                              w=1. / f0e[:] ** 2, cov=True)
            del0_i = [del0_coefs[1], np.sqrt(gam0_cov[1, 1])]
            del0_b = [del0_coefs[0], np.sqrt(gam0_cov[0, 0])]

            ind_g0 = 0
            ind_g2 = 1
            ind_f0 = 2
            ind_i0 = 3

            ax[ind_f0].errorbar(p_self, f0[:] - del0_i[0],
                                xerr=None, yerr=f0e[:],
                                fmt=points_style)
            ax[ind_f0].plot(p_self, p_self * del0_b[0], line_style)
            ax[ind_f0].set_xlabel('P_self')
            ax[ind_f0].set_ylabel('Center freq. shift, MHz')
            ax[ind_f0].text(0.5, 0.8, 'Self shift: %.3f(%0.f) Mhz/Torr' % (del0_b[0], del0_b[1]*1.E3),
                            ha='center', va='center', transform=ax[ind_f0].transAxes)

            ax[ind_g0].errorbar(p_self, g0[:],
                                xerr=None, yerr=g0e[:],
                                fmt=points_style)
            ax[ind_g0].plot(p_self, gam0_i[0] + p_self * gam0_b[0], line_style)
            ax[ind_g0].set_xlabel('P_self')
            ax[ind_g0].set_ylabel('Gamma_0, MHz/Torr')
            ax[ind_g0].text(0.5, 0.8, 'Self g0: %.3f(%0.f) MHz/Torr' % (gam0_b[0], gam0_b[1] * 1.E3),
                            ha='center', va='center', transform=ax[ind_g0].transAxes)

            ax[ind_g2].errorbar(p_self, g2[:],
                                xerr=None, yerr=g2e[:],
                                fmt=points_style)
            ax[ind_g2].plot(p_self, gam2_i[0] + p_self * gam2_b[0], line_style)
            ax[ind_g2].set_xlabel('P_self')
            ax[ind_g2].set_ylabel('Gamma_2, MHz/Torr')
            ax[ind_g2].text(0.5, 0.8, 'Self g2: %.3f(%0.f) MHz/Torr' % (gam2_b[0], gam2_b[1] * 1.E3),
                            ha='center', va='center', transform=ax[ind_g2].transAxes)

            ax[ind_i0].errorbar(p_self[:] / (kB * tmpr[:]), i0, xerr=None, yerr=i0e, fmt=points_style)
            ax[ind_i0].set_xlabel('Absorber concentration, 1/m^3')
            ax[ind_i0].set_ylabel('Intensity correction, unitless')


        return fig

    ### MULTIFIT SECTION ###

    # this shows the preview of the loaded spectra and multifit model funtion with initial parameters
    # when "Preview multifit button" is pressed
    @output
    @render.plot(alt="Preview of the loaded spectra")
    @reactive.event(input.b_preview_multifit)
    def preview_spectrum_multifit():
        if input.input_files() is None:
             fig, ax = plt.subplots()
             ax = plt.text(0.4, 0.4, "Please upload your spectra first")
             return fig
        if input.aux_data() is None:
            fig, ax = plt.subplots()
            ax = plt.text(0.4, 0.4, "Please upload data on the recordings pressures, temperatures, etc")
            return fig
        if (input.partition() is None) and not input.s_nopart():
            fig, ax = plt.subplots()
            ax = plt.text(0.4, 0.4, "Please upload your partition data")
            return fig
        # print('Preview start')
        aux_df.set(aux_df.get().reindex(names.get()))  # set order of the lines in metadata in accordance with the filenames order
        # read the list of the loaded filenames etc
        f: list[FileInfo] = input.input_files()
        
        aux_data = aux_df.get().to_numpy(dtype=float)
        p_self = aux_data[:, 0]  # self pressure, Torr
        p_for = aux_data[:, 1]  # foreign pressure, Torr
        tmpr = [t+273.5 if abs(t) < 100. else t for t in aux_data[:, 2]]  # temperature, K
        dev = aux_data[:, 3]  # frq deviation where applicable, MHz (for RAD and VID)
        rtype = aux_data[:, 4]  # record type (0=CAV, 1=RAD+dev, 2=VID+dev, 3=VID natural)
        clen = aux_data[:, 5]  # cell length where applicable, cm (for RAD and VID)

        tmp_recs = [np.copy(rec) for rec in recording.get()]        

        for ifil in range(len(f)):            
            rec_cur = tmp_recs[ifil]  # recording
            f_cur = rec_cur[:, 0]  # frequency
            s_cur = rec_cur[:, 1]  # signal
            fr_factor = 1.  # to recalc GHz to MHz when necessary

            aux_cur = np.zeros((rec_cur.shape[0], aux_data.shape[1]+3))            
            for i_aux in range(aux_data.shape[1]):
                aux_cur[:, i_aux] = aux_data[ifil, i_aux]

            aux_cur[:, 2] = [t+273.5 if abs(t) < 100. else t for t in aux_cur[:, 2]]

            aux_cur[:, -1] = ifil            
                
            partition_factor = 1.
            if not input.s_nopart:
                part_t, part_q = part_data_global.get()
                partition_factor = qpart(tmpr[ifil], t_ref, part_t, part_q)
            strength = partition_factor \
                       * math.exp(-c2 * input.melow() / tmpr[ifil]) \
                       * (1 - math.exp(-c2 * input.mf0() / (tmpr[ifil] * mhz_in_cm))) \
                       / math.exp(-c2 * input.melow() / t_ref) \
                       / (1 - math.exp(-c2 * input.mf0() / (t_ref * mhz_in_cm)))
            # molecule and thermodynamics-related factors
            strength *= input.mint() * p_self[ifil] * k_st * 1.e-25 / tmpr[ifil]
            aux_cur[:, -2] = strength

            aux_cur[:, -3] = (input.mf0() / clight) \
                       * math.sqrt(2 * math.log(2.) * tmpr[ifil] / (input.mmolm() * k_vs_aem))

            if rtype[ifil] == 0:
                s_cur[:] *= 1.E5
            # if rtype[ifil] in [1, 2, 3]:
            #     tmp_scales[ifil] = strength * clen[ifil] / max(s_cur)

            if f_cur[0] < 1000.:
                fr_factor = 1000.            
            if ifil == 0:
                frqs = np.copy(f_cur)
                sgnl = np.copy(s_cur)
                aux_list = np.copy(aux_cur)
            else:
                frqs = np.concatenate((frqs, f_cur))
                sgnl = np.concatenate((sgnl, s_cur))
                aux_list = np.concatenate((aux_list, aux_cur))

        frqs[:] *= fr_factor                
        # pnum = np.arange(len(frqs))

        mnpar = n_const_par + len(f) * n_add_par  # number of params for multifit
        params = [0.] * mnpar # np.zeros(mnpar)
        params[multi_params_indx['mint']] = 1.
        params[multi_params_indx['mf0']] = input.mf0()
        params[multi_params_indx['mg0s']] = float(input.mg0s())
        params[multi_params_indx['mg0f']] = input.mg0f()
        params[multi_params_indx['mg2s']] = input.mg2s()
        params[multi_params_indx['mg2f']] = input.mg2f()
        params[multi_params_indx['md0s']] = input.md0s()
        params[multi_params_indx['md0f']] = input.md0f()
        params[multi_params_indx['md2s']] = input.md2s()
        params[multi_params_indx['md2f']] = input.md2f()
        params[multi_params_indx['my0s']] = input.my0s()
        params[multi_params_indx['my0f']] = input.my0f()
        params[multi_params_indx['mnuvcs']] = input.mnuvcs()
        params[multi_params_indx['mnuvcf']] = input.mnuvcf()
        params[multi_params_indx['mcs']] = input.mcs()
        params[multi_params_indx['mcf']] = input.mcf()
        params[multi_params_indx['mfrab']] = input.mfrab()
        params[multi_params_indx['mng0s']] = input.mng0s()
        params[multi_params_indx['mng0f']] = input.mng0f()
        params[multi_params_indx['mng2s']] = input.mng2s()
        params[multi_params_indx['mng2f']] = input.mng2f()
        params[multi_params_indx['mnd0s']] = input.mnd0s()
        params[multi_params_indx['mnd0f']] = input.mnd0f()
        params[multi_params_indx['mnd2s']] = input.mnd2s()
        params[multi_params_indx['mnd2f']] = input.mnd2f()
        params[multi_params_indx['mny0s']] = input.mny0s()
        params[multi_params_indx['mny0f']] = input.mny0f()
        params[multi_params_indx['mnnuvcs']] = input.mnnuvcs()
        params[multi_params_indx['mnnuvcf']] = input.mnnuvcf()
        params[multi_params_indx['mncs']] = input.mncs()
        params[multi_params_indx['mncf']] = input.mncf()
        
        for ifil in range(len(f)):
            #tmp_where = aux_list[:, -1] == ifil
            istart = n_const_par + ifil * n_add_par
            params[istart] = input.mscl()  # integral intensity correction
            params[istart+1] = input.mpow()  # radiation source power vs frequency correction
            
            params[istart+2] = input.mbl0()  # bl0
            params[istart+3] = input.mbl1()  # bl1
            params[istart+4] = input.mbl2()  # bl2
            params[istart+5] = input.mbl3()  # bl3
        
        params = np.asarray(params)
        # print(params)

        model0 = mdl_multi(frqs, params, aux_list)

        for ifil in range(len(f)):                        
            tmp_where = aux_list[:, -1] == ifil
            m_cur = model0[tmp_where]
            s_cur = sgnl[tmp_where]
            istart = n_const_par + ifil * n_add_par
            params[istart] = np.max(s_cur) * input.mscl() / np.max(m_cur)
            tmp_bl0 = 0.
            npt_bl = 5
            for i_bl in range (npt_bl):
                tmp_bl0 += s_cur[i_bl] + s_cur[-1 - i_bl]
            tmp_bl0 *= 1./(2. * npt_bl)
            # print(tmp_bl0)
            params[istart+2] = tmp_bl0  # bl0
            # print('Scale', ifil, params[istart])
        
        model0 = mdl_multi(frqs, params, aux_list)
        
        # fig, ax = plt.subplots()
        fig, ax = plt.subplots(1, len(f), sharey=True)
        for ifil in range(len(f)):
            tmp_where = aux_list[:, -1] == ifil
            f_cur = frqs[tmp_where]
            s_cur = sgnl[tmp_where]
            m_cur = model0[tmp_where]
            ax[ifil].plot(f_cur - input.mf0(), s_cur, 'ro')
            ax[ifil].plot(f_cur - input.mf0(), m_cur, 'b-')
            ax[ifil].plot(f_cur - input.mf0(), (s_cur - m_cur)*10., 'k-')
            # ax[ifil].set_xlabel('Frequency, GHz')
            ax[ifil].text(0.2, 0.8, 'P_self = '+str(p_self[ifil]), ha='left', va='center', transform=ax[ifil].transAxes)
            ax[ifil].text(0.2, 0.7, 'P_foreign = '+str(p_for[ifil]), ha='left', va='center', transform=ax[ifil].transAxes)
        # print('Preview finish')
        return fig

    
    # Preview of the multifit residuals
    @output
    @render.plot(alt="Multifit residuals")
    def preview_multifit():
        if f_fit_multi.get():
            fig, ax = plt.subplots(1, nrecs.get(), sharey=True)
            for ifil in range(nrecs.get()):
                f_cur = recording.get()[ifil][:, 0]
                r_cur = residuals_multi.get()[ifil]                
                if input.s_f0_offset_multifit():
                    frq0 = input.f0_shift_multifit()
                else:
                    frq0 = 0.
                q_factor = max(recording.get()[ifil][:, 1]) / np.std(r_cur)
                
                ax[ifil].plot(f_cur - frq0, r_cur, 'b-')
                ax[ifil].text(0.2, 0.8, 'Q = ' + str(round(q_factor)), ha='left', va='center', transform=ax[ifil].transAxes)
                ax[ifil].hlines(y=0., 
                                xmin=min(f_cur - frq0),
                                xmax=max(f_cur - frq0),
                                linestyles='dashed', colors='grey')
        else:
            fig, ax = plt.subplots()
            ax.text(0.5, 0.5, 'No residuals to show', transform=ax.transAxes)
            return fig

    
    # Processing the "Run multifit" button pressed    
    @reactive.Effect
    @reactive.event(input.b_fit_multifit)  
    def _():
        f_fit_multi.set(False)        
        aux_df.set(aux_df.get().reindex(names.get()))
        # read the list of the loaded filenames etc
        f: list[FileInfo] = input.input_files()        
        aux_data = aux_df.get().to_numpy(dtype=float)
        p_self = aux_data[:, 0]  # self pressure, Torr
        p_for = aux_data[:, 1]  # foreign pressure, Torr
        tmpr = [t+273.5 if abs(t) < 100. else t for t in aux_data[:, 2]]  # temperature, K
        dev = aux_data[:, 3]  # frq deviation where applicable, MHz (for RAD and VID)
        rtype = aux_data[:, 4]  # record type (0=CAV, 1=RAD+dev, 2=VID+dev, 3=VID natural)
        clen = aux_data[:, 5]  # cell length where applicable, cm (for RAD and VID)

        tmp_recs = [np.copy(rec) for rec in recording.get()]
        
        # prepare data for multifit processing
        for ifil in range(len(f)):            
            rec_cur = tmp_recs[ifil]  # get recording number ifil
            f_cur = rec_cur[:, 0]  # frequency
            s_cur = rec_cur[:, 1]  # signal
            fr_factor = 1.  # to recalc GHz to MHz when necessary

            aux_cur = np.zeros((rec_cur.shape[0], aux_data.shape[1]+3))  # array of aux data        
            for i_aux in range(aux_data.shape[1]):
                aux_cur[:, i_aux] = aux_data[ifil, i_aux]

            aux_cur[:, 2] = [t+273.5 if abs(t) < 100. else t for t in aux_cur[:, 2]]  # recalc temperature in aux data

            aux_cur[:, -1] = ifil  # the last col of aux_data is number of recording, it is necessary as all data are stacked together

            # calculate integral intensity with respect to partition etc
            # currently the program works with only one line profile in the recording, when expanding it to arbitrary number of lines, smth to be done here
            
            partition_factor = 1.
            if not input.s_nopart:
                part_t, part_q = part_data_global.get()  # partition function (info from HITRAN is recommended)
                partition_factor = qpart(tmpr[ifil], t_ref, part_t, part_q)  # there is an option for absence of partition
            strength = partition_factor \
                       * math.exp(-c2 * input.melow() / tmpr[ifil]) \
                       * (1 - math.exp(-c2 * input.mf0() / (tmpr[ifil] * mhz_in_cm))) \
                       / math.exp(-c2 * input.melow() / t_ref) \
                       / (1 - math.exp(-c2 * input.mf0() / (t_ref * mhz_in_cm)))
            # molecule and thermodynamics-related factors
            strength *= input.mint() * p_self[ifil] * k_st * 1.e-25 / tmpr[ifil]  # the intensity at p_self
            aux_cur[:, -2] = strength

            aux_cur[:, -3] = (input.mf0() / clight) \
                       * math.sqrt(2 * math.log(2.) * tmpr[ifil] / (input.mmolm() * k_vs_aem))

            s_factor = 1.
            if rtype[ifil] == 0:
                s_factor = 1.E5
            # if rtype[ifil] in [1, 2, 3]:
            #     tmp_scales[ifil] = strength * clen[ifil] / max(s_cur)

            if f_cur[0] < 1000.:
                fr_factor = 1000.            
            if ifil == 0:
                frqs = np.copy(f_cur)
                sgnl = np.copy(s_cur * s_factor)
                aux_list = np.copy(aux_cur)
            else:
                frqs = np.concatenate((frqs, f_cur))
                sgnl = np.concatenate((sgnl, s_cur * s_factor))
                aux_list = np.concatenate((aux_list, aux_cur))

        frqs[:] *= fr_factor                
        # pnum = np.arange(len(frqs))

        mnpar = n_const_par + len(f) * n_add_par  # number of params for multifit
        params0 = [0.] * mnpar # np.zeros(mnpar)
        tmp_err = np.zeros(mnpar)
        params0[multi_params_indx['mint']] = 1.
        params0[multi_params_indx['mf0']] = input.mf0()
        params0[multi_params_indx['mg0s']] = input.mg0s()
        params0[multi_params_indx['mg0f']] = input.mg0f()
        params0[multi_params_indx['mg2s']] = input.mg2s()
        params0[multi_params_indx['mg2f']] = input.mg2f()
        params0[multi_params_indx['md0s']] = input.md0s()
        params0[multi_params_indx['md0f']] = input.md0f()
        params0[multi_params_indx['md2s']] = input.md2s()
        params0[multi_params_indx['md2f']] = input.md2f()
        params0[multi_params_indx['my0s']] = input.my0s()
        params0[multi_params_indx['my0f']] = input.my0f()
        params0[multi_params_indx['mnuvcs']] = input.mnuvcs()
        params0[multi_params_indx['mnuvcf']] = input.mnuvcf()
        params0[multi_params_indx['mcs']] = input.mcs()
        params0[multi_params_indx['mcf']] = input.mcf()
        params0[multi_params_indx['mfrab']] = input.mfrab()
        params0[multi_params_indx['mng0s']] = input.mng0s()
        params0[multi_params_indx['mng0f']] = input.mng0f()
        params0[multi_params_indx['mng2s']] = input.mng2s()
        params0[multi_params_indx['mng2f']] = input.mng2f()
        params0[multi_params_indx['mnd0s']] = input.mnd0s()
        params0[multi_params_indx['mnd0f']] = input.mnd0f()
        params0[multi_params_indx['mnd2s']] = input.mnd2s()
        params0[multi_params_indx['mnd2f']] = input.mnd2f()
        params0[multi_params_indx['mny0s']] = input.mny0s()
        params0[multi_params_indx['mny0f']] = input.mny0f()
        params0[multi_params_indx['mnnuvcs']] = input.mnnuvcs()
        params0[multi_params_indx['mnnuvcf']] = input.mnnuvcf()
        params0[multi_params_indx['mncs']] = input.mncs()
        params0[multi_params_indx['mncf']] = input.mncf()        
        
        for ifil in range(nrecs.get()):
            istart = n_const_par + ifil * n_add_par
            
            params0[istart] = 1.  # integral intensity correction
            params0[istart+1] = 0.  # radiation source power vs frequency correction            
            params0[istart+2] = 0.  # bl0
            params0[istart+3] = 0.  # bl1
            params0[istart+4] = 0.  # bl2
            params0[istart+5] = 0.  # bl3
        
        params0 = np.asarray(params0)         

        model0 = mdl_multi(frqs, params0, aux_list)
        # first time we ran the model with scale factors equal to 1.
        # then we correct them based on the experimental data and calcs with scl=1.
        for ifil in range(nrecs.get()):                        
            tmp_where = aux_list[:, -1] == ifil            
            m_cur = model0[tmp_where]
            s_cur = sgnl[tmp_where]
            istart = n_const_par + ifil * n_add_par
            params0[istart] = np.max(s_cur)/np.max(m_cur)  # here the scale is fixed to the proper value
            tmp_bl0 = 0.
            npt_bl = 5
            if rtype[ifil] != 0:
                for i_bl in range (npt_bl):
                    tmp_bl0 += s_cur[i_bl] + s_cur[-1 - i_bl]
                tmp_bl0 *= 1./(2. * npt_bl)
                params0[istart+2] = tmp_bl0  # bl0

        jac_flag_multi = [int(elem in input.jac_check_multi()) for elem in multi_params_dict]

        jac_flag_part1 = jac_flag_multi[0:n_const_par]
        jac_flag_part2 = jac_flag_multi[n_const_par:]
        jac_flag_multi_uncert = jac_flag_part1 + nrecs.get() * jac_flag_part2
        # print(jac_flag_multi)              
        # print(jac_flag_multi_uncert)

        params1, str_status = fit_params(frqs, sgnl, params0, mdl_multi, mdljac_multi, jac_flag_multi, 
                                         1.e-4, 1.e-12, 1.e4, 10, aux_params=aux_list, f_verbose_fit=input.s_verbose())
        
        text_status_multifit.set(str_status)  # write the status string to the global variable

        uncert_1 = fit_uncertainties(frqs, sgnl,
                                     params1, mdl_multi, mdljac_multi, jac_flag_multi_uncert,
                                         aux_params=aux_list)


        model1 = mdl_multi(frqs, params1, aux_list)

        tmp_resid_multi = []
        for ifil in range(nrecs.get()):
            tmp_where = aux_list[:, -1] == ifil
            m_cur = model1[tmp_where]
            s_cur = sgnl[tmp_where]
            tmp_resid = s_cur - m_cur
            tmp_resid_multi.append(tmp_resid)

        residuals_multi.set(tmp_resid_multi)
        params_fit_multi.set(params1)
        uncert_fit_multi.set(uncert_1)

        f_fit_multi.set(True)

    # Showing the status of the latest multifit
    @output
    @render.text
    def fit_status_multifit():
        if f_fit_multi.get():
            return text_status_multifit.get()


    # Showing the results of the multifit routine
    @output
    @render.ui
    def multifit_result_table():
        if f_fit_multi.get():
            headers = ['Value', 'Error', 'Parameter']
            col_params = params_fit_multi.get()
            col_errors = uncert_fit_multi.get()
            col_headers_list = list(multi_params_dict.values())
            col_comments = col_headers_list[0:n_const_par]
            # multi_params_dict.values()[0:n_const_par]
            col_comments_add = col_headers_list[n_const_par:]            
            # multi_params_dict.values()[n_const_par:]
            nfil = len(recording.get())
            for ifil in range(nfil):
                col_comments = col_comments + col_comments_add            
            col_params = pd.Series(col_params)
            col_errors = pd.Series(col_errors)
            col_comments = pd.Series(col_comments)
            params_show = pd.DataFrame({headers[0]: col_params, 
                                        headers[1]: col_errors, 
                                        headers[2]: col_comments})
            return ui.HTML(params_show.to_html(classes="table table-striped"))

    # Download the results of the multifit (lineshape parameters)
    @render.download(filename='fit_params_multi.txt')
    def download_params_multifit():
        if f_fit_multi.get():
            col_headers_list = list(multi_params_dict.values())
            col_comments = col_headers_list[0:n_const_par]
            col_comments_add = col_headers_list[n_const_par:]
            nfil = len(recording.get())
            for ifil in range(nfil):
                col_comments = col_comments + col_comments_add
            headers = ('#Parameter', 'Uncertainty', 'Description')
            headers = [elem.replace(' ', '_') for elem in headers]
            yield '   '.join([('%15s' % headers[0]).ljust(15), 
                              ('%15s' % headers[1]).ljust(15), 
                              ('%s' % headers[2]).ljust(15)] + ['\n'])  # this puts line of headers into the output file
            params_array = np.asarray(params_fit_multi.get())
            uncert_array = np.asarray(uncert_fit_multi.get())
            for i in range(params_array.shape[0]):                                
                yield '   '.join([('%15.10g' % params_array[i]).ljust(15), 
                                  ('%15.10g' % uncert_array[i]).ljust(15),
                                  ('%s' % '#' + col_comments[i]).ljust(15)] + ['\n'])  # this puts lines with parameters, errors and descrpition
                
    @render.download(filename='residuals_multi.zip')
    def download_resid_multifit():
        if f_fit_multi.get():
            with io.BytesIO() as buf:
                # this creates in-memory buffer where zip archive is stored
                # without storing file on the disk
                testzip = zipfile.ZipFile(buf, 'w')  # create a ZipFile object, instead of filename use buffer object name
                tmp_names = names.get()
                tmp_resid = residuals_multi.get()
                for i in range(len(tmp_names)):
                    current_resid_name = input.resid_prefix_multifit()+tmp_names[i]
                    current_df = pd.DataFrame(tmp_resid[i])
                    # create a file within zip archive
                    testzip.writestr(current_resid_name, current_df.to_string(header=False, index=False))
                testzip.close()
                yield buf.getvalue()  # return the buffer content as file to download
    
    # @render_widget
    # def spectra_info():
    #     if input.input_files is None:
    #         return
    #     f: list[FileInfo] = input.input_files()
    #     fnames = list()
    #     for item in f:
    #         fnames.append(item["name"])
    #     for fname in fnames:
    #         items = [widgets.Label(fname),
    #                  widgets.FloatText(value=0., description='P self, Torr'),
    #                  widgets.FloatText(value=0., description='P foreign, Torr'),
    #                  widgets.FloatText(value=0., description='T, K'),
    #                  widgets.Dropdown(options=["VID", "RAD", "CAV"], description="Type")]
    #         yield widgets.HBox(items)


app = App(app_ui, server)
