from shiny.types import FileInfo
from shiny import *  # App, Inputs, Outputs, Session, render, ui
# from shinywidgets import *
# import ipywidgets as widgets
# import ipysheet
import numpy as np
from pathlib import Path
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
    ui.panel_title("This is my shiny Line-by-line fit!"),
    # this sidebar layout for marking which parameters to adjust and initial parameters values
    ui.layout_sidebar(
        ui.panel_sidebar(ui.h4("Mark adjustable"),
                         ui.input_checkbox_group("jac_check", "",
                                                 single_params_dict, selected=list(single_params_dict.keys())),
                         ui.input_action_button("b_check", "Check"),
                         ui.output_text_verbatim("jac_flag_out", placeholder=True)
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
                              ui.column(3, ui.input_numeric("g0", "g0 self, MHz/Torr", value=3.375)),
                              ui.column(3, ui.input_numeric("g0f", "foreign, MHz/Torr", value=2.750))
                          ),
                          ui.row(
                              ui.column(3, ui.input_numeric("g2", "g2 self, MHz/Torr", value=0.33)),
                              ui.column(3, ui.input_numeric("g2f", "foreign, MHz/Torr", value=0.28))
                          ),
                          ui.row(
                              ui.column(3, ui.input_numeric("d2", "d2 self, MHz/Torr", value=0.00)),
                              ui.column(3, ui.input_numeric("d2f", "foreign, MHz/Torr", value=0.00))
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

    # this layout is for uploading "markup" for recordings, i.e. pressures, temperature, recording type etc

    ui.layout_sidebar(
        ui.panel_sidebar(ui.input_file("aux_data", "Load recordings info", accept=["*.*"], multiple=False),
                         ui.input_file("partition", "Load partition function", accept=["*.*"], multiple=False)),
        ui.panel_main(ui.output_ui("aux_data_show"))

    ),

    # an interface for loading recordings
    ui.row(
        ui.column(4, ui.input_file("input_files", "Load your spectra", accept=["*.*"], multiple=True)),
        ui.column(4, ui.input_text("text_comment", "Commented lines", value="//")),
        ui.column(4, ui.input_action_button("b_preview", "Preview"))
    ),
    ui.output_plot("preview_data"),
    ui.output_text("preview_status"),

    # an interface for fitting the model to the spectra and visualisation of the result

    ui.input_action_button("b_fit", "Fit model to recordings"),
    ui.output_text("fit_status"),
    ui.output_plot("preview_fit"),
    ui.output_ui("fit_result_table"),
    ui.row(ui.column(3, ui.input_switch('s_foreign', "Foreign pressure")),
           ui.column(3, ui.input_switch('s_normself', "Norm. to self P")),
           ui.column(4, ui.input_numeric('f0_shift', "Unshifted center", value=115271.202)),
           ui.column(3, ui.input_action_button("b_coefs", "Plot vs pressure")),
           ui.column(3, ui.download_button("download_params", "Download results"))
           ),

    ui.output_plot("preview_params", height="1200px")
)


def server(input, output, session):

    # some global variables which are visible from any part of the server code

    residuals: reactive.Value[list] = reactive.Value([])
    recording: reactive.Value[list] = reactive.Value([])
    f_aux = reactive.Value(False)  # flag that aux info loaded
    f_part = reactive.Value(False)  # flag that partition function is loaded
    f_preview = reactive.Value(False)  # flag that spectra are loaded and preview is available
    f_fit = reactive.Value(False)  # flag that fit is finished
    params_fit: reactive.Value[list] = reactive.Value([])  # params from fit
    uncert_fit: reactive.Value[list] = reactive.Value([])  # uncertainties
    params_fit_aux: reactive.Value[list] = reactive.Value([])  # aux params (see model function)

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

    @reactive.Effect
    @reactive.event(input.input_files)
    def _():
        tmp_recs = []
        f: list[FileInfo] = input.input_files()
        for ifil in range(len(f)):
            cur_data = np.loadtxt(f[ifil]["datapath"], comments=input.text_comment())
            tmp_recs.append(cur_data)
        recording.set(tmp_recs)
        if f_aux.get() and f_part.get():
            f_preview.set(True)

    @reactive.Effect
    @reactive.event(input.aux_data)
    def _():
        if not (input.aux_data() is None):
            f_aux.set(True)
        if f_aux.get() and f_part.get() and not (input.input_files() is None):
            f_preview.set(True)

    @reactive.Effect
    @reactive.event(input.partition)
    def _():
        if not (input.partition() is None):
            f_part.set(True)
        if f_aux.get() and f_part.get() and not (input.input_files() is None):
            f_preview.set(True)

    # @output
    # @render.text
    # def jac_flag_out():
    #     # return str(f_preview.get()) + str(f_aux.get()) + str(f_part.get())
    #     if f_fit.get():
    #         return str(params_fit.get()[:, 0])

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
        aux_out = pd.read_csv(aux_file[0]["datapath"], header=0, delim_whitespace=True, index_col=0)
        return ui.HTML(aux_out.to_html(classes="table table-striped"))

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
    def preview_data():
        if input.input_files() is None:
            fig, ax = plt.subplots()
            ax = plt.text(0.4, 0.4, "Please upload your spectra first")
            return fig
        if input.aux_data() is None:
            fig, ax = plt.subplots()
            ax = plt.text(0.4, 0.4, "Please upload data on the recordings pressures, temperatures, etc")
            return fig
        if input.partition() is None:
            fig, ax = plt.subplots()
            ax = plt.text(0.4, 0.4, "Please upload your partition data")
            return fig
        # read the list of the loaded filenames etc
        f: list[FileInfo] = input.input_files()
        fig, ax = plt.subplots(1, len(f), sharey=True)
        # print('length of files dataset is ', len(f), ' files')
        # f has fields: name, size, type, datapath
        aux_df = pd.read_csv(input.aux_data()[0]["datapath"], header=0, delim_whitespace=True, index_col=0)
        aux_data = aux_df.to_numpy(dtype=float)
        p_self = aux_data[:, 0]  # self pressure, Torr
        p_for = aux_data[:, 1]  # foreign pressure, Torr
        tmpr = [t+273.5 if abs(t) < 100. else t for t in aux_data[:, 2]]  # temperature, K
        dev = aux_data[:, 3]  # frq deviation where applicable, MHz (for RAD and VID)
        rtype = aux_data[:, 4]  # record type (0=CAV, 1=RAD+dev, 2=VID+dev, 3=VID natural)
        clen = aux_data[:, 5]  # cell length where applicable, cm (for RAD and VIC)

        # read data for partition sum calculations for line strength
        part_data = np.loadtxt(input.partition()[0]["datapath"],
                               dtype='float', comments='#', delimiter=None,
                               unpack=False)  # partition data for intensity

        for ifil in range(len(f)):
            cur_data = np.loadtxt(f[ifil]["datapath"], comments=input.text_comment())
            if rtype[ifil] == 0:  # cur_data[0, 0] < 1000: # CAV recordings is usually in GHz and 1/cm
                # we convert it to MHz and 1/km
                cur_data[:, 0] = cur_data[:, 0] * 1000.
                cur_data[:, 1] = cur_data[:, 1] * 1.e+5

            # calculate initial parameters based on the values specified above, pressures and temperatures
            params0 = np.empty(npar)
            aux_params0 = np.empty(nauxpar)
            params0[0] = input.f0()  # center
            params0[1] = 1.
            t_dep = (t_ref / tmpr[ifil]) ** input.ngam()  # T-factor for collisional params
            params0[2] = (input.g0() * p_self[ifil] + input.g0f() * p_for[ifil]) * t_dep  # G0
            params0[3] = (input.g2() * p_self[ifil] + input.g2f() * p_for[ifil]) * t_dep  # G2
            params0[4] = (input.d2() * p_self[ifil] + input.d2f() * p_for[ifil]) * t_dep  # D2
            params0[5] = (input.y0() * p_self[ifil] + input.y0f() * p_for[ifil]) * t_dep  # Y
            params0[6] = (input.nu() * p_self[ifil] + input.nuf() * p_for[ifil])  # Dicke parameter
            params0[7:] = 0.  # some hardware-related params (and non-resonant, if applicable

            # auxilary params
            # statistics-related part of the strength
            strength = qpart(tmpr[ifil], t_ref, part_data[:, 0], part_data[:, 1]) \
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

            ax[ifil].plot(cur_data[:, 0] * 0.001, cur_data[:, 1], 'ro')
            ax[ifil].plot(cur_data[:, 0] * 0.001, model0, 'b-')
            ax[ifil].text(0.8, 0.8, str(p_self[ifil]), ha='center', va='center', transform=ax[ifil].transAxes)
            # ax = plt.subplot(ifil, cur_data[:,0]*0.001, cur_data[:,1]/np.max(cur_data[:,1]))
        # print(fnames)
        # spectra_info()
        return fig

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

    @reactive.Effect
    @reactive.event(input.b_fit)  # within this func we get the parameters fitted to each loaded profile
    def _():
        f_fit.set(False)
        # load the list of files with spectra
        f: list[FileInfo] = input.input_files()
        # f has fields: name, size, type, datapath
        # load aux info and form aux params arrays (pressure, temperature etc)
        aux_df = pd.read_csv(input.aux_data()[0]["datapath"], header=0, delim_whitespace=True, index_col=0)
        aux_data = aux_df.to_numpy(dtype=float)
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
        part_data = np.loadtxt(input.partition()[0]["datapath"],
                               dtype='float', comments='#', delimiter=None,
                               unpack=False)  # partition data for intensity
        p = ui.Progress(min=0, max=len(f))
        for ifil in range(len(f)):
            p.set(ifil, message='Fit in progress')
            cur_data = np.loadtxt(f[ifil]["datapath"], comments=input.text_comment())  # recording read from file
            if rtype[ifil] == 0:  # cur_data[0, 0] < 1000: # CAV recordings is usually in GHz and 1/cm
                # we convert it to MHz and 1/km
                cur_data[:, 0] = cur_data[:, 0] * 1000.
                cur_data[:, 1] = cur_data[:, 1] * 1.e+5

            # calculate initial parameters based on the values specified above, pressures and temperatures
            params0 = np.empty(npar)
            params_aux0 = np.empty(nauxpar)
            params0[0] = input.f0()  # center
            params0[1] = 1.  # I scale
            t_dep = (t_ref / tmpr[ifil]) ** input.ngam()  # T-factor for collisional params
            params0[2] = (input.g0() * p_self[ifil] + input.g0f() * p_for[ifil]) * t_dep  # G0
            params0[3] = (input.g2() * p_self[ifil] + input.g2f() * p_for[ifil]) * t_dep  # G2
            params0[4] = (input.d2() * p_self[ifil] + input.d2f() * p_for[ifil]) * t_dep  # D2
            params0[5] = (input.y0() * p_self[ifil] + input.y0f() * p_for[ifil]) * t_dep  # Y
            params0[6] = (input.nu() * p_self[ifil] + input.nuf() * p_for[ifil])  # Dicke parameter
            params0[7:] = 0.  # some hardware-related params (and non-resonant, if applicable

            # auxilary params
            strength = qpart(tmpr[ifil], t_ref, part_data[:, 0], part_data[:, 1]) \
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
            params1 = fit_params(cur_data[:, 0], cur_data[:, 1], params0, mdl, mdljac, jac_flag,
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

    @session.download(filename='fit_params.txt')
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
        y0 = params_fit.get()[:, 8]
        cc = params_fit.get()[:, 13]

        f0e = uncert_fit.get()[:, 3] + 1.E-10
        i0e = uncert_fit.get()[:, 4] + 1.E-10
        g0e = uncert_fit.get()[:, 5] + 1.E-10
        g2e = uncert_fit.get()[:, 6] + 1.E-15
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
            gam0_self = [gam0_coefs[1], np.sqrt(gam0_cov[1, 1])]
            gam0_for = [gam0_coefs[0], np.sqrt(gam0_cov[0, 0])]

            gam2_coefs, gam2_cov = np.polyfit(p_foreign, g2[:] / p_self[:], 1, rcond=None, full=False, w=1. / g2e[:] ** 2,
                                              cov=True)
            gam2_self = [gam2_coefs[1], np.sqrt(gam2_cov[1, 1])]
            gam2_for = [gam2_coefs[0], np.sqrt(gam2_cov[0, 0])]

            del0_coefs, del0_cov = np.polyfit(p_foreign, f0[:] / p_self[:], deg=1, rcond=None, full=False,
                                              w=1. / f0e[:] ** 2, cov=True)
            del0_self = [del0_coefs[1], np.sqrt(gam0_cov[1, 1])]
            del0_for = [del0_coefs[0], np.sqrt(gam0_cov[0, 0])]

            ind_g0 = 0
            ind_g2 = 1
            ind_f0 = 2
            ind_i0 = 3

        return fig

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
