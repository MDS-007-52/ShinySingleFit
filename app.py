from shiny.types import FileInfo
from shiny import * #App, Inputs, Outputs, Session, render, ui
from shinywidgets import *
import ipywidgets as widgets
import ipysheet
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

from SDRP import SDRP
from fitmodels import *
from subr_fit import fit_params, fit_uncertainties
from constants import *

global aux_out
aux_out = pd.DataFrame([0., 0., 0., 0., 0, 1.], [1., 1., 1., 1., 0, 1.])

app_ui = ui.page_fluid(
    ui.head_content(
        ui.tags.script(
            src="https://polyfill.io/v3/polyfill.min.js?features=es6"
        ),
        ui.tags.script(
            "if (window.MathJax) MathJax.Hub.Queue(['Typeset', MathJax.Hub]);"
        ),
    ),
    ui.panel_title("This is my shiny Line-by-line fit!"),
    # this sidebar layout for marking which parameters to adjust and initial parameters values
    ui.layout_sidebar(
        ui.panel_sidebar(ui.h4("Mark adjustable"),
            ui.input_checkbox_group("jac_check", "",
                                    single_params_dict, selected=list(single_params_dict.keys())),
                        ui.input_action_button("b_check", "Check"),
                         ui.output_text("jac_flag_out")
                        ),
        ui.panel_main(ui.h4("Parameters and coefficients for initial values"),
            ui.row(
                ui.column(3,  ui.input_numeric("f0", "Line center, MHz", value=115271.)),
                ui.column(4, ui.input_numeric("I0", "Intensity, 1e-25 cm/mol", value=33.)),
                ui.column(4, ui.input_numeric("elow", "Lower level energy, 1/cm", value=0.))
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
    ui.output_text("preview_status")
)

#

def server(input, output, session):
    # this shows the recordings markup info
    # aux_out = pd.DataFrame()

    @output
    @render.ui
    def aux_data_show():
        if input.aux_data() is None:
            return "Upload a file with recordings conditions"
        aux_file: list[FileInfo] = input.aux_data()
        #print(aux_file[0]["datapath"])
        aux_out = pd.read_csv(aux_file[0]["datapath"], header=0, delim_whitespace=True, index_col=0)
        #print(aux_out.columns)
        return ui.HTML(aux_out.to_html(classes="table table-striped"))

    # check literally anything which is going inside my code
    @output
    @render.text
    @reactive.event(input.b_check)
    def jac_flag_out():
        #print([int(elem in input.jac_check()) for elem in single_params_dict])
        #return [int(elem in input.jac_check()) for elem in single_params_dict]
        #aux_out = pd.read_csv(input.aux_data()[0]["datapath"], header=0, delim_whitespace=True, index_col=0)
        return aux_out.to_numpy()[:, 0]

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
        print('length of files dataset is ', len(f), ' files')
        # f has fields: name, size, type, datapath
        aux_df = pd.read_csv(input.aux_data()[0]["datapath"], header=0, delim_whitespace=True, index_col=0)
        aux_data = aux_df.to_numpy(dtype=float)
        p_self = aux_data[:, 0]
        p_for = aux_data[:, 1]
        tmpr = aux_data[:, 2]
        dev = aux_data[:, 3]
        rtype = aux_data[:, 4]
        clen = aux_data[:, 5]
        for ifil in range(len(f)):
            cur_data = np.loadtxt(f[ifil]["datapath"], comments=input.text_comment())
            if rtype[ifil] == 0: #cur_data[0, 0] < 1000:
                cur_data[:, 0] = cur_data[:, 0] * 1000.
            ax[ifil].plot(cur_data[:,0]*0.001, cur_data[:,1]/np.max(cur_data[:,1]), 'ro')
            ax[ifil].text(0.8, 0.8, str(p_self[ifil]), ha='center', va='center', transform=ax[ifil].transAxes)
            #ax = plt.subplot(ifil, cur_data[:,0]*0.001, cur_data[:,1]/np.max(cur_data[:,1]))
        #print(fnames)
        #spectra_info()
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
        if len(text_status)==0:
            text_status = 'All data present, see your preview'
        return text_status




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
