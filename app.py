from shiny.types import FileInfo
from shiny import * #App, Inputs, Outputs, Session, render, ui
from shinywidgets import *
import ipywidgets as widgets
import ipysheet
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd


from SDRP import SDRP
from subr_fit import fit_params, fit_uncertainties
from constants import *

app_ui = ui.page_fluid(
    ui.panel_title("This is my shiny Line-by-line fit!"),
    # this sidebar layout for marking which parameters to adjust and initial parameters values
    ui.layout_sidebar(
        ui.panel_sidebar(ui.row(ui.h4("Mark adjustable")), ui.h4(""),
            ui.input_checkbox_group("jac_check", "",
                         {"frq":"Central freq",
                          "int":"Intensity",
                           "g0":"Gamma 0",
                           "g2":"Gamma 2",
                           "d2":"Delta 2",
                            "y":"Mixing",
                         "nuvc":"Vel.chng.rate",
                          "pow":"Power factor",
                          "bl0":"Baseline 0",
                          "bl1":"Baseline 1",
                          "bl2":"Baseline 2",
                          "bl3":"Baseline 3"}),
                        ui.input_action_button("b_check", "Check")
                        ),
        ui.panel_main(ui.row(ui.h4(" "), ui.h4(" "), ui.h4("Parameters and coefficients for initial values")),
            ui.row(
                ui.column(3,  ui.input_numeric("f0", "Line center, MHz", value=115271.)),
                ui.column(4, ui.input_numeric("I0", "Intensity, 1e-25 cm/mol", value=33.)),
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
        ui.panel_sidebar(ui.input_file("aux_data", "Load recordings info", accept=["*.*"], multiple=False)),
        ui.panel_main(ui.output_ui("aux_data_show"))

    ),

    # an interface for loading recordings
    ui.row(
        ui.column(4, ui.input_file("input_files", "Load your spectra", accept=["*.*"], multiple=True)),
        ui.column(4, ui.input_text("text_comment", "Commented lines", value="//")),
        ui.column(4, ui.input_action_button("b_preview", "Preview"))
        ),
    ui.output_plot("preview_data")
)

#

def server(input, output, session):
    # this shows the recordings markup info
    @output
    @render.ui
    def aux_data_show():
        if input.aux_data() is None:
            return "Upload a file with recordings conditions"
        aux_file: list[FileInfo] = input.aux_data()
        print(aux_file[0]["datapath"])
        aux_out = pd.read_csv(aux_file[0]["datapath"], header=0, delim_whitespace=True, index_col=0)
        print(aux_out.columns)
        return ui.HTML(aux_out.to_html(classes="table table-striped"))

    # this shows the preview of the loaded spectra ("Preview" button)
    @output
    @render.plot(alt="Preview of the loaded spectra")
    @reactive.event(input.b_preview)
    def preview_data():
        if input.input_files() is None:
            fig, ax = plt.subplots()
            ax = plt.text(0.4, 0.4, "Please upload your spectra first")
            return fig
        # read the list of the loaded filenames etc
        f: list[FileInfo] = input.input_files()
        fig, ax = plt.subplots(1, len(f), sharey=True)
        print('length of files dataset is ', len(f), ' files')
        # f has fields: name, size, type, datapath
        # global fnames
        # fnames = list()
        # for item in f:
        #     fnames.append(item["name"])
        for ifil in range(len(f)):
            cur_data = np.loadtxt(f[ifil]["datapath"], comments=input.text_comment())
            if cur_data[0, 0] < 1000:
                cur_data[:, 0] = cur_data[:, 0] * 1000.
            ax[ifil].plot(cur_data[:,0]*0.001, cur_data[:,1]/np.max(cur_data[:,1]), 'ro')
            #ax = plt.subplot(ifil, cur_data[:,0]*0.001, cur_data[:,1]/np.max(cur_data[:,1]))
        #print(fnames)
        #spectra_info()
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
