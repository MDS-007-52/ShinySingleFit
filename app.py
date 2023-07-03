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
    ui.h2("This is my shiny Line-by-line fit!"),
    ui.layout_sidebar(
        ui.panel_sidebar(ui.input_file("aux_data", "Load recordings info", accept=["*.*"], multiple=False)),
        ui.panel_main(ui.output_ui("aux_data_show"))

    ),
#    ui.output_ui("aux_data_table"),
    ui.input_file("input_files", "Load your spectra", accept=["*.*"], multiple=True),
    ui.input_action_button("b_preview", "Preview"),
    ui.output_plot("preview_data")
)

#

def server(input, output, session):
    @output
    @render.ui
    def aux_data_show():
        if input.aux_data() is None:
            return "Upload a file"
        aux_file: list[FileInfo] = input.aux_data()
        print(aux_file[0]["datapath"])
        aux_out = pd.read_csv(aux_file[0]["datapath"], header=0, delim_whitespace=True, index_col=0)
        print(aux_out.columns)
        return ui.HTML(aux_out.to_html(classes = "table table-striped"))

    @output
    @render.plot(alt="Preview of the loaded spectra")
    @reactive.event(input.b_preview)
    def preview_data():
        if input.input_files() is None:
            fig, ax = plt.subplots()
            ax = plt.text(0.4, 0.4, "Please upload your spectra first")
            return fig
        f: list[FileInfo] = input.input_files()
        fig, ax = plt.subplots(1, len(f), sharey=True)
        print('length of files dataset is ', len(f), ' files')
        # f has fields: name, size, type, datapath
        # global fnames
        # fnames = list()
        # for item in f:
        #     fnames.append(item["name"])
        for ifil in range(len(f)):
            cur_data = np.loadtxt(f[ifil]["datapath"], comments='//')
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
