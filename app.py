from shiny import *
from shinywidgets import output_widget, render_widget
import ipywidgets as widgets
import numpy as np
import matplotlib.pyplot as plt

from SDRP import SDRP
from subr_fit import fit_params, fit_uncertainties
from constants import *

app_ui = ui.page_fluid(
    ui.h2("Hello Shiny!"),
    ui.input_file("input_files", "Load your spectra", accept=["*.*"], multiple=True),
    ui.input_action_button("b_preview", "Preview"),
    output_widget("spectra_info"),
    ui.output_plot("preview_data"),
)

#

def server(input, output, session):
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
