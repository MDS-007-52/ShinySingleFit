from shiny import *
import numpy as np
import matplotlib.pyplot as plt
from SDRP import SDRP
from subr_fit import fit_params, fit_uncertainties
from constants import *

app_ui = ui.page_fluid(
    ui.h2("Hello Shiny!"),
    ui.input_file("input_files", "Load your spectra", accept=["*.*"], multiple=True),
    ui.input_action_button("b_preview", "Preview"),
    ui.output_plot("preview_data"),
)


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
        for ifil in range(len(f)):
            cur_data = np.loadtxt(f[ifil]["datapath"], comments='//')
            ax[ifil].plot(cur_data[:,0]*0.001, cur_data[:,1]/np.max(cur_data[:,1]), 'ro')
            #ax = plt.subplot(ifil, cur_data[:,0]*0.001, cur_data[:,1]/np.max(cur_data[:,1]))
        return fig


app = App(app_ui, server)
