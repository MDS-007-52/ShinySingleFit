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
        #if input.input_files() is None:
        fig, ax = plt.subplots()
        ax = plt.text(0.4, 0.4, "Please upload your spectra first")
        #    return fig
        #f: list[FileInfo] = input.file1()
        #fig, ax = plt.subplots()
        return fig


app = App(app_ui, server)
