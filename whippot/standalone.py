"""
WHIPPOT from a terminal
"""
import ipywidgets

import shiny.express 
from shinywidgets import render_widget

from whippot import whippot_tools

@render_widget
def return_ui():
    ui = whippot_tools.ComputePositions(initial_values = {}).ui
    return ui

# @render_widget
# def widget():
#     return ipywidgets.IntSlider(1, 1, 100, 1)

