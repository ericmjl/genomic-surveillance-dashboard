"""
A microservice that provides the interface for making predictions.

The front-end interface looks something like this:
1. The user selects a given protein and a single drug.
1. With one click of a button:
    1. makes prediction of the sequence pasted in, using the appropriate model.
"""
from flask import Flask, render_template, request
from gsdash.sequence_transformer import to_numeric_rep, standardize_sequence
from sklearn.externals import joblib
from sklearn.externals.joblib import Parallel, delayed
from bokeh.charts import Bar, BoxPlot
from bokeh.resources import INLINE
from bokeh.embed import components
from bokeh.models import HoverTool, ResetTool, WheelZoomTool, PanTool, SaveTool
from bokeh.plotting import figure
from gsdash.predutils import predictions, load_model
from gsdash.bokehutils import yerrorbars

import numpy as np


drugs = ['FPV', 'ATV', 'IDV', 'LPV', 'NFV', 'SQV', 'TPV', 'DRV',]

predictor = Flask(__name__)


models_and_drugs = Parallel(n_jobs=-1)(delayed(load_model)(drug)
                                       for drug in drugs)

models = [mdl for mdl, drug in models_and_drugs]
drugs = [drug for mdl, drug in models_and_drugs]

@predictor.route('/')
def home():
    return render_template('predictor/index.html')

@predictor.route('/predict', methods=['POST'])
def predict():
    input_sequence = request.form['sequence']
    seq = standardize_sequence(to_numeric_rep(input_sequence, 'mw'),
                               'protease').reshape(1, -1)


    preds = predictions(drugs, models, seq)

    TOOLS = [PanTool(), ResetTool(), WheelZoomTool(), SaveTool()]

    plot = BoxPlot(data=preds, values='log10(DR)', label='drug',
                   color='drug',
                   title="protease drug resistance", plot_width=600,
                   plot_height=400, legend=False, tools=TOOLS)

    js_resources = INLINE.render_js()
    css_resources = INLINE.render_css()
    script, div = components(plot, INLINE)

    return render_template('predictor/predictions.html', plot_script=script,
                           plot_div=div, js_resources=js_resources,
                           css_resources=css_resources,)

if __name__ == '__main__':
    predictor.run(debug=True, host='0.0.0.0', port=5550)
