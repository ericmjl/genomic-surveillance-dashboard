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
from bokeh.charts import Bar
from bokeh.resources import INLINE
from bokeh.embed import components
drugs = ['FPV', 'ATV', 'IDV', 'LPV', 'NFV', 'SQV', 'TPV', 'DRV',]

predictor = Flask(__name__)



@predictor.route('/')
def home():
    return render_template('predictor/index.html')


@predictor.route('/predict', methods=['POST'])
def predict():
    input_sequence = request.form['sequence']
    seq = standardize_sequence(to_numeric_rep(input_sequence, 'mw'), 'protease')

    preds = dict()
    preds['drug'] = list()
    preds['log10(DR)'] = list()

    for drug in drugs:
        print(drug)
        mdl = joblib.load("../models/base/{drug}/{drug}.pkl".format(drug=drug))
        pred = mdl.predict(seq)[0]
        # preds.append(dict(drug=drug, pred=pred))
        preds['drug'].append(drug)
        preds['log10(DR)'].append(pred)
        print(pred)

    bar = Bar(data=preds, values='log10(DR)', label='drug', title="protease drug resistance",
              plot_width=600, plot_height=400, legend=False)
    js_resources = INLINE.render_js()
    css_resources = INLINE.render_css()
    script, div = components(bar, INLINE)

    return render_template('predictor/predictions.html', plot_script=script,
                           plot_div=div, js_resources=js_resources,
                           css_resources=css_resources,)

if __name__ == '__main__':
    predictor.run(debug=True, host='0.0.0.0', port=5550)
