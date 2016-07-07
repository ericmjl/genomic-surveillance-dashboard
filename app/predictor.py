"""
A microservice that provides the interface for making predictions.

The front-end interface looks something like this:
1. The user selects a given protein and a single drug.
1. With one click of a button:
    1. makes prediction of the sequence pasted in, using the appropriate model.
"""
from flask import Flask, render_template

predictor = Flask(__name__)


@predictor.route('/')
def home():
    return render_template('predictor/index.html')


@predictor.route('/predict')
def predict():
    return render_template('predictor/predict.html')

if __name__ == '__main__':
    predictor.run(debug=True, host='0.0.0.0', port=5550)
