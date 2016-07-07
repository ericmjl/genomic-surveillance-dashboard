"""
A microservice that provides the interface for model training.

The front-end interface looks something like this:
1. The user selects a given protein and a single drug.
1. With one click of a button:
    1. a series of ensemble ML models will be trained
    1. the best one will automatically be stored as a Python pickle file in
       the model_store directory.
"""
from flask import Flask, render_template

model_trainer = Flask(__name__)


@model_trainer.route('/')
def home():
    return render_template('model_trainer/index.html')


@model_trainer.route('/train')
def train():
    return render_template('model_trainer/train.html')

if __name__ == '__main__':
    model_trainer.run(debug=True, host='0.0.0.0', port=5550)
