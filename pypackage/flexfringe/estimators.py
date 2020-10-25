"""
This is an estimator based on the DFASAT library developed at the University of Luxembourg
"""
import numpy as np
import scipy as sp
import numbers
from sklearn.base import BaseEstimator
import flexfringe.lib.flexfringe as flexfringe
import flexfringe.eval as e
import sys
from multiprocessing import Process, cpu_count
import random
import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

import os,sys,select

# output capture based on
# https://stackoverflow.com/a/9489139
def more_data(elm):
    r, _, _ = select.select([elm],[],[],0)
    return bool(r)

def read_pipe(elm):
    out = ''
    while more_data(elm):
        out += " " + str(os.read(elm, 1024))
    return out

has_graphviz = True
try:
    from graphviz import Source
except:
    has_graphviz = False


class StateMachinePlot():
    """ A helper class implementing a _repr_svg_ method

    Parameters
    ----------
    svg : string
        The plot in SVG format
    """
    _svg = None

    def __init__(self, svg):
        self._svg = svg

    def _repr_svg_(self):
        return self._svg


num_cores = cpu_count()

class BaggingClassifier:
    estimator = None
    number = 10
    estimators = []
    file_prefix = "dfa"

    def __init__(self, estimator=estimator, number=10, output_file="dfa", random_seed=False, random_counts=None, **kwargs):
        self.estimator = estimator
        self.number = number
        self.file_prefix = output_file

        for _ in range(number):
            if random_seed:
                kwargs['seed'] = random.randint(0, 5000000)
            if random_counts is not None:
                kwargs['state_count'] = random.choice(random_counts)
                kwargs['symbol_count'] = random.choice(random_counts)

            self.estimators.append(estimator(**kwargs))

    def calc_model(self, i, estimator, X, y):
        estimator.fit(X, y)
        with open(self.file_prefix + '_' + str(i) + '.model.dot', 'w') as r:
            r.write(estimator.merger.dot_output)


    def fit(self, X, y, subset=True):
        processes = []
        for i, estimator in enumerate(self.estimators):
            if subset:
                print('randomly sampling')
                s_X=[]
                for x in X:
                    if random.randint(1, 5) != 1:
                        s_X.append(x)
                s_y = [1]*len(s_X)
            else:
                s_X = X
                s_y = y

            p = Process(target=self.calc_model, args=(i, estimator, s_X, s_y))
            processes.append(p)
            print('starting process '+str(i))
            p.start()
            if (i+1) % num_cores == 0:
                for p in processes:
                    p.join()
                    processes.remove(p)
        for p in processes:
            p.join()
            processes.remove(p)


    def predict(self, prefix):
        predictions = []
        for estimator in self.estimators:
            m = e.load_model(estimator.merger.dot_output, normalize=True)
            prediction = e.predict(prefix, m)
            predictions.append(prediction)

        sorted_predictions = {}
        for sample in predictions:
            for p in sample:
                if p[0] not in sorted_predictions:
                    sorted_predictions[p[0]] = p[1][1]
                else:
                    sorted_predictions[p[0]] += p[1][1]

        sorted_predictions = sorted(sorted_predictions.items(), key=lambda p: p[1], reverse=True)

        return sorted_predictions

   def predict(self, X):
#        self._verify_x(X)
        """ A reference implementation of a predicting function.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.

        Returns
        -------
        y : array of shape = [n_samples]
            Returns :math:`x^2` where :math:`x` is the first column of `X`.
        """
        if not self.result_dot:
            raise Exception("Run the fit function first!")

        model = e.load_model(self.result_dot)

        return e.predict(X, model)

    def plot(self):
        """ Plot the state machine after fitting

        Returns
        -------
        A StateMachinePlot object implementing a _repr_svg_ method, used
        by Jupyter/IPython notebook to plot a graph
        """
        if has_graphviz:
            if self.result_dot is not None:
                with open("pre_"+self.result_dot) as result_dot_file:
                    return StateMachinePlot(
                       Source(result_dot_file.read())._repr_svg_()
                    )
            else:
                print("Run fit first!")
        else:
            print("This feature requires the graphviz python library")


class flexfringeEstimator(BaseEstimator):
    """ An estimator that uses the flexfringe library

       import code
       code.interact(local=locals())

    """
    result_dot = None
    parameters = None
    output_file = "dfa"
    output_counter = 1

    def __init__(self, **kwargs):
        """hData="likelihood_data", hName="likelihoodratio",
            tries=1, sat_program="", dot_file="dfa", sinkson=1, seed=12345678,
            apta_bound=2000, dfa_bound=50, merge_sinks_p=1, merge_sinks_d=0,
            lower_bound=-1.0, offset=5, method=1, heuristic=1, extend=1, symbol_count=10,
            state_count=25, correction=1.0, parameters=0.5, extra_states=0,
            target_rejecting=0, symmetry=1, forcing=0"""
        
        for key, value in kwargs.items():
            if key == 'output_file':
                self.output_file = value
                continue
            if not hasattr(self.parameters, key):
                raise Exception(key + " not a valid argument")
            setattr(self.parameters, key, value)
            
    def _header(self, X):
        return str(len(X)) + " " + str(self._alphabet_size(X))

    def _alphabet_size(self, X):
        s = set()
        for row in X:
            for value in row[1:]:
                s.add(value.split("/")[0])
        print(s)
        return len(s)
    
    def fit_(self, dfa_file=None, dfa_data=None, output_file="dfa"):
        """A reference implementation of a fitting function
        """
        
    def fit(data):
        return
    
    def predict_(data):
        return
        
    def predict(data):
        return
    
    
