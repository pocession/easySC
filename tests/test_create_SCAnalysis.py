import argparse
from pathlib import Path
import inspect

from src.easySC import SCAnalysis
def test_unittest():

    class Args():
        def __init__(self):
            self.data = '../data'

    args = list()
    exp = SCAnalysis(args)
    
    assert hasattr(exp, 'args')
    assert hasattr(exp, 'data')