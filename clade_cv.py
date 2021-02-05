from sklearn.model_selection import BaseCrossValidator
from collections.abc import Iterable
import warnings
from itertools import chain, combinations
from math import ceil, floor
import numbers
from abc import ABCMeta, abstractmethod
from inspect import signature

import numpy as np
from scipy.special import comb


from sklearn.utils import indexable, check_random_state, _safe_indexing
from sklearn.utils import _approximate_mode
from sklearn.utils.validation import _num_samples, column_or_1d
from sklearn.utils.validation import check_array
from sklearn.utils.validation import _deprecate_positional_args
from sklearn.utils.multiclass import type_of_target

__all__ = ['LeaveOneCladeOut']

class LeaveOneCladeOut(BaseCrossValidator):
    """Leave one taxid labeled clade out.
    One taxid is drawn from both binary classes.
    """
    def __init__(self, n_splits = 7):
        self.n_splits = n_splits

    def _iter_test_masks(self, X, y, groups, target):
        if groups is None:
            raise ValueError("The 'groups' parameter should not be None.")
        # We make a copy of groups to avoid side-effects during iteration
        groups = check_array(groups, copy=True, ensure_2d=False, dtype=None)
        unique_groups = np.unique(groups[y == target])
        for i in unique_groups:
            yield groups == i

    def get_n_splits(self, X=None, y=None, groups=None):
        """Returns the number of splitting iterations in the cross-validator
        Parameters
        ----------
        X : object
            Always ignored, exists for compatibility.
        y : object
            Always ignored, exists for compatibility.
        groups : array-like of shape (n_samples,)
            Group labels for the samples used while splitting the dataset into
            train/test set. This 'groups' parameter must always be specified to
            calculate the number of splits, though the other parameters can be
            omitted.
        Returns
        -------
        n_splits : int
            Returns the number of splitting iterations in the cross-validator.
        """
        if groups is None:
            raise ValueError("The 'groups' parameter should not be None.")
        groups = check_array(groups, ensure_2d=False, dtype=None)
        return self.n_splits

    def split(self, X, y, groups=None):
        """Generate indices to split data into training and test set.
        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Training data, where n_samples is the number of samples
            and n_features is the number of features.
        y : array-like of shape (n_samples,), default=None
            The target variable for supervised learning problems.
        groups : array-like of shape (n_samples,)
            Group labels for the samples used while splitting the dataset into
            train/test set.
        Yields
        ------
        train : ndarray
            The training set indices for that split.
        test : ndarray
            The testing set indices for that split.
        """
        X, y, groups = indexable(X, y, groups)
        indices = np.arange(_num_samples(X))
        train_targets = {t:[] for t in np.unique(y)}
        test_targets = {t:[] for t in np.unique(y)}
        for target in targets.keys():
            for test_i in self._iter_test_masks(X, y, groups, target):
                train_targets[target].append(indices[np.logical_not(test_i)])
                test_targets[target].append(indices[test_i])
        train_index = np.concat(tuple(map(np.array,train_targets.values())))
        test_index = np.concat(tuple(map(np.array,test_targets.values())))
        yield train_index, test_index