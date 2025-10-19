#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np

def chi2_rd_prior(rdrag_model, rdrag_mean, rdrag_sigma):
    """Gaussian prior on r_d (in Mpc)."""
    if rdrag_sigma <= 0:
        return 0.0
    return ((rdrag_model - rdrag_mean) / rdrag_sigma) ** 2