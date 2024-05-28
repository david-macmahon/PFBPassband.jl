# PFB Response

This repository contains a Jupyter notebook that describes how to calculate the
frequency response of a *polyphase filterbank* (PFB).  Along the way it
describes the theory of operation of a PFB, how the PFB frequency response is
observed in practice, and how aliasing plays a role in the PFB response.

After that the notebook presents a detailed step by step walkthrough of the
analytic computaion of the response of a specific PFB in use at the Green Bank
Telescope.  To demonstrate the correctness of the analysis, the computed
response is used to correct the passband shape of a single coarse channel of
actual GBT data included in this repository as a test case.

The [notebook](https://github.com/david-macmahon/PFBResponse.jl/blob/main/pfbresponse.ipynb)
might be viewable directly on GitHub.
