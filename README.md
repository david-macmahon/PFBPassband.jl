# PFBPassband.jl

`PFBPassband.jl` is a Julia package for calculating the passband response of a
polyphase filter bank (PFB).  Currently only critically sampled PFBs are
supported.  Any arbitrary PFB can be analyzed by supplying its filter
coefficients, but PFBs designed using the CASPER tool flow can be analyzed
based on just their design parameters.  The filter coefficients for CASPER PFBs
are computed automatically from the design parameters (including a parameter
that indicates whether the CASPER PFB being modeled includes a specific bug in
its implementation).

## CASPER Polyphase Filter Banks

For analysis purposes, a CASPER polyphase filter bank has 5 main design
parameters.  They are:

- Number of channels (`nchan`)
- Number of taps (`ntaps`)
- Windowing function (`window`)
- Low pass filter function (`lpf`)
- Bug flag (whether to model bug) (`bug`)

The `CasperPolyphaseFilterbank` structure contains fields for each of these five
parameters.  Three constructors are provided to allow flexibility in creating
such structures.

```julia
CasperPolyphaseFilterbank(nchan, ntaps, [window[, lpf[, bug]]])
CasperPolyphaseFilterbank(nchan, ntaps; [window,] [lpf,] [bug])
CasperPolyphaseFilterbank(; nchan, ntaps, [window,] [lpf,] [bug])
```

The main difference between these is which parameters are given as positional
arguments vs keyword arguments.  If not given, `window` defaults to `hamming`;
`lpf` defaults to `sinc`; and `bug` defaults to `false`.

`PFBPassband.jl` provides `hamming` and `hanning` window functions, but does not
export them to avoid name clashes with `DSP.jl`, which also provides these
as well as several other window functions suitable for use with
`PFBPassband.jl`.  The default window function is `hamming` (specifically,
`PFBPassband.hamming`).

### Examples

To create a CasperPolyphaseFilterbank with 16 channels and 8 taps (and defaults
for `window`, `lpf` and `bug`):

```julia-repl
julia> using PFBPassband

julia> pfb = CasperPolyphaseFilterbank(16, 8)
CasperPolyphaseFilterbank(;nchan=16,ntaps=8,window=hamming,lpf=sinc,bug=false)
```

## PFB Filter Coefficients

The PFB filter coefficients for a CasperPolyphaseFilterbank can be obtained by
passing it to `coefs`, which will return a `Vector{Float64}` of length
`nchan*ntaps`:

```
coefs(pfb::CasperPolyphaseFilterbank; normalize=true)
coefs(::Type{T}, pfb::CasperPolyphaseFilterbank; normalize=true) where {T<:Base.IEEEFloat}
coefs!(dest::AbstractArray, pfb::CasperPolyphaseFilterbank; normalize=true)
```

By default, the coefficients will be normalized such that their sum is 1.0.
This makes the response of the PFB filter unity at 0 Hz.  To prevent
normalization, pass keyword argument `normalize=false`.

These can be plotted using any Julia plotting package (UnicodePlots shown
here):

```
julia> using PFBPassband, UnicodePlots

julia> pfb = CasperPolyphaseFilterbank(16, 8)
CasperPolyphaseFilterbank(;nchan=16,ntaps=8,window=hamming,lpf=sinc,bug=false)

julia> h = coefs(pfb);

julia> lineplot(h, xlim=(1,length(h)))
         ┌────────────────────────────────────────┐
    0.07 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
         │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡞⢣⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
         │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⢇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
         │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
         │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
         │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡎⠀⠀⠀⠀⢱⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
         │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢠⠇⠀⠀⠀⠀⠸⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
         │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
         │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
         │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢰⠁⠀⠀⠀⠀⠀⠀⠈⡆⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
         │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡜⠀⠀⠀⠀⠀⠀⠀⠀⢣⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
         │⠤⢤⣤⣤⡤⠴⠶⠭⠭⠦⡤⠤⠤⠤⢤⠧⠤⠤⠤⠤⠤⠤⠤⠤⠼⡤⠤⠤⠤⢤⠴⠭⠭⠶⠦⢤⣤⣤⡤⠤│
         │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⢆⠀⠀⡎⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢱⠀⠀⡰⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
         │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠑⠊⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⠊⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
   -0.02 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
         └────────────────────────────────────────┘
         ⠀1⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀128
```

## PFB filter response

The non-aliased (pre-decimation) frequency response of the PFB filter can be
obtained by using `DSP.jl`.  To minimize dependencies, `PFBPassband` does not
depend on `DSP`, so you will need to add `DSP` and `PFBPassband` to your own (or
a temporary) `Project.toml` file.  `DSP.freqresp` computes the frequency
response of a given filter for a range of normalized frequencies, which range
from DC (0 rad/sample) to Nyquist (π rad/sample).  Often the response is only
calculated over a relatively small subset of the entire bandwidth.

Using the PFB from the previous example, here is how to compute and plot the
pre-decimated response over half a coarse channel and all of its neighboring
channel:

```
julia> using DSP

julia> w = range(0, step=1/(32*pfb.nchan), length=32*2) .* pi;

julia> H = freqresp(PolynomialRatio(h, [1]), w);

julia> p=lineplot(10log10.(abs2.(H)), grid=false, ylabel="dBc"); hline!(p, [-6]); vline!(p, [33])
           ┌────────────────────────────────────────┐
        10 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
           │⠠⠒⠒⠒⠒⠒⠒⠒⠒⠒⠒⠤⠤⠤⠤⢄⣀⡀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
           │⠒⠒⠒⠒⠒⠒⠒⠒⠒⠒⠒⠒⠒⠒⠒⠒⠒⠚⢻⢖⡒⠒⠒⠒⠒⠒⠒⠒⠒⠒⠒⠒⠒⠒⠒⠒⠒⠒⠒⠒│
           │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠈⠢⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
           │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠘⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
           │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠘⣆⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
           │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠘⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
   dBc     │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⢱⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
           │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
           │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⡇⡔⠑⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
           │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⢱⡇⠀⠈⢆⠀⠀⠀⠀⠀⠀⠀⠀⠀│
           │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⢸⠇⠀⠀⠈⢆⠀⠀⠀⡠⠊⠀⠀⠀│
           │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠘⠀⠀⠀⠀⠀⠳⡀⢠⠃⠀⠀⠀⠀│
           │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢣⡸⠀⠀⠀⠀⠀│
       -90 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠇⠀⠀⠀⠀⠀│
           └────────────────────────────────────────┘
           ⠀0⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀70
```

## PFB Passband

Use `passband` to get the aliased (post-decimation) frequency response of the
PFB filter coefficients.  This returns the response as normalized power (i.e.
magnitude squared) relative to the center frequency of the coarse channel.
Because it is post-decimation, this includes power that has aliased in from
other channels.

```
julia> pb = passband(pfb, 128);

julia> p = lineplot(-64:63, 10log10.(pb), xlim=(-64,63), ylabel="dBc"); hline!(p, [-3])
          ┌────────────────────────────────────────┐
        1 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
          │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
          │⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣇⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀⣀│
          │⠀⠀⠀⠀⠀⠀⡰⠊⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠢⡀⠀⠀⠀⠀⠀│
          │⠀⠀⠀⠀⠀⡜⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠱⡀⠀⠀⠀⠀│
          │⠀⠀⠀⠀⡸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢱⠀⠀⠀⠀│
          │⠀⠀⠀⢠⠃⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢇⠀⠀⠀│
   dBc    │⠀⠀⠀⡎⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠸⡀⠀⠀│
          │⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢣⠀⠀│
          │⠀⠀⡏⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠘⡄⠀│
          │⠀⡸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢱⠀│
          │⣰⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢣│
          │⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⡏⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉│
          │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
       -4 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
          └────────────────────────────────────────┘
          ⠀-64⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀63
```

## Performance tips

`PFBPassband.jl` uses `FFTW.jl` to perform FFTs.  FFTW can improve throughput
by using multiple threads which takes advantage of multi-core CPUs.  To utilize
this feature, call `PFBPassband.set_num_threads([num_threads])`.  If the
optional `num_threads` argument is not specified it defaults to
`Sys.CPU_THREADS`.  Note that FFTW's worker threads are completely separate
from Julia's thread pool.  FFTW can run multithreaded even when Julia is
started with `--threads=1`.

## Notebook

This repository contains a Jupyter notebook that describes how to calculate the
frequency response of a *polyphase filter bank* (PFB).  Along the way it
describes the theory of operation of a PFB, how the PFB frequency response is
observed in practice, and how aliasing plays a role in the PFB response.

After that the notebook presents a detailed step by step walk through the
analytic computation of the response of a specific PFB in use at the Green Bank
Telescope.  To demonstrate the correctness of the analysis, the computed
response is used to correct the passband shape of a single coarse channel of
actual GBT data included in this repository as a test case.

The [notebook](https://github.com/david-macmahon/PFBPassband.jl/blob/main/notebooks/01_pfb_response.ipynb)
might be viewable directly on GitHub.
