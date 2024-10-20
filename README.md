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
- Normalized channel width (`width`)
- Windowing function (`window`)
- Low pass filter function (`lpf`)
- Bug flag (whether to model bug) (`bug`)

The `CasperPolyphaseFilterbank` structure contains fields for each of these
parameters.  Three constructors are provided to allow flexibility in creating
such structures.

```julia
CasperPolyphaseFilterbank(nchan, ntaps, width, [window[, lpf[, bug]]])
CasperPolyphaseFilterbank(nchan, ntaps; [width,] [window,] [lpf,] [bug])
CasperPolyphaseFilterbank(; nchan, ntaps, [width,] [window,] [lpf,] [bug])
```

The main difference between these is which parameters are given as positional
arguments vs keyword arguments.  If not given, `width` defaults to 1, window`
defaults to `hamming`; `lpf` defaults to `sinc`; and `bug` defaults to `false`.

`PFBPassband.jl` provides `hamming` and `hanning` window functions, but does not
export them to avoid name clashes with `DSP.jl`, which also provides these
as well as several other window functions suitable for use with
`PFBPassband.jl`.  The default window function is `hamming` (specifically,
`PFBPassband.hamming`).

The `width` parameter should be between 0 and 1.  It specifies where the
frequency response will be -6 dB in power relative to the channel width.  A
value of 1 puts this point at the channel edge.  Slightly smaller values can be
used to have the filter roll off start earlier so that the response at the
channel edge will be somewhat reduced.  This reduces aliasing, but at the
expense of attenuating the channel edges.

### Examples

To create a CasperPolyphaseFilterbank with 16 channels and 8 taps (and defaults
for `width`, `window`, `lpf` and `bug`):

```julia-repl
julia> using PFBPassband

julia> pfb = CasperPolyphaseFilterbank(16, 8)
CasperPolyphaseFilterbank(;nchan=16,ntaps=8,width=1.0,window=hamming,lpf=sinc,bug=false)
```

### Pre-defined PFBs

The `CasperFPBs` submodule contains pre-defined `CasperPolyphaseFilterbank`
objects for some CASPER PFB implementations that are known to be in active use.

| name         | nchan | ntaps | width | window  | lpf  |  bug  |
|:-------------|------:|------:|------:|:-------:|:----:|:-----:|
| `ATA1K`      |  2048 |     4 |   1.0 | hamming | sinc | true  |
| `COSMIC1K`   |  2048 |     4 |   1.0 | hamming | sinc | true  |
| `GBT512`     |  1024 |    12 |   1.0 | hamming | sinc | true  |
| `MEERKAT1K`  |  2047 |    16 |  0.91 | hanning | sinc | false |
| `MEERKAT4K`  |  8192 |    16 |  1.00 | hanning | sinc | false |
| `MEERKAT32K` | 65536 |     4 |  1.00 | hanning | sinc | false |

Pull requests are welcome for additional PFB implementations.

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

## Notebooks

This repository contains two Jupyter notebooks in the `notebooks` directory.
They are runnable in any Jupyter installation that has a modern Julia kernel
available, but they are perhaps most easily used from within Visual Studio Code
with the Julia extension.

* [01_pfb_response.ipynb](notebooks/01_pfb_response.ipynb)

  This notebook describes how to compute the frequency response of a
  *polyphase filter bank* (PFB).  Along the way it
  presents the theory of operation of a PFB, how the PFB frequency response is
  observed in practice, and how aliasing plays a role in the PFB response.

  After that the notebook presents a detailed step by step walk through the
  analytic computation of the response of a specific PFB in use at the Green
  Bank Telescope.  To demonstrate the correctness of the analysis, the computed
  response is used to correct the passband shape of a single coarse channel of
  actual GBT data included in this repository as a test case.

  If the above link doesn't work, this [notebook](
  https://github.com/david-macmahon/PFBPassband.jl/blob/main/notebooks/01_pfb_response.ipynb)
  may be viewed directly on GitHub.

* [02_PFBPassband.ipynb](notebooks/02_PFBPassband.ipynb)

  This notebook introduces the `PFBPassband,jl` package.  It covers the
  following topics:

  1. How to use `PFBPassband.jl` to model CASPER polyphase fitlerbanks
  1. How to generate the PFB filter coefficients for a CASPER polyphase filterbank
  1. How to use `DSP.jl` to compute the (unaliased) freqeuncy response of PFB
     filter coefficients
  1. How to use the `passband` function to compute the (aliased) PFB passband
  1. Presents the PFB passband shape of PFBs in active use at several different
     radio telescopes

  If the above link doesn't work, this [notebook](
  https://github.com/david-macmahon/PFBPassband.jl/blob/main/notebooks/02_PFBPassband.ipynb)
  may be viewed directly on GitHub.
