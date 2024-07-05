module PFBPassband

using FFTW

export AbstractPolyphaseFilterbank, CasperPolyphaseFilterbank, coefs, passband

abstract type AbstractPolyphaseFilterbank end

"""
    struct CasperPolyphaseFilterbank <: AbstractPolyphaseFilterbank
        nchan::Int
        ntaps::Int
        window::Function
        lpf::Function
        bug::Bool
    end

Structure to hold parameters of a CASPER polyphase filterbank.  `nchan` is the
total number of frequency channels, including redundant chanels for negative
frequencies that may be omitted in the output of a real input (i.e. non-complex
input) PFB implementaion.  `ntaps` is the number of "taps" (overlapped windows)
of the PFB implementation.  `coeffs` is a `Vector` of length `nchan * ntaps`
containing the coefficients of the PFB filter.  The coefficients are often
computed by the constructor rather being explicitly specified by by the user.
"""
struct CasperPolyphaseFilterbank <: AbstractPolyphaseFilterbank
    nchan::Int
    ntaps::Int
    window::Function
    lpf::Function
    bug::Bool
end

"""
    CasperPolyphaseFilterbank(nchan, ntaps, [window[, lpf[, bug]]])
    CasperPolyphaseFilterbank(nchan, ntaps; [window,] [lpf,] [bug])
    CasperPolyphaseFilterbank(; nchan, ntaps, [window,] [lpf,] [bug])

Constructs a PolyphaseFilterbank object with windowing function window (defaults
to `hamming`) and low pass filter function `lpf` (defaults to `sinc`).
`casperbug` should be given as `true` if the PFB implementation being modeled
contains a specific bug from the CASPER PFB implementation, otherwise it should
be `false` (its default value).

`nchan` is the total number of frequency channels, including redundant chanels
for negative frequencies that may be omitted by real input (i.e. non-complex
input) PFB implementaion.  `ntaps` is the number of "taps" (overlapped windows)
of the PFB implementation.
"""
function CasperPolyphaseFilterbank(nchan, ntaps; window=hamming, lpf=sinc, bug=false)
    CasperPolyphaseFilterbank(nchan, ntaps, window, lpf, bug)
end

function CasperPolyphaseFilterbank(; nchan, ntaps, window=hamming, lpf=sinc, bug=false)
    CasperPolyphaseFilterbank(nchan, ntaps, window, lpf, bug)
end

function Base.show(io::IO, pfb::CasperPolyphaseFilterbank)
    print(io, "CasperPolyphaseFilterbank(;")
    print(io, "nchan=$(pfb.nchan),")
    print(io, "ntaps=$(pfb.ntaps),")
    print(io, "window=$(pfb.window),")
    print(io, "lpf=$(pfb.lpf),")
    print(io, "bug=$(pfb.bug))")
end

"""
    hamming(n)

Built-in `hamming` function.  For more control, use `DSP.Windows.hamming`.
"""
hamming(n) = 0.54 .+ 0.46 .* cospi.(range(-1, stop=1, length=n))

"""
    hanning(n)

Built-in `hanning` function.  For more control, use `DSP.Windows.hanning`.
"""
hanning(n) = 0.5 .* (1 .+ cospi.(range(-1, stop=1, length=n)))

function coefs(pfb::CasperPolyphaseFilterbank)
    n = pfb.nchan * pfb.ntaps
    window = pfb.window(n)
    # LPF function (exclude half-step offset when pfb.bug is true)
    lpf = pfb.lpf.(((0:n-1).+0.5*(1-pfb.bug)) ./ pfb.nchan .- (pfb.ntaps/2))
    window .* lpf
end

"""
    passband(h::AbstractVector, nchan, nfine)

Compute the passband response of a channel of a PFB of `nchan` channels whose
filter coefficients are given in `h`.  The response is computed at `nfine`
points evenly spaced across the coarse channel.
"""
function passband(h::AbstractVector, nchan, nfine)
    # Make sure nfine is large enough
    nchan * nfine > length(h) || error("nfine is too small")

    # Make matrix zeros(nchan, nfine)
    m = zeros(nchan, nfine)

    # Copy h to m (fills the first `ntaps` fine channels)
    copyto!(m, h)

    # Compute 2D real FFT of m
    fm = rfft(m)

    # Response is magnitude squared of fm summed along the first dimension.
    # Take adjoint and sum along second dimension to make result Nx1 rather than
    # 1xN.
    response = dropdims(sum(abs2, fm', dims=2); dims=2)
    # Don't forget to include the symmetrically redundant parts!
    response[1] += sum(abs2, @view(fm[2:end-1,1]))
    response[2:end] .+= reverse(sum(abs2, (@view(fm[2:end-1,2:end]))'; dims=2); dims=1)

    # Normalize to DC channel
    response ./= response[1]

    # Put in proper (i.e. "unshifted") order
    circshift!(response, nfine÷2)
end

"""
    passband(pfb::CasperPolyphaseFilterbank, nfine)

Compute the passband response of a channel of `pfb`.  The response is computed
at `nfine` points evenly spaced across the PFB channel.
"""
function passband(pfb::CasperPolyphaseFilterbank, nfine)
    passband(coefs(pfb), pfb.nchan, nfine)
end

"""
    get_num_threads()

Get the current number of threads that FFTW will utilize.
"""
get_num_threads() = FFTW.get_num_threads()

"""
    set_num_threads(num_threads=Sys.CPU_THREADS)

Set the current number of threads that FFTW will utilize.  If `num_threads` is
not specified, it will default to `Sys.CPU_THREADS`.
"""
function set_num_threads(num_threads::Integer=Sys.CPU_THREADS)
    FFTW.set_num_threads(num_threads)
end

end # module PFBPassband