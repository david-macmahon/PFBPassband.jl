module CasperPFBs

const ATA1K = CasperPolyphaseFilterbank(;
    nchan=2^11, ntaps=4, width=1.0, window=hamming, lpf=sinc, bug=true
)

# COSMIC uses same PFB as ATA
const COSMIC1K = ATA1K

const GBT512 = CasperPolyphaseFilterbank(;
    nchan=2^10, ntaps=12, width=1.0, window=hamming, lpf=sinc, bug=true
)

const MEERKAT1K = CasperPolyphaseFilterbank(;
    nchan=2^11, ntaps=16, width=0.91, window=hanning, lpf=sinc, bug=false
)

const MEERKAT4K = CasperPolyphaseFilterbank(;
    nchan=2^13, ntaps=16, width=1.00, window=hanning, lpf=sinc, bug=false
)

const MEERKAT32K = CasperPolyphaseFilterbank(;
    nchan=2^16, ntaps= 4, width=1.00, window=hanning, lpf=sinc, bug=false
)

end # module CasperPFBs
