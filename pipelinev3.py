"""Written by ravioli1369@gmail.com"""


"""----------------------------------------------------------------------------------------------------------------"""
import subprocess as sp
import glob
import argparse


"""-----Automatic Pipeline-----------------------------------------------------------------------------------------"""


"""-----functions for automatic pipeline---------------------------------------------------------------------------"""


# cztgtigen
def gtigen(evt, mkf, mkf_threshold):
    sp.call(
        [
            "cztgtigen",
            "eventfile=" + evt,
            "mkffile=" + mkf,
            "thresholdfile=" + mkf_threshold,
            "outfile=" + evt.replace(".evt", ".gti"),
            "usergtifile=-",
            "clobber=YES",
            "history=YES",
        ]
    )


# cztgaas
def gaas(evt, mkf):
    sp.call(
        [
            "cztgaas",
            "par_evtfile=" + evt,
            "par_mkffile=" + mkf,
            "par_outAspectfile=" + evt.replace("_bc.evt", ".aspect"),
            "par_clobber=YES",
            "par_history=YES",
        ]
    )


# cztdatasel
def datasel(evt):
    sp.call(
        [
            "cztdatasel",
            "infile=" + evt,
            "gtitype=QUAD",
            "gtifile=" + evt.replace(".evt", ".gti"),
            "outfile=" + evt.replace("bc.evt", "bc_ds.evt"),
            "clobber=YES",
            "history=YES",
        ]
    )


# cztpixclean
def pixclean(evt, livetime):
    sp.call(
        [
            "cztpixclean",
            "infile=" + evt.replace("bc.evt", "bc_ds.evt"),
            "inlivetimefile=" + livetime,
            "outfile1=" + evt.replace("bc.evt", "quad_pc.evt"),
            "outlivetimefile=" + evt.replace("bc.evt", "quad_livetime.fits"),
            "badpixfile=" + evt.replace("bc.evt", "quad_badpix.fits"),
            "det_tbinsize=1",
            "pix_tbinsize=1",
            "det_count_thresh=100",
            "pix_count_thresh=1000",
        ]
    )


# cztevtclean
def evtclean(evt):
    sp.call(
        [
            "cztevtclean",
            "infile=" + evt.replace("bc.evt", "quad_pc.evt"),
            "outfile=" + evt.replace("bc.evt", "quad_clean.evt"),
            "alphaval=" + "0",
            "vetorange=" + "0-0",
            "clobber=YES",
            "history=YES",
        ]
    )


# cztflagbadpix
def flagbadpix(evt):
    sp.call(
        [
            "cztflagbadpix",
            "nbadpixFiles=1",
            "badpixfile=" + evt.replace("bc.evt", "quad_badpix.fits"),
            "outfile=" + evt.replace("bc.evt", "quad_badpix_out.fits"),
            "clobber=YES",
            "history=YES",
        ]
    )


# cztbindata
def bindata(evt, mkf, timebinsize, emin, emax):
    sp.call(
        [
            "cztbindata",
            "inevtfile=" + evt.replace("bc.evt", "quad_clean.evt"),
            "mkffile=" + mkf,
            "badpixfile=" + evt.replace("bc.evt", "quad_badpix_out.fits"),
            "badpixthreshold=0",
            "livetimefile=" + evt.replace("bc.evt", "quad_livetime.fits"),
            "outputtype=lc",
            "emin=" + emin,
            "emax=" + emax,
            "timebinsize=" + str(timebinsize),
            "outfile=" + evt.replace("bc.evt", "quad_clean"),
            "outevtfile=" + evt.replace("bc.evt", "quad_clean_weight.evt"),
            "maskWeight=n",
            "rasrc=0",
            "decsrc=0",
            "clobber=YES",
            "history=YES",
        ]
    )


# cztdpigen
def dpigen(evt):
    sp.call(
        [
            "cztdpigen",
            "par_infile=" + evt.replace("bc.evt", "quad_clean.evt"),
            "par_badpixFile=" + evt.replace("bc.evt", "quad_badpix_out.fits"),
            "par_badpixThreshold=0",
            "par_outDPHfile=" + evt.replace("bc.evt", "quad_clean.dpi"),
            "par_outDPIfile=" + evt.replace("bc.evt", "quad_clean.dph"),
            "par_quadsToProcess=0,1,2,3",
            "par_timerange=-",
            "par_ebins=-",
            "par_clobber=YES",
            "par_history=YES",
        ]
    )


# cztimage
def image(evt):
    sp.call(
        [
            "cztimage",
            "par_intype=dpi",
            "par_infile=" + evt.replace("bc.evt", "quad_clean.dpi"),
            "par_aspectfileQ0=" + evt.replace("_bc.evt", ".aspect") + "_Q0",
            "par_aspectfileQ1=" + evt.replace("_bc.evt", ".aspect") + "_Q1",
            "par_aspectfileQ2=" + evt.replace("_bc.evt", ".aspect") + "_Q2",
            "par_aspectfileQ3=" + evt.replace("_bc.evt", ".aspect") + "_Q3",
            "par_outImgFile=" + evt.replace("bc.evt", "quad_clean.img"),
            "par_quadsToProcess=0,1,2,3",
            "par_clobber=YES",
            "par_history=YES",
        ]
    )


# cztrspgen
def rspgen(evt):
    quad_clean = evt.replace("bc.evt", "quad_clean")
    quad_cleanevt = evt.replace("bc.evt", "quad_clean.evt")
    quad_badpix_outfits = evt.replace("bc.evt", "quad_badpix_out.fits")
    Q0pha = quad_clean + "_Q0.pha"
    Q1pha = quad_clean + "_Q1.pha"
    Q2pha = quad_clean + "_Q2.pha"
    Q3pha = quad_clean + "_Q3.pha"
    Q0rsp = quad_clean + "_Q0.rsp"
    Q1rsp = quad_clean + "_Q1.rsp"
    Q2rsp = quad_clean + "_Q2.rsp"
    Q3rsp = quad_clean + "_Q3.rsp"
    sp.call(
        [
            "cztrspgen",
            "phafile=" + Q0pha,
            "rspfile=" + Q0rsp,
            "evtfile=" + quad_cleanevt,
            "badpixfile=" + quad_badpix_outfits,
            "clobber=YES",
            "history=YES",
        ]
    )
    sp.call(
        [
            "cztrspgen",
            "phafile=" + Q1pha,
            "rspfile=" + Q1rsp,
            "evtfile=" + quad_cleanevt,
            "badpixfile=" + quad_badpix_outfits,
            "clobber=YES",
            "history=YES",
        ]
    )
    sp.call(
        [
            "cztrspgen",
            "phafile=" + Q2pha,
            "rspfile=" + Q2rsp,
            "evtfile=" + quad_cleanevt,
            "badpixfile=" + quad_badpix_outfits,
            "clobber=YES",
            "history=YES",
        ]
    )
    sp.call(
        [
            "cztrspgen",
            "phafile=" + Q3pha,
            "rspfile=" + Q3rsp,
            "evtfile=" + quad_cleanevt,
            "badpixfile=" + quad_badpix_outfits,
            "clobber=YES",
            "history=YES",
        ]
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="CZTI Pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-d", help="Input directory", type=str, default="")
    parser.add_argument("-e", help="Input event file", type=str, default="")
    parser.add_argument("-m", help="Input mkf file", type=str, default="")
    parser.add_argument("-t", help="Input mkf threshold", type=str, default="")
    parser.add_argument("-l", help="Input livetime file", type=str, default="")
    parser.add_argument("-time", help="Input timebinsize", type=float, default=1)
    parser.add_argument("-emin", help="Energy lower limit", type=float, default=20)
    parser.add_argument("-emax", help="Energy upper limit", type=float, default=200)

    """-------------variable creation-----------------------------------------------------------------------------------"""
    args = parser.parse_args()
    if args.d == "":
        evt = args.e
        mkf = args.m
        mkf_threshold = args.t
        livetime = args.l
    else:
        evt = glob.glob(f"{args.d}/*bc.evt")[0]
        mkf = glob.glob(f"{args.d}/*.mkf")[0]
        mkf_threshold = glob.glob(f"{args.d}/*.txt")[0]
        livetime = glob.glob(f"{args.d}/*bc_livetime.fits")[0]
    emin = args.emin
    emax = args.emax
    timebinsize = args.time
    print(evt, mkf, mkf_threshold, livetime, timebinsize)

    """-------------calling functions-----------------------------------------------------------------------------------"""
    gtigen(evt, mkf, mkf_threshold)
    datasel(evt)
    pixclean(evt, livetime)
    evtclean(evt)
    flagbadpix(evt)
    bindata(evt, mkf, timebinsize, emin, emax)
