'''Written by ravioli1369@gmail.com and shreyasbharadwaj04@gmail.com'''


'''----------------------------------------------------------------------------------------------------------------'''
import subprocess as sp
import os
from astropy.io import fits
import matplotlib.pyplot as plt
import argparse

'''----passing arguments----------------------------------------------------------------------------------------------'''
parser = argparse.ArgumentParser(description='CZTI Pipeline', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--e', '--evt', help='Input event file', type=str, default='')
parser.add_argument('--m', '--mkf', help='Input mkf file', type=str, default='')
parser.add_argument('--t', '--threshold', help='Input mkf threshold', type=str, default='')
parser.add_argument('--l', '--livetime', help='Input livetime file', type=str, default='')
parser.add_argument('--timebinsize', help='Input timebinsize', type=float, default=1)
'''-------------variable creation-----------------------------------------------------------------------------------'''
args = parser.parse_args()
evt = args.e
mkf = args.m
mkf_threshold = args.t
livetime = args.l
timebinsize = args.timebinsize

'''-----Automatic Pipeline-----------------------------------------------------------------------------------------'''


'''-----functions for automatic pipeline---------------------------------------------------------------------------'''

#cztgtigen
def gtigen():
    sp.call(['cztgtigen', 
             'eventfile='+evt, 
             'mkffile='+mkf, 
             'thresholdfile='+mkf_threshold, 
             'outfile='+evt.replace(".evt", ".gti"),
             'usergtifile=-',
             'clobber=YES',
             'history=YES'])

#cztgaas
def gaas():
    sp.call(['cztgaas', 
             'par_evtfile='+evt, 
             'par_mkffile='+mkf, 
             'par_outAspectfile='+evt.replace("_bc.evt", ".aspect"), 
             'par_clobber=YES', 
             'par_history=YES'])

#cztdatasel
def datasel():
    sp.call(['cztdatasel', 
             'infile='+evt,
             'gtitype=QUAD',
             'gtifile='+evt.replace(".evt", ".gti"), 
             'outfile='+evt.replace("bc.evt", "bc_ds.evt"),
             'clobber=YES',
             'history=YES'])

#cztpixclean
def pixclean():
    sp.call(['cztpixclean', 
             'par_infile='+evt.replace("bc.evt", "bc_ds.evt"), 
             'par_inlivetimefile='+livetime, 
             'par_outfile1='+evt.replace("bc.evt", "quad_pc.evt"), 
             'par_outlivetimefile='+evt.replace("bc.evt", "quad_livetime.fits"), 
             'par_badpixfile='+evt.replace("bc.evt", "quad_badpix.fits"), 
             'par_det_tbinsize=1', 
             'par_pix_tbinsize=1', 
             'par_det_count_thresh=100', 
             'par_pix_count_thresh=1000'])

#cztevtclean
def evtclean():
    sp.call(['cztevtclean', 
             'infile='+evt.replace("bc.evt", "quad_pc.evt"), 
             'outfile='+evt.replace("bc.evt", "quad_clean.evt"),
             'alphaval='+'0', 
             'vetorange='+'0-0',
             'clobber=YES',
             'history=YES'])
    
#cztflagbadpix
def flagbadpix():
    sp.call(['cztflagbadpix', 
             'nbadpixFiles=1', 
             'badpixfile='+evt.replace("bc.evt", "quad_badpix.fits"),
             'outfile='+evt.replace("bc.evt", "quad_badpix_out.fits"),
             'clobber=YES', 
             'history=YES'])

#cztbindata
def bindata():
    sp.call(['cztbindata', 
             'inevtfile='+evt.replace("bc.evt", "quad_clean.evt"), 
             'mkffile='+mkf, 
             'badpixfile='+evt.replace("bc.evt", "quad_badpix_out.fits"), 
             'badpixthreshold=0', 
             'livetimefile='+evt.replace("bc.evt", "quad_livetime.fits"),
             'outputtype=both',
             'energyrange=-',
             'timebinsize='+str(timebinsize),
             'outfile='+evt.replace("bc.evt", "quad_clean"),
             'outevtfile='+evt.replace("bc.evt", "quad_clean_weight.evt"),
             'maskWeight=n',
             'rasrc=0',
             'decsrc=0',
             'clobber=YES',
             'history=YES'])

#cztdpigen
def dpigen():
    sp.call(['cztdpigen',
             'par_infile='+evt.replace("bc.evt", "quad_clean.evt"),
             'par_badpixFile='+evt.replace("bc.evt", "quad_badpix_out.fits"),
             'par_badpixThreshold=0',
             'par_outDPHfile='+evt.replace("bc.evt", "quad_clean.dpi"),
             'par_outDPIfile='+evt.replace("bc.evt", "quad_clean.dph"),
             'par_quadsToProcess=0,1,2,3',
             'par_timerange=-',
             'par_ebins=-',
             'par_clobber=YES',
             'par_history=YES'])

#cztimage
def image():
    sp.call(['cztimage', 
             'par_intype=dpi',
             'par_infile='+evt.replace("bc.evt", "quad_clean.dpi"), 
             'par_aspectfileQ0='+evt.replace("_bc.evt", ".aspect")+"_Q0", 
             'par_aspectfileQ1='+evt.replace("_bc.evt", ".aspect")+"_Q1",
             'par_aspectfileQ2='+evt.replace("_bc.evt", ".aspect")+"_Q2",
             'par_aspectfileQ3='+evt.replace("_bc.evt", ".aspect")+"_Q3",
             'par_outImgFile='+evt.replace("bc.evt", "quad_clean.img"),
             'par_quadsToProcess=0,1,2,3', 
             'par_clobber=YES', 
             'par_history=YES'])

#cztrspgen
def rspgen():
    quad_clean = evt.replace("bc.evt", "quad_clean")
    quad_cleanevt = evt.replace("bc.evt", "quad_clean.evt")
    quad_badpix_outfits = evt.replace("bc.evt", "quad_badpix_out.fits")
    Q0pha = quad_clean+"_Q0.pha"
    Q1pha = quad_clean+"_Q1.pha"
    Q2pha = quad_clean+"_Q2.pha"
    Q3pha = quad_clean+"_Q3.pha"
    Q0rsp = quad_clean+"_Q0.rsp"
    Q1rsp = quad_clean+"_Q1.rsp"
    Q2rsp = quad_clean+"_Q2.rsp"
    Q3rsp = quad_clean+"_Q3.rsp"
    sp.call(['cztrspgen', 
             'phafile='+Q0pha,
             'rspfile='+Q0rsp,
             'evtfile='+quad_cleanevt,
             'badpixfile='+quad_badpix_outfits,
             'clobber=YES',
             'history=YES'])
    sp.call(['cztrspgen', 
             'phafile='+Q1pha,
             'rspfile='+Q1rsp,
             'evtfile='+quad_cleanevt,
             'badpixfile='+quad_badpix_outfits,
             'clobber=YES',
             'history=YES'])
    sp.call(['cztrspgen', 
             'phafile='+Q2pha,
             'rspfile='+Q2rsp,
             'evtfile='+quad_cleanevt,
             'badpixfile='+quad_badpix_outfits,
             'clobber=YES',
             'history=YES'])
    sp.call(['cztrspgen', 
             'phafile='+Q3pha,
             'rspfile='+Q3rsp,
             'evtfile='+quad_cleanevt,
             'badpixfile='+quad_badpix_outfits,
             'clobber=YES',
             'history=YES'])

if __name__ == '__main__':
    gtigen()
    gaas()
    datasel()
    pixclean()
    evtclean()
    flagbadpix()
    bindata()
    # dpigen()
    # image()
    # rspgen()