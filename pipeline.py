'''----------------------------------------------------------------------------------------------------------------'''
import subprocess as sp
from tkinter import *
from tkinter import filedialog 
import os
from astropy.io import fits
import matplotlib.pyplot as plt




'''----function for choosing the modueles---------------------------------------------------------------------------'''


def homegui():
    home_win = Tk()
    home_win.title('CZTI Pipeline')
    home_win.geometry("1280x720")
    home_win.configure(bg='#282a36')
    home_label = Label(home_win, text="Welcome to CZTI Pipeline", font='Courier 36 bold', bg='#282a36', fg='#f8f8f2')
    home_label.pack(pady=20)
    czti_pipeline_frame = Frame(home_win, bg='#282a36')
    czti_labels_frame = Frame(home_win, bg='#282a36')
    czti_pipeline_button = Button(czti_pipeline_frame, text='cztpipeline', font='Courier 20 bold', command=automaticgui, bg='#6272a4', fg='#f8f8f2', width=15)
   
    czti_modules_frame = Frame(home_win, bg='#282a36')
    
    cztgtigen_button = Button(czti_modules_frame, text="cztgtigen", font='Courier 20 bold', command=gtigui, bg='#6272a4', fg='#f8f8f2', width=15)
    cztgaas_button = Button(czti_modules_frame, text="cztgaas", font='Courier 20 bold', command=gaasgui, bg='#6272a4', fg='#f8f8f2', width=15)
    cztdatasel_button = Button(czti_modules_frame, text="cztdatasel", font='Courier 20 bold', command=dataselgui, bg='#6272a4', fg='#f8f8f2', width=15)
    cztpixclean_button = Button(czti_modules_frame, text="cztpixclean", font='Courier 20 bold', command=pixcleangui, bg='#6272a4', fg='#f8f8f2', width=15)
    cztevtclean_button = Button(czti_modules_frame, text="cztevtclean", font='Courier 20 bold', command=evtcleangui, bg='#6272a4', fg='#f8f8f2', width=15)
    cztflagbadpix_button = Button(czti_modules_frame, text="cztflagbadpix", font='Courier 20 bold', command=flagbadpixgui, bg='#6272a4', fg='#f8f8f2', width=15)
    cztbindata_button = Button(czti_modules_frame, text="cztbindata", font='Courier 20 bold', command=bindatagui, bg='#6272a4', fg='#f8f8f2', width=15)
    cztdpigen_button = Button(czti_modules_frame, text="cztdpigen", font='Courier 20 bold', command=dpigengui, bg='#6272a4', fg='#f8f8f2', width=15)
    cztimage_button = Button(czti_modules_frame, text="cztimage", font='Courier 20 bold', command=imagegui, bg='#6272a4', fg='#f8f8f2', width=15)
    cztrspgen_button = Button(czti_modules_frame, text="cztrspgen", font='Courier 20 bold', command=rspgen, bg='#6272a4', fg='#f8f8f2', width=15)

    cztgtigen_button.grid(row=1, column=0, padx=20, pady=20)
    cztgaas_button.grid(row=2, column=0,padx=20, pady=20)
    cztdatasel_button.grid(row=3, column=0,padx=20, pady=20)
    cztpixclean_button.grid(row=4, column=0,padx=20, pady=20)
    cztevtclean_button.grid(row=5, column=0,padx=20, pady=20)

    cztflagbadpix_button.grid(row=1, column=1,padx=0, pady=20)
    cztdpigen_button.grid(row=2, column=1,padx=0, pady=20)
    cztbindata_button.grid(row=3, column=1,padx=0, pady=20)
    cztimage_button.grid(row=4, column=1,padx=0, pady=20)
    cztrspgen_button.grid(row=5, column=1,padx=0, pady=20)

    czti_labels_frame.pack(side=TOP, padx=200)
    czti_pipeline_button.pack(side=TOP, padx=100, pady=250)
    czti_pipeline_frame.pack(side=LEFT, padx=30)

    czti_modules_frame.pack(side=RIGHT, padx=100, pady=75)
    home_win.mainloop()




'''----manual pipeline----------------------------------------------------------------------------------------------'''


'''-------------variable creation-----------------------------------------------------------------------------------'''

evt_label = ''
mkf_label = ''
mkf_thres_label = ''
evt = ''
mkf = ''
mkf_threshold = ''
gti = ''
out = ''
livetime = ''
badpix = ''
dpiordph = ''
imginp = ''

directory_label = ''
gtitype_label = ''
gtitype = 'QUAD'




'''----function definition for manual pipeline----------------------------------------------------------------------------------------------'''



'''--------------cztgtigen------------------------------------------------------------------------------------------'''

def man_gtigen(gti_bcevt, gti_mkf, gti_mkf_threshold, gti_usrgti='-'):
    sp.call(['cztgtigen', 
             'usergtifile='+gti_usrgti,
             'eventfile='+gti_bcevt,
             'mkffile='+gti_mkf,
             'thresholdfile='+gti_mkf_threshold, 
             'outfile='+out+'gti_out.gti',
             'clobber=YES', 
             'history=YES'])



'''--------------cztgaas--------------------------------------------------------------------------------------------'''

def man_gaas(gaas_bcevt, gaas_mkf):
    sp.call(['cztgaas', 
             'par_evtfile='+gaas_bcevt, 
             'par_mkffile='+gaas_mkf, 
             'par_outAspectfile='+out+'gaas_out.aspect', 
             'par_clobber=YES', 
             'par_history=YES'])



'''--------------cztdatasel-----------------------------------------------------------------------------------------'''

def man_datasel(datasel_bcevt, datasel_bcgti):
    sp.call(['cztdatasel', 
             'gtitype='+gtitype,
             'infile='+datasel_bcevt, 
             'gtifile='+datasel_bcgti, 
             'outfile='+out+'datasel_out.evt',
             'clobber=YES',
             'history=YES'])



'''--------------cztpixclean---------------------------------------------------------------------------------------'''

def man_pixclean(pixclean_bc_dsevt, pixclean_bc_livetimefits):
    sp.call(['cztpixclean', 
             'par_infile='+pixclean_bc_dsevt, 
             'par_inlivetimefile='+pixclean_bc_livetimefits, 
             'par_outfile1='+out+'pixclean_out.evt', 
             'par_outlivetimefile='+out+'pixclean_out_livetime.fits', 
             'par_badpixfile='+out+'pixclean_out_badpix.fits', 
             'par_det_tbinsize=1', 
             'par_pix_tbinsize=1', 
             'par_det_count_thresh=100', 
             'par_pix_count_thresh=1000'])



'''--------------cztevtclean---------------------------------------------------------------------------------------'''

def man_evtclean(evtclean_quad_pcevt):
    sp.call(['cztevtclean', 
             'infile='+evtclean_quad_pcevt, 
             'outfile='+out+'evtclean_out.evt', 
             'alphaval=0', 
             'vetorange=0-0',
             'clobber=YES',
             'history=YES'])



'''--------------cztflagbadpix-------------------------------------------------------------------------------------'''

def man_flagbadpix(flagbadpix_quad_badpixfits):
    sp.call(['cztflagbadpix', 
             'nbadpixFiles=1', 
             'badpixfile='+flagbadpix_quad_badpixfits, 
             'outfile='+out+'flagbadpix_out.fits', 
             'clobber=YES', 
             'history=YES'])



'''--------------cztbindata----------------------------------------------------------------------------------------'''

def man_bindata(bindata_quad_cleanevt, bindata_mkf, bindata_quad_badpix_outfits, bindata_quad_livetimefits):
    sp.call(['cztbindata', 
             'inevtfile='+bindata_quad_cleanevt, 
             'mkffile='+bindata_mkf, 
             'badpixfile='+bindata_quad_badpix_outfits, 
             'badpixthreshold=0', 
             'livetimefile='+bindata_quad_livetimefits,
             'outputtype=both',
             'energyrange=-',
             'timebinsize=1',
             'outfile='+out+'bindata',
             'outevtfile='+out+'bindata_weights.evt',
             'maskWeight=n',
             'rasrc=0',
             'decsrc=0',
             'clobber=YES',
             'history=YES'])



'''--------------cztdpigen-----------------------------------------------------------------------------------------'''

def man_dpigen(dpigen_quad_cleanevt, dpigen_quad_badpix_outfits):
    sp.call(['cztdpigen',
             'par_infile='+dpigen_quad_cleanevt,
             'par_badpixFile='+dpigen_quad_badpix_outfits,
             'par_badpixThreshold=0',
             'par_outDPHfile='+out+'dpigen_out.dph',
             'par_outDPIfile='+out+'dpigen_out.dpi',
             'par_quadsToProcess=0,1,2,3',
             'par_timerange=-',
             'par_ebins=-',
             'par_clobber=YES',
             'par_history=YES'])



'''--------------cztimage------------------------------------------------------------------------------------------'''

def man_image(image_dpiordph, image_imginp, image_aspectQ0, image_aspectQ1, image_aspectQ2, image_aspectQ3):
    sp.call(['cztimage', 
             'par_intype='+image_dpiordph,
             'par_infile='+image_imginp, 
             'par_outImgFile='+out+'image_out.img',
             'par_aspectfileQ0='+image_aspectQ0,
             'par_aspectfileQ1='+image_aspectQ1,
             'par_aspectfileQ2='+image_aspectQ2,
             'par_aspectfileQ3='+image_aspectQ3,
             'par_quadsToProcess=0,1,2,3', 
             'par_clobber=YES', 
             'par_history=YES'])



'''--------------cztrspgen-----------------------------------------------------------------------------------------'''

def man_rspgen():
    sp.call(['cztrspgen', 
             'phafile='+directory+Q0pha,
             'rspfile='+directory+Q0rsp,
             'evtfile='+directory+quad_cleanevt,
             'badpixfile='+directory+quad_badpix_outfits,
             'clobber=YES'
             'history=YES'])
    sp.call(['cztrspgen', 
             'phafile='+directory+Q1pha,
             'rspfile='+directory+Q1rsp,
             'evtfile='+directory+quad_cleanevt,
             'badpixfile='+directory+quad_badpix_outfits,
             'clobber=YES'
             'history=YES'])
    sp.call(['cztrspgen', 
             'phafile='+directory+Q2pha,
             'rspfile='+directory+Q2rsp,
             'evtfile='+directory+quad_cleanevt,
             'badpixfile='+directory+quad_badpix_outfits,
             'clobber=YES'
             'history=YES'])
    sp.call(['cztrspgen', 
             'phafile='+directory+Q3pha,
             'rspfile='+directory+Q3rsp,
             'evtfile='+directory+quad_cleanevt,
             'badpixfile='+directory+quad_badpix_outfits,
             'clobber=YES'
             'history=YES'])




'''----gui for manual pipeline-------------------------------------------------------------------------------------'''


'''--------------gtigen--------------------------------------------------------------------------------------------'''

def gtigui():
    global evt_label
    global gtigen_win
    global mkf_label
    global mkf_thres_label
    global directory_label
    gtigen_win = Tk() #creating the window
    gtigen_win.title('cztgtigen')
    gtigen_win.geometry("1280x720")
    gtigen_win.configure(bg='#282a36')
    input_frame = Frame(gtigen_win,bg='#282a36') #creating the frame
    ilabel = Label(input_frame, text="Choose the required input files", font='Courier 23 bold', bg='#282a36', fg='#f8f8f2')
    ilabel.pack(side=TOP, pady=20)
    b1 = Button(input_frame, text='Browse event file', command=open_evt_file, bg='#6272a4', fg='#f8f8f2',font='Courier 12 bold')
    b2 = Button(input_frame, text='Browse mkf file', command=open_mkf_file, bg='#6272a4', fg='#f8f8f2',font='Courier 12 bold')
    b3 = Button(input_frame, text='Browse mkf threshold file', command=open_mkf_thres_file, bg='#6272a4', fg='#f8f8f2',font='Courier 12 bold')
    evt_label = Label(input_frame, text="Event file: ", bg='#282a36', fg='#f8f8f2',font='Courier 14 bold')
    mkf_label = Label(input_frame, text="Mkf file: ", bg='#282a36',fg='#f8f8f2',font='Courier 14 bold')
    mkf_thres_label = Label(input_frame, text="Mkf Threshold file: ", bg='#282a36', fg='#f8f8f2',font='Courier 14')
    evt_label.pack(pady=10)
    mkf_label.pack(pady=10)
    mkf_thres_label.pack(pady=10)

    confirm_button = Button(gtigen_win, text='Confirm', command=confirm_gtigen, bg='#6272a4', fg='#f8f8f2',font='Courier 20 bold')
    button_frame = Frame(gtigen_win)
    
    b1.pack(side=LEFT, padx=10, pady=20)
    b2.pack(side=LEFT, padx=10, pady=20)
    b3.pack(side=LEFT, padx=10, pady=20)
    button_frame.pack(side=BOTTOM)

    output_frame = Frame(gtigen_win, bg='#282a36')
    olabel = Label(output_frame, text="Choose the required output directory", font='Courier 23 bold', bg='#282a36', fg='#f8f8f2')
    directory_label = Label(output_frame,  width=100, text="Output Directory: ", font='Courier 14', bg='#282a36', fg='#f8f8f2')
    directory_button = Button(output_frame, text="Open", command=output, bg='#6272a4', fg='#f8f8f2', font='Courier 12 bold')
    directory_button.pack(side=BOTTOM, pady=10)
    directory_label.pack(side=BOTTOM, padx=100, pady=45)
    olabel.pack(side=TOP, pady=10)

    confirm_button.pack(side=BOTTOM, pady=60)
    input_frame.pack(side=LEFT, padx=40)
    output_frame.pack(side=RIGHT, padx=0)

    gtigen_win.mainloop()



'''--------------gaas--------------------------------------------------------------------------------------------'''

def gaasgui():
    global evt_label
    global gaas_win
    global mkf_label
    global directory_label
    gaas_win = Tk() #creating the window
    gaas_win.title('cztgaas')
    gaas_win.geometry("1280x720")
    gaas_win.configure(bg='#282a36')
    input_frame = Frame(gaas_win,bg='#282a36') #creating the frame
    ilabel = Label(input_frame, text="Choose the required input files", font='Courier 23 bold', bg='#282a36', fg='#f8f8f2')
    ilabel.pack(side=TOP, pady=20)
    b1 = Button(input_frame, text='Browse event file', command=open_evt_file, bg='#6272a4', fg='#f8f8f2', font='Courier 12 bold')
    b2 = Button(input_frame, text='Browse mkf file', command=open_mkf_file, bg='#6272a4', fg='#f8f8f2', font='Courier 12 bold')
    evt_label = Label(input_frame, text="Event file: ", bg='#282a36', fg='#f8f8f2', font='Courier 14')
    mkf_label = Label(input_frame, text="Mkf file: ", bg='#282a36',fg='#f8f8f2', font='Courier 14')
    evt_label.pack(pady=10)
    mkf_label.pack(pady=10)

    confirm_button = Button(gaas_win, text='Confirm', command=confirm_gaas, bg='#6272a4', fg='#f8f8f2', font='Courier 20 bold')
    button_frame = Frame(gaas_win)
    button_frame.pack(side=BOTTOM)

    output_frame = Frame(gaas_win, bg='#282a36')
    olabel = Label(output_frame, text="Choose the required output directory", font='Courier 23 bold', bg='#282a36', fg='#f8f8f2')
    directory_label = Label(output_frame,  width=100, text="Output Directory: ", font='Courier 14', bg='#282a36', fg='#f8f8f2')
    directory_button = Button(output_frame, text="Open", command=output, bg='#6272a4', fg='#f8f8f2')
    directory_button.pack(side=BOTTOM, pady=10)
    directory_label.pack(side=BOTTOM, padx=100, pady=45)
    olabel.pack(side=TOP, pady=10)

    b1.pack(side=LEFT, padx=40, pady=20)
    b2.pack(side=LEFT, padx=40, pady=20)
    confirm_button.pack(side=BOTTOM, pady=60)
    input_frame.pack(side=LEFT, padx=80)
    output_frame.pack(side=RIGHT, padx=0)
    
    gaas_win.mainloop()



'''--------------datasel------------------------------------------------------------------------------------------'''

def dataselgui():
    global evt_label
    global datasel_win
    global gti_label
    global directory_label
    datasel_win = Tk() #creating the window
    datasel_win.title('cztdatasel')
    datasel_win.geometry("1280x720")
    datasel_win.configure(bg='#282a36')
    input_frame = Frame(datasel_win,bg='#282a36') #creating the frame
    ilabel = Label(input_frame, text="Choose the required input files", font='Courier 23 bold', bg='#282a36', fg='#f8f8f2')
    ilabel.pack(side=TOP, pady=20)
    b1 = Button(input_frame, text='Browse event file', command=open_evt_file, bg='#6272a4', fg='#f8f8f2', font='Courier 12 bold')
    b2 = Button(input_frame, text='Browse gti file', command=open_gti_file, bg='#6272a4', fg='#f8f8f2', font='Courier 12 bold')
    evt_label = Label(input_frame, text="Event file: ", bg='#282a36', fg='#f8f8f2',font='Courier 14')
    gti_label = Label(input_frame, text="GTI file: ", bg='#282a36',fg='#f8f8f2',font='Courier 14')
    evt_label.pack(pady=10)
    gti_label.pack(pady=10)

    confirm_button = Button(datasel_win, text='Confirm', command=confirm_datasel, bg='#6272a4', fg='#f8f8f2', font='Courier 12 bold')
    button_frame = Frame(datasel_win)
    button_frame.pack(side=BOTTOM)

    output_frame = Frame(datasel_win, bg='#282a36')
    olabel = Label(output_frame, text="Choose the required output directory", font='Courier 23 bold', bg='#282a36', fg='#f8f8f2')
    directory_label = Label(output_frame,  width=100, text="Output Directory: ",font='Courier 14', bg='#282a36', fg='#f8f8f2')
    directory_button = Button(output_frame, text="Open", command=output, bg='#6272a4', fg='#f8f8f2', font='Courier 12 bold')
    directory_button.pack(side=BOTTOM, pady=10)
    directory_label.pack(side=BOTTOM, padx=100, pady=45)
    olabel.pack(side=TOP, pady=10)

    b1.pack(side=LEFT, padx=40, pady=20)
    b2.pack(side=LEFT, padx=40, pady=20)
    confirm_button.pack(side=BOTTOM, pady=40)
    input_frame.pack(side=LEFT, padx=80)
    output_frame.pack(side=RIGHT, padx=0)

    datasel_win.mainloop()



'''--------------pixclean-----------------------------------------------------------------------------------------'''

def pixcleangui():
    global evt_label
    global pixclean_win
    global livetime_label
    global directory_label
    pixclean_win = Tk() #creating the window
    pixclean_win.title('cztpixclean')
    pixclean_win.geometry("1280x720")
    pixclean_win.configure(bg='#282a36')
    input_frame = Frame(pixclean_win,bg='#282a36') #creating the frame
    ilabel = Label(input_frame, text="Choose the required input files", font='Courier 23 bold', bg='#282a36', fg='#f8f8f2')
    ilabel.pack(side=TOP, pady=20)
    b1 = Button(input_frame, text='Browse event file', command=open_evt_file, bg='#6272a4', fg='#f8f8f2', font='Courier 12 bold')
    b2 = Button(input_frame, text='Browse livetime file', command=open_livetime_file, bg='#6272a4', fg='#f8f8f2', font='Courier 12 bold')
    evt_label = Label(input_frame, text="Event file: ", bg='#282a36', fg='#f8f8f2',font='Courier 14')
    livetime_label = Label(input_frame, text="Livetime file: ", bg='#282a36',fg='#f8f8f2',font='Courier 14')
    evt_label.pack(pady=10)
    livetime_label.pack(pady=10)

    confirm_button = Button(pixclean_win, text='Confirm', command=confirm_pixclean, bg='#6272a4', fg='#f8f8f2', font='Courier 20 bold')
    button_frame = Frame(pixclean_win)
    button_frame.pack(side=BOTTOM)

    output_frame = Frame(pixclean_win, bg='#282a36')
    olabel = Label(output_frame, text="Choose the required output directory", font='Courier 23 bold', bg='#282a36', fg='#f8f8f2')
    directory_label = Label(output_frame,  width=100, text="Output Directory: ",font='Courier 14', bg='#282a36', fg='#f8f8f2')
    directory_button = Button(output_frame, text="Open", command=output, bg='#6272a4', fg='#f8f8f2', font='Courier 12 bold')
    directory_button.pack(side=BOTTOM, pady=10)
    directory_label.pack(side=BOTTOM, padx=100, pady=45)
    olabel.pack(side=TOP, pady=10)
    
    b1.pack(side=LEFT, padx=40, pady=20)
    b2.pack(side=LEFT, padx=40, pady=20)
    confirm_button.pack(side=BOTTOM, pady=20)
    input_frame.pack(side=LEFT, padx=80)
    output_frame.pack(side=RIGHT, padx=0)

    pixclean_win.mainloop()



'''--------------evtclean-----------------------------------------------------------------------------------------'''

def evtcleangui():
    global evt_label
    global evtclean_win
    global directory_label
    evtclean_win = Tk() #creating the window
    evtclean_win.title('cztevtclean')
    evtclean_win.geometry("1280x720")
    evtclean_win.configure(bg='#282a36')
    input_frame = Frame(evtclean_win,bg='#282a36') #creating the frame
    ilabel = Label(input_frame, text="Choose the required input files", font='Courier 23 bold', bg='#282a36', fg='#f8f8f2')
    ilabel.pack(side=TOP, pady=20)
    b1 = Button(input_frame, text='Browse event file', command=open_evt_file, bg='#6272a4', fg='#f8f8f2', font='Courier 12 bold')
    evt_label = Label(input_frame, text="Event file: ", bg='#282a36', fg='#f8f8f2',font='Courier 14')
    evt_label.pack(pady=10)

    confirm_button = Button(evtclean_win, text='Confirm', command=confirm_evtclean, bg='#6272a4', fg='#f8f8f2', font='Courier 20 bold')
    button_frame = Frame(evtclean_win)
    button_frame.pack(side=BOTTOM)

    output_frame = Frame(evtclean_win, bg='#282a36')
    olabel = Label(output_frame, text="Choose the required output directory", font='Courier 23 bold', bg='#282a36', fg='#f8f8f2')
    directory_label = Label(output_frame,  width=100, text="Output Directory: ",font='Courier 14', bg='#282a36', fg='#f8f8f2')
    directory_button = Button(output_frame, text="Open", command=output, bg='#6272a4', fg='#f8f8f2', font='Courier 12 bold')
    directory_button.pack(side=BOTTOM, pady=10)
    directory_label.pack(side=BOTTOM, padx=100, pady=30)
    olabel.pack(side=TOP, pady=10)

    b1.pack(side=LEFT, padx=140, pady=20)
    confirm_button.pack(side=BOTTOM, pady=60)
    input_frame.pack(side=LEFT, padx=80)
    output_frame.pack(side=RIGHT, padx=0)

    evtclean_win.mainloop()



'''--------------flagbadpix---------------------------------------------------------------------------------------'''

def flagbadpixgui():
    global evt_label
    global flagbadpix_win
    global badpix_label
    global directory_label
    flagbadpix_win = Tk() #creating the window
    flagbadpix_win.title('cztflagbadpix')
    flagbadpix_win.geometry("1280x720")
    flagbadpix_win.configure(bg='#282a36')
    input_frame = Frame(flagbadpix_win,bg='#282a36') #creating the frame
    ilabel = Label(input_frame, text="Choose the required input files", font='Courier 23 bold', bg='#282a36', fg='#f8f8f2')
    ilabel.pack(side=TOP, pady=20)
    b1 = Button(input_frame, text='Browse event file', command=open_evt_file, bg='#6272a4', fg='#f8f8f2', font='Courier 12 bold')
    b2 = Button(input_frame, text='Browse badpix file', command=open_badpix_file, bg='#6272a4', fg='#f8f8f2', font='Courier 12 bold')
    evt_label = Label(input_frame, text="Event file: ", bg='#282a36', fg='#f8f8f2',font='Courier 14')
    badpix_label = Label(input_frame, text="Badpix file: ", bg='#282a36',fg='#f8f8f2',font='Courier 14')
    evt_label.pack(pady=10)
    badpix_label.pack(pady=10)

    confirm_button = Button(flagbadpix_win, text='Confirm', command=confirm_flagbadpix, bg='#6272a4', fg='#f8f8f2', font='Courier 20 bold')
    button_frame = Frame(flagbadpix_win)
    button_frame.pack(side=BOTTOM)
    
    output_frame = Frame(flagbadpix_win, bg='#282a36')
    olabel = Label(output_frame, text="Choose the required output directory", font='Courier 23 bold', bg='#282a36', fg='#f8f8f2')
    directory_label = Label(output_frame,  width=100, text="Output Directory: ", font='Courier 14', bg='#282a36', fg='#f8f8f2')
    directory_button = Button(output_frame, text="Open", command=output, bg='#6272a4', fg='#f8f8f2', font='Courier 12 bold')
    directory_button.pack(side=BOTTOM, pady=10)
    directory_label.pack(side=BOTTOM, padx=100, pady=45)
    olabel.pack(side=TOP, pady=10)

    b1.pack(side=LEFT, padx=40, pady=20)
    b2.pack(side=LEFT, padx=40, pady=20)
    confirm_button.pack(side=BOTTOM, pady=60)
    input_frame.pack(side=LEFT, padx=80)
    output_frame.pack(side=RIGHT, padx=0)
    
    flagbadpix_win.mainloop()



'''--------------bindata-------------------------------------------------------------------------------------------'''

def bindatagui():
    global evt_label
    global bindata_win
    global mkf_label
    global badpix_label
    global livetime_label
    global directory_label
    bindata_win = Tk() #creating the window
    bindata_win.title('cztbindata')
    bindata_win.geometry("1400x720")
    bindata_win.configure(bg='#282a36')
    input_frame = Frame(bindata_win,bg='#282a36') #creating the frame
    ilabel = Label(input_frame, text="Choose the required input files", font='Courier 23 bold', bg='#282a36', fg='#f8f8f2')
    ilabel.pack(side=TOP, pady=20)
    b1 = Button(input_frame, text='Browse event file', command=open_evt_file, bg='#6272a4', fg='#f8f8f2', font='Courier 12 bold')
    b2 = Button(input_frame, text='Browse mkf file', command=open_mkf_file, bg='#6272a4', fg='#f8f8f2', font='Courier 12 bold')
    b3 = Button(input_frame, text='Browse badpix file', command=open_badpix_file, bg='#6272a4', fg='#f8f8f2', font='Courier 12 bold')
    b4 = Button(input_frame, text='Browse livetime file', command=open_livetime_file, bg='#6272a4', fg='#f8f8f2', font='Courier 12 bold')
    evt_label = Label(input_frame, text="Event file: ", bg='#282a36', fg='#f8f8f2',font='Courier 14')
    mkf_label = Label(input_frame, text="Mkf file: ", bg='#282a36',fg='#f8f8f2',font='Courier 14')
    badpix_label = Label(input_frame, text="Badpix file: ", bg='#282a36',fg='#f8f8f2',font='Courier 14')
    
    livetime_label = Label(input_frame, text="Livetime file: ", bg='#282a36',fg='#f8f8f2', font='Courier 14')
    evt_label.pack(pady=10)
    mkf_label.pack(pady=10)
    badpix_label.pack(pady=10)
    livetime_label.pack(pady=10)

    confirm_button = Button(bindata_win, text='Confirm', command=confirm_bindata, bg='#6272a4', fg='#f8f8f2', font='Courier 20 bold')
    button_frame = Frame(bindata_win)
    button_frame.pack(side=BOTTOM)
    
    output_frame = Frame(bindata_win, bg='#282a36')
    olabel = Label(output_frame, text="Choose the required output directory", font='Courier 23 bold', bg='#282a36', fg='#f8f8f2')
    directory_label = Label(output_frame,  width=100, text="Output Directory: ",font='Courier 14', bg='#282a36', fg='#f8f8f2')
    directory_button = Button(output_frame, text="Open", command=output, bg='#6272a4', fg='#f8f8f2', font='Courier 12 bold')
    directory_button.pack(side=BOTTOM, pady=10)
    directory_label.pack(side=BOTTOM, padx=100, pady=45)
    olabel.pack(side=TOP, pady=10)
    
    b1.pack(side=LEFT, padx=10, pady=20)
    b2.pack(side=LEFT, padx=10, pady=20)
    b3.pack(side=LEFT, padx=10, pady=20)
    b4.pack(side=LEFT, padx=10, pady=20)
    confirm_button.pack(side=BOTTOM, pady=60)
    input_frame.pack(side=LEFT, padx=20)
    output_frame.pack(side=RIGHT, padx=0)

    bindata_win.mainloop()



'''--------------dpigen-------------------------------------------------------------------------------------------'''

def dpigengui():
    global evt_label
    global dpigen_win
    global badpix_label
    global directory_label
    dpigen_win = Tk() #creating the window
    dpigen_win.title('cztdpigen')
    dpigen_win.geometry("1280x720")
    dpigen_win.configure(bg='#282a36')
    input_frame = Frame(dpigen_win,bg='#282a36') #creating the frame
    ilabel = Label(input_frame, text="Choose the required input files", font='Courier 23 bold', bg='#282a36', fg='#f8f8f2')
    ilabel.pack(side=TOP, pady=20)
    b1 = Button(input_frame, text='Browse event file', command=open_evt_file, bg='#6272a4', fg='#f8f8f2', font='Courier 12 bold')
    b2 = Button(input_frame, text='Browse badpix file', command=open_badpix_file, bg='#6272a4', fg='#f8f8f2', font='Courier 12 bold')
    evt_label = Label(input_frame, text="Event file: ", bg='#282a36', fg='#f8f8f2',font='Courier 14')
    badpix_label = Label(input_frame, text="Badpix file: ", bg='#282a36',fg='#f8f8f2',font='Courier 14')
    evt_label.pack(pady=10)
    badpix_label.pack(pady=10)

    confirm_button = Button(dpigen_win, text='Confirm', command=confirm_dpigen, bg='#6272a4', fg='#f8f8f2', font='Courier 20 bold')
    button_frame = Frame(dpigen_win)
    button_frame.pack(side=BOTTOM)

    output_frame = Frame(dpigen_win, bg='#282a36')
    olabel = Label(output_frame, text="Choose the required output directory", font='Courier 23 bold', bg='#282a36', fg='#f8f8f2')
    directory_label = Label(output_frame,  width=100, text="Output Directory: ",font='Courier 14', bg='#282a36', fg='#f8f8f2')
    directory_button = Button(output_frame, text="Open", command=output, bg='#6272a4', fg='#f8f8f2', font='Courier 12 bold')
    directory_button.pack(side=BOTTOM, pady=10)
    directory_label.pack(side=BOTTOM, padx=100, pady=45)
    olabel.pack(side=TOP, pady=10)

    b1.pack(side=LEFT, padx=40, pady=20)
    b2.pack(side=LEFT, padx=10, pady=20)
    confirm_button.pack(side=BOTTOM, pady=40)
    input_frame.pack(side=LEFT, padx=80)
    output_frame.pack(side=RIGHT, padx=0)

    dpigen_win.mainloop()



'''--------------image--------------------------------------------------------------------------------------------'''

def imagegui():
    global dpiordph_label
    global image_win
    global imginp_label
    global aspect_label
    global dpiordph
    global imginp
    global aspectQ0
    global aspectQ1
    global aspectQ2
    global aspectQ3
    global directory_label
    global out
    image_win = Tk() #creating the window
    image_win.title('cztimage')
    image_win.geometry("1280x720")
    image_win.configure(bg='#282a36')
    input_frame = Frame(image_win,bg='#282a36') #creating the frame
    ilabel = Label(input_frame, text="Choose the required input files", font='Courier 23 bold', bg='#282a36', fg='#f8f8f2')
    ilabel.pack(side=TOP, pady=20)
    b1 = Button(input_frame, text='Browse DPI/DPH file', command=open_dpi_or_dph_file, bg='#6272a4', fg='#f8f8f2', font='Courier 12 bold')
    dpiordph_label = Label(input_frame, text="DPI/DPH file: ", bg='#282a36', fg='#f8f8f2',font='Courier 14')
    dpiordph_label.pack(pady=10)
    b2 = Button(input_frame, text="Choose aspect files directory", command=aspect_dir, bg='#6272a4', fg='#f8f8f2', font='Courier 12 bold')
    aspect_label = Label(input_frame, text='Aspect file directory: ', bg='#282a36', fg='#f8f8f2',font='Courier 14')
    aspectQ0 = out+evt.replace('_bc.evt','_aspect')+'_Q0'
    aspectQ1 = out+evt.replace('_bc.evt','_aspect')+'_Q1'
    aspectQ2 = out+evt.replace('_bc.evt','_aspect')+'_Q2'
    aspectQ3 = out+evt.replace('_bc.evt','_aspect')+'_Q3'
    confirm_button = Button(image_win, text='Confirm', command=confirm_image, bg='#6272a4', fg='#f8f8f2', font='Courier 20 bold')
    button_frame = Frame(image_win)
    button_frame.pack(side=BOTTOM)
        
    output_frame = Frame(image_win, bg='#282a36')
    olabel = Label(output_frame, text="Choose the required output directory", font='Courier 23 bold', bg='#282a36', fg='#f8f8f2')
    directory_label = Label(output_frame,  width=100, text="Output Directory: ",font='Courier 14', bg='#282a36', fg='#f8f8f2')
    directory_button = Button(output_frame, text="Open", command=output, bg='#6272a4', fg='#f8f8f2', font='Courier 12 bold')
    directory_button.pack(side=BOTTOM, pady=10)
    directory_label.pack(side=BOTTOM, padx=100, pady=45)
    olabel.pack(side=TOP, pady=10)

    b1.pack(side=LEFT, padx=120, pady=20)
    b2.pack(side=LEFT, padx=40, pady=20)
    aspect_label.pack(side=BOTTOM)
    confirm_button.pack(side=BOTTOM, pady=60)
    input_frame.pack(side=LEFT, padx=20)
    output_frame.pack(side=RIGHT, padx=0)

    image_win.mainloop()



'''--------------rspgen-------------------------------------------------------------------------------------------'''

def rspgengui():
    global rspgen_win
    global evt_label
    global mkf_label
    global mkf_thres_label
    global gti_label
    global livetime_label
    global badpix_label
    global directory_label
    rspgen_win = Tk() #creating the window
    rspgen_win.title('cztrspgen')
    rspgen_win.geometry("1280x720")
    rspgen_win.configure(bg='#282a36')
    input_frame = Frame(rspgen_win,bg='#282a36') #creating the frame
    ilabel = Label(input_frame, text="Choose the required input files", font='Courier 23 bold', bg='#282a36', fg='#f8f8f2')
    ilabel.pack(side=TOP, pady=20)
    
    b1 = Button(input_frame, text='Browse event file', command=open_evt_file, bg='#6272a4', fg='#f8f8f2', font='Courier 12 bold')
    b2 = Button(input_frame, text='Browse mkf file', command=open_mkf_file, bg='#6272a4', fg='#f8f8f2', font='Courier 12 bold')
    b3 = Button(input_frame, text='Browse mkf threshold file', command=open_mkf_thres_file, bg='#6272a4', fg='#f8f8f2', font='Courier 12 bold')
    b4 = Button(input_frame, text='Browse gti file', command=open_gti_file, bg='#6272a4', fg='#f8f8f2', font='Courier 12 bold')
    b5 = Button(input_frame, text='Browse livetime file', command=open_livetime_file, bg='#6272a4', fg='#f8f8f2', font='Courier 12 bold')
    b6 = Button(input_frame, text='Browse badpix file', command=open_badpix_file, bg='#6272a4', fg='#f8f8f2', font='Courier 12 bold')

    evt_label = Label(input_frame, text="Event file: ", bg='#282a36', fg='#f8f8f2',font='Courier 14')
    mkf_label = Label(input_frame, text="Mkf file: ", bg='#282a36',fg='#f8f8f2',font='Courier 14')
    mkf_thres_label = Label(input_frame, text="Mkf Threshold file: ", bg='#282a36',fg='#f8f8f2',font='Courier 14')
    gti_label = Label(input_frame, text="GTI file: ", bg='#282a36',fg='#f8f8f2',font='Courier 14')
    livetime_label = Label(input_frame, text="Livetime file: ", bg='#282a36',fg='#f8f8f2',font='Courier 14')
    badpix_label = Label(input_frame, text="Badpix file: ", bg='#282a36',fg='#f8f8f2',font='Courier 14')

    evt_label.pack(pady=10)
    mkf_label.pack(pady=10)
    mkf_thres_label.pack(pady=10)
    gti_label.pack(pady=10)
    livetime_label.pack(pady=10)
    badpix_label.pack(pady=10)
    
    confirm_button = Button(rspgen_win, text='Confirm', command=confirm_rspgen, bg='#6272a4', fg='#f8f8f2', font='Courier 20 bold')
    button_frame = Frame(rspgen_win)
    button_frame.pack(side=BOTTOM)

    output_frame = Frame(rspgen_win, bg='#282a36')
    olabel = Label(output_frame, text="Choose the required output directory", font='Courier 23 bold', bg='#282a36', fg='#f8f8f2')
    directory_label = Label(output_frame,  width=100, text="Output Directory: ",font='Courier 14', bg='#282a36', fg='#f8f8f2')
    directory_button = Button(output_frame, text="Open", command=output, bg='#6272a4', fg='#f8f8f2', font='Courier 12 bold')
    directory_button.pack(side=BOTTOM, pady=10)
    directory_label.pack(side=BOTTOM, padx=100, pady=45)
    olabel.pack(side=TOP, pady=10)

    b1.pack(side=LEFT, padx=10, pady=20)
    b2.pack(side=LEFT, padx=10, pady=20)
    b3.pack(side=LEFT, padx=10, pady=20)
    b4.pack(side=LEFT, padx=10, pady=20)
    b5.pack(side=LEFT, padx=10, pady=20)
    b6.pack(side=LEFT, padx=10, pady=20)
    confirm_button.pack(side=BOTTOM, pady=60)
    input_frame.pack(side=LEFT, padx=20)
    output_frame.pack(side=RIGHT, padx=0)

    rspgen_win.mainloop()






'''------functions used by gui-----------------------------------------------------------------------------------------'''


'''--------------opening and storing location of evt files----------------------------------------------------------------------------------------'''''

def open_evt_file():
    file = filedialog.askopenfile(mode='r', title="Open File")
    global evt
    evt =  (os.path.abspath(file.name))
    evt_label.config(text="Event file: "+ os.path.basename(evt), font='Courier 12 bold')



'''--------------opening and storing location of mkf files----------------------------------------------------------------------------------------'''

def open_mkf_file():
    file = filedialog.askopenfile(mode='r', title="Open File")
    global mkf
    mkf =  (os.path.abspath(file.name))
    mkf_label.config(text='Mkf file: '+ os.path.basename(mkf), font='Courier 12 bold')



'''--------------opening and storing location of mkf threshold files----------------------------------------------------------------------------------------'''

def open_mkf_thres_file():
    file = filedialog.askopenfile(mode='r', title="Open File")
    global mkf_threshold
    mkf_threshold =  (os.path.abspath(file.name))
    mkf_thres_label.config(text='Mkf Threshold file: '+ os.path.basename(mkf_threshold), font='Courier 12 bold')



'''--------------opening and storing location of gti files----------------------------------------------------------------------------------------'''

def open_gti_file():
    file = filedialog.askopenfile(mode='r', title="Open File")
    global gti
    gti =  (os.path.abspath(file.name))
    gti_label.config(text='GTI file: '+ os.path.basename(gti), font='Courier 12 bold')



'''--------------asking for the gtitype------------------------------------------------------------------------------------------------------------'''

def gti_type_common():
    global gtitype
    global gtitype_label
    gtitype = 'COMMON'
    gtitype_label.config(text='GTI Type: '+gtitype, font='Courier 12 bold')



'''--------------opening and storing location of livetime files----------------------------------------------------------------------------------------'''

def open_livetime_file():
    file = filedialog.askopenfile(mode='r', title="Open File")
    global livetime
    livetime =  (os.path.abspath(file.name))
    livetime_label.config(text='Livetime file: '+ os.path.basename(livetime), font='Courier 12 bold')



'''--------------opening and storing location of badpix files----------------------------------------------------------------------------------------'''

def open_badpix_file():
    file = filedialog.askopenfile(mode='r', title="Open File")
    global badpix
    badpix =  (os.path.abspath(file.name))
    badpix_label.config(text='Badpix file: '+ os.path.basename(badpix), font='Courier 12 bold')



'''--------------opening and storing location of dpi or dph files----------------------------------------------------------------------------------------'''

def open_dpi_or_dph_file():
    file = filedialog.askopenfile(mode='r', title="Open File")
    global imginp
    global dpiordph
    imginp =  (os.path.abspath(file.name))
    dpiordph_label.config(text='DPI/DPH file: '+ os.path.basename(imginp), font='Courier 12 bold')
    dpiordph = imginp.split('/')[-1].split('.')[-1]



'''--------------opening and storing location of aspect files----------------------------------------------------------------------------------------'''

def aspect_dir():
    global out
    global aspect_label
    out = filedialog.askdirectory() +'/'
    aspect_label.config(text="Aspect Files Directory: "+out, font='Courier 12 bold')



'''--------------opening and storing location of output files----------------------------------------------------------------------------------------'''

def output():
    global out
    global dir
    out = filedialog.askdirectory() +'/'
    directory_label.config(text="Output Directory: "+out, font='Courier 12 bold')



'''--------------functions for starting the modules---------------------------------------------------------------------------------------'''

def confirm_gtigen():
    global gtigen_win
    gtigen_win.destroy()
    man_gtigen(evt, mkf, mkf_threshold)

def confirm_gaas():
    global gaas_win
    gaas_win.destroy()
    man_gaas(evt, mkf)

def confirm_datasel():
    global datasel_win
    datasel_win.destroy()
    man_datasel(evt, gti)

def confirm_pixclean():
    global pixclean_win
    pixclean_win.destroy()
    man_pixclean(evt, livetime)

def confirm_evtclean():
    global evtclean_win
    evtclean_win.destroy()
    man_evtclean(evt)

def confirm_flagbadpix():
    global flagbadpix_win
    flagbadpix_win.destroy()
    man_flagbadpix(badpix)

def confirm_bindata():
    global bindata_win
    bindata_win.destroy()
    man_bindata(evt, mkf, badpix, livetime)

def confirm_dpigen():
    global dpigen_win
    dpigen_win.destroy()
    man_dpigen(evt, badpix)

def confirm_image():
    global image_win
    image_win.destroy()
    man_image(dpiordph, imginp, aspectQ0, aspectQ1, aspectQ2, aspectQ3)

def confirm_rspgen():
    global rspgen_win
    rspgen_win.destroy()
    man_rspgen()






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
             'timebinsize=1',
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



'''-----Functions to generate lightcurve, spectral plots and images--------------------------------------------------------------'''

def lightcurves():
    global evt
    lc1 = fits.open(evt.replace('bc.evt','quad_clean_Q0.lc'))
    lc2 = fits.open(evt.replace('bc.evt','quad_clean_Q1.lc'))
    lc3 = fits.open(evt.replace('bc.evt','quad_clean_Q2.lc'))
    lc4 = fits.open(evt.replace('bc.evt','quad_clean_Q3.lc'))
    lc1_data = lc1[1].data
    lc2_data = lc2[1].data
    lc3_data = lc3[1].data
    lc4_data = lc4[1].data
    clean_lc = [lc1_data, lc2_data, lc3_data, lc4_data]
    for i in range (0,4):
        plt.figure(dpi=250, figsize=(15,5))
        plt.title('Q'+str(i)+ ' lightcurve')
        plt.xlabel('Time (sec)')
        plt.ylabel('Counts/s')
        plt.plot(clean_lc[i].field('TIME'), clean_lc[i].field('RATE'), color='mediumslateblue')
        plt.savefig(evt.replace("bc.evt", "quad_clean_")+'Q'+str(i)+'_lightcurve.png')
        plt.close()




'''-----GUI for automatic pipeline-----------------------------------------------------------------------------------------'''


def automaticgui():
    global auto_win
    global evt_label
    global mkf_label
    global mkf_thres_label
    global livetime_label
    global evt
    global mkf
    global mkf_threshold
    global livetime
    auto_win = Tk() #creating the window
    auto_win.title('Automatic Pipeline')
    auto_win.configure(bg='#282a36')
    auto_win.geometry("1280x720")
    input_frame = Frame(auto_win,bg='#282a36') #creating the frame
    ilabel = Label(input_frame, text="Choose the required input files", font='Courier 25 bold', bg='#282a36', fg='#f8f8f2')
    ilabel.pack(side=TOP, pady=20)
    b1 = Button(input_frame, text='Browse event file', command=open_evt_file,bg='#6272a4', fg='#f8f8f2', font='Courier 12 bold')
    b2 = Button(input_frame, text='Browse mkf file', command=open_mkf_file,bg='#6272a4', fg='#f8f8f2', font='Courier 12 bold')
    b3 = Button(input_frame, text='Browse mkf threshold file', command=open_mkf_thres_file,bg='#6272a4', fg='#f8f8f2', font='Courier 12 bold')
    b4 = Button(input_frame, text='Browse livetime file', command=open_livetime_file,bg='#6272a4', fg='#f8f8f2', font='Courier 12 bold')
    evt_label = Label(input_frame, text="Event file: "+evt, bg='#282a36',fg='#f8f8f2',font='Courier 14')
    mkf_label = Label(input_frame, text="Mkf file: "+mkf, bg='#282a36',fg='#f8f8f2',font='Courier 14')
    mkf_thres_label = Label(input_frame, text="Mkf Threshold file: "+mkf_threshold, bg='#282a36',fg='#f8f8f2',font='Courier 14')
    livetime_label = Label(input_frame, text="Livetime file: "+livetime, bg='#282a36',fg='#f8f8f2',font='Courier 14')
    evt_label.pack(pady=10)
    mkf_label.pack(pady=10)
    mkf_thres_label.pack(pady=10)
    livetime_label.pack(pady=10)

    confirm_button = Button(auto_win, text='Confirm', command=confirm_auto, bg='#6272a4', fg='#f8f8f2', font='Courier 20 bold')
    button_frame = Frame(auto_win)

    b1.pack(side=LEFT, padx=20, pady=10)
    b2.pack(side=LEFT, padx=20, pady=10)
    b3.pack(side=LEFT, padx=20, pady=10)
    b4.pack(side=LEFT, padx=20, pady=10)
    button_frame.pack(side=BOTTOM)
    confirm_button.pack(side=BOTTOM, pady=40)
    input_frame.pack(side=TOP, pady=120)
    
    auto_win.mainloop()



'''-----function to run the automatic pipeline-----------------------------------------------------------------------------------------'''

def confirm_auto():
    global auto_win
    auto_win.destroy()
    gtigen()
    gaas()
    datasel()
    pixclean()
    evtclean()
    flagbadpix()
    bindata()
    dpigen()
    image()
    rspgen()
    lightcurves()


'''-----calling the home gui-----------------------------------------------------------------------------------------'''

homegui()