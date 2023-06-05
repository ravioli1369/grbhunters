'''----------------------------------------------------------------------------------------------------------------'''
import subprocess as sp
from tkinter import *
from tkinter import filedialog 
import os




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
             'vetorange=0-0'])
    
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
def man_image(image_dpiordph, image_imginp, image_aspectQ0, image_aspectQ1, image_aspectQ2, image_aspectQ3): #imginp=quad_clean+dpiordph  take dpiordph from user
    sp.call(['cztimage', 
             'par_intype='+image_dpiordph,
             'par_infile='+image_imginp, 
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
    global evt
    global mkf
    global mkf_threshold
    global dir
    gtigen_win = Tk() #creating the window
    gtigen_win.title('cztgtigen')
    gtigen_win.geometry("1200x400")
    input_frame = Frame(gtigen_win) #creating the frame
    ilabel = Label(input_frame, text="Choose the required input files", font='Arial 16 bold')
    ilabel.pack(side=TOP)
    b1 = Button(input_frame, text='Browse event file', command=open_evt_file)
    b2 = Button(input_frame, text='Browse mkf file', command=open_mkf_file)
    b3 = Button(input_frame, text='Browse mkf threshold file', command=open_mkf_thres_file)
    dir = StringVar()
    evt_label = Label(input_frame, text="Event file: "+evt)
    mkf_label = Label(input_frame, text="Mkf file: "+mkf)
    mkf_thres_label = Label(input_frame, text="Mkf Threshold file: "+mkf_threshold)
    evt_label.pack(pady=10)
    mkf_label.pack(pady=10)
    mkf_thres_label.pack(pady=10)

    confirm_button = Button(gtigen_win, text='Confirm', command=confirm_gtigen)
    button_frame = Frame(gtigen_win)
    directory_button = Button(button_frame, text="Open", command=output)
    directory_label = Label(gtigen_win,  width=200, bg='white', fg='black', textvariable=dir)

    directory_label.pack(side=TOP, padx=10, pady=20)
    directory_button.pack(side=TOP)
    button_frame.pack(side=BOTTOM)

    output_frame = Frame(gtigen_win)
    olabel = Label(output_frame, text="Choose the required output directory", font='Arial 16 bold')
    olabel.pack(side=TOP, pady=10)

    b1.pack(side=LEFT)
    b2.pack(side=LEFT)
    b3.pack(side=LEFT)
    confirm_button.pack(side=BOTTOM)
    input_frame.pack(side=LEFT)
    output_frame.pack(side=RIGHT, padx=40)

    gtigen_win.mainloop()

'''--------------gaas--------------------------------------------------------------------------------------------'''
def gaasgui():
    global gaas_win
    global evt_label
    global mkf_label
    global evt
    global mkf
    global dir
    gaas_win = Tk() #creating the window
    gaas_win.title('cztgaas')
    gaas_win.geometry("1200x400")
    input_frame = Frame(gaas_win) #creating the frame
    ilabel = Label(input_frame, text="Choose the required input files", font='Arial 16 bold')
    ilabel.pack(side=TOP)
    b1 = Button(input_frame, text='Browse event file', command=open_evt_file) #browse button for event file
    b2 = Button(input_frame, text='Browse mkf file', command=open_mkf_file) #browse button for mkf file
    dir = StringVar()
    evt_label = Label(input_frame, text="Event file: "+evt)
    mkf_label = Label(input_frame, text="Mkf file: "+mkf)
    evt_label.pack(pady=10)
    mkf_label.pack(pady=10)

    confirm_button = Button(gaas_win, text='Confirm', command=confirm_gaas)
    button_frame = Frame(gaas_win)
    directory_button = Button(button_frame, text="Open", command=output)
    directory_label = Label(gaas_win,  width=200, bg='white', fg='black', textvariable=dir)

    directory_label.pack(side=TOP, padx=10, pady=20)
    directory_button.pack(side=TOP)
    button_frame.pack(side=BOTTOM)

    output_frame = Frame(gaas_win)
    olabel = Label(output_frame, text="Choose the required output directory", font='Arial 16 bold')
    olabel.pack(side=TOP, pady=10)

    b1.pack(side=LEFT)
    b2.pack(side=LEFT)
    confirm_button.pack(side=BOTTOM)
    input_frame.pack(side=LEFT)
    output_frame.pack(side=RIGHT, padx=40)

    gaas_win.mainloop()

'''--------------datasel------------------------------------------------------------------------------------------'''
def dataselgui():
    global datasel_win
    global evt_label
    global gti_label
    global evt
    global gti
    global dir
    datasel_win = Tk() #creating the window
    datasel_win.title('cztdatasel')
    datasel_win.geometry("1200x400")
    input_frame = Frame(datasel_win) #creating the frame
    ilabel = Label(input_frame, text="Choose the required input files", font='Arial 16 bold')
    ilabel.pack(side=TOP)
    b1 = Button(input_frame, text='Browse event file', command=open_evt_file)
    b2 = Button(input_frame, text='Browse GTI file', command=open_gti_file)
    dir = StringVar()
    evt_label = Label(input_frame, text="Event file: "+evt)
    gti_label = Label(input_frame, text="GTI file: "+gti)
    evt_label.pack(pady=10)
    gti_label.pack(pady=10)

    confirm_button = Button(datasel_win, text='Confirm', command=confirm_datasel)
    button_frame = Frame(datasel_win)
    directory_button = Button(button_frame, text="Open", command=output)
    directory_label = Label(datasel_win,  width=200, bg='white', fg='black', textvariable=dir)

    directory_label.pack(side=TOP, padx=10, pady=20)
    directory_button.pack(side=TOP)
    button_frame.pack(side=BOTTOM)

    output_frame = Frame(datasel_win)
    olabel = Label(output_frame, text="Choose the required output directory", font='Arial 16 bold')
    olabel.pack(side=TOP, pady=10)

    b1.pack(side=LEFT)
    b2.pack(side=LEFT)
    confirm_button.pack(side=BOTTOM)
    input_frame.pack(side=LEFT)
    output_frame.pack(side=RIGHT, padx=40)

    datasel_win.mainloop()

'''--------------pixclean-----------------------------------------------------------------------------------------'''
def pixcleangui():
    global pixclean_win
    global evt_label
    global fits_label
    global evt
    global fits
    global dir
    pixclean_win = Tk() #creating the window
    pixclean_win.title('cztpixclean')
    pixclean_win.geometry("1200x400")
    input_frame = Frame(pixclean_win) #creating the frame
    ilabel = Label(input_frame, text="Choose the required input files", font='Arial 16 bold')
    ilabel.pack(side=TOP)
    b1 = Button(input_frame, text='Browse event file', command=open_evt_file)
    b2 = Button(input_frame, text='Browse livetime file', command=open_livetime_file)
    dir = StringVar()
    evt_label = Label(input_frame, text="Event file: "+evt)
    fits_label = Label(input_frame, text="Livetime file: "+fits)
    evt_label.pack(pady=10)
    fits_label.pack(pady=10)

    confirm_button = Button(pixclean_win, text='Confirm', command=confirm_pixclean)
    button_frame = Frame(pixclean_win)
    directory_button = Button(button_frame, text="Open", command=output)
    directory_label = Label(pixclean_win,  width=200, bg='white', fg='black', textvariable=dir)

    directory_label.pack(side=TOP, padx=10, pady=20)
    directory_button.pack(side=TOP)
    button_frame.pack(side=BOTTOM)

    output_frame = Frame(pixclean_win)
    olabel = Label(output_frame, text="Choose the required output directory", font='Arial 16 bold')
    olabel.pack(side=TOP, pady=10)

    b1.pack(side=LEFT)
    b2.pack(side=LEFT)
    confirm_button.pack(side=BOTTOM)
    input_frame.pack(side=LEFT)
    output_frame.pack(side=RIGHT, padx=40)

    pixclean_win.mainloop()

'''--------------evtclean-----------------------------------------------------------------------------------------'''
def evtcleangui():
    global evtclean_win
    global evt_label
    global evt
    global dir
    evtclean_win = Tk() #creating the window
    evtclean_win.title('cztevtclean')
    evtclean_win.geometry("1200x400")
    input_frame = Frame(evtclean_win) #creating the frame
    ilabel = Label(input_frame, text="Choose the required input files", font='Arial 16 bold')
    ilabel.pack(side=TOP)
    b1 = Button(input_frame, text='Browse event file', command=open_evt_file)
    dir = StringVar()
    evt_label = Label(input_frame, text="Event file: "+evt)
    evt_label.pack(pady=10)

    confirm_button = Button(evtclean_win, text='Confirm', command=confirm_evtclean)
    button_frame = Frame(evtclean_win)
    directory_button = Button(button_frame, text="Open", command=output)
    directory_label = Label(evtclean_win,  width=200, bg='white', fg='black', textvariable=dir)

    directory_label.pack(side=TOP, padx=10, pady=20)
    directory_button.pack(side=TOP)
    button_frame.pack(side=BOTTOM)

    output_frame = Frame(evtclean_win)
    olabel = Label(output_frame, text="Choose the required output directory", font='Arial 16 bold')
    olabel.pack(side=TOP, pady=10)

    b1.pack(side=LEFT)
    confirm_button.pack(side=BOTTOM)
    input_frame.pack(side=LEFT)
    output_frame.pack(side=RIGHT, padx=40)

    evtclean_win.mainloop()

'''--------------flagbadpix---------------------------------------------------------------------------------------'''
def flagbadpixgui():
    global flagbadpix_win
    global fits_label
    global fits
    global dir
    flagbadpix_win = Tk() #creating the window
    flagbadpix_win.title('cztflagbadpix')
    flagbadpix_win.geometry("1200x400")
    input_frame = Frame(flagbadpix_win) #creating the frame
    ilabel = Label(input_frame, text="Choose the required input files", font='Arial 16 bold')
    ilabel.pack(side=TOP)
    b1 = Button(input_frame, text='Browse badpix file', command=open_badpix_file)
    dir = StringVar()
    fits_label = Label(input_frame, text="Badpix file: "+fits)
    fits_label.pack(pady=10)

    confirm_button = Button(flagbadpix_win, text='Confirm', command=confirm_flagbadpix)
    button_frame = Frame(flagbadpix_win)
    directory_button = Button(button_frame, text="Open", command=output)
    directory_label = Label(flagbadpix_win,  width=200, bg='white', fg='black', textvariable=dir)

    directory_label.pack(side=TOP, padx=10, pady=20)
    directory_button.pack(side=TOP)
    button_frame.pack(side=BOTTOM)

    output_frame = Frame(flagbadpix_win)
    olabel = Label(output_frame, text="Choose the required output directory", font='Arial 16 bold')
    olabel.pack(side=TOP, pady=10)

    b1.pack(side=LEFT)
    confirm_button.pack(side=BOTTOM)
    input_frame.pack(side=LEFT)
    output_frame.pack(side=RIGHT, padx=40)

    flagbadpix_win.mainloop()

'''--------------bindata------------------------------------------------------------------------------------------'''
def bindatagui():
    global bindata_win
    global evt_label
    global fits_label
    global evt
    global livetime
    global badpix
    global dir
    bindata_win = Tk() #creating the window
    bindata_win.title('cztbindata')
    bindata_win.geometry("1200x400")
    input_frame = Frame(bindata_win) #creating the frame
    ilabel = Label(input_frame, text="Choose the required input files", font='Arial 16 bold')
    ilabel.pack(side=TOP)
    b1 = Button(input_frame, text='Browse event file', command=open_evt_file)
    b2 = Button(input_frame, text='Browse livetime file', command=open_livetime_file)
    b3 = Button(input_frame, text='Browse badpix file', command=open_badpix_file)
    dir = StringVar()
    evt_label = Label(input_frame, text="Event file: "+evt)
    fits_label = Label(input_frame, text="Livetime file: "+fits)
    evt_label.pack(pady=10)
    fits_label.pack(pady=10)

    confirm_button = Button(bindata_win, text='Confirm', command=confirm_bindata)
    button_frame = Frame(bindata_win)
    directory_button = Button(button_frame, text="Open", command=output)
    directory_label = Label(bindata_win,  width=200, bg='white', fg='black', textvariable=dir)

    directory_label.pack(side=TOP, padx=10, pady=20)
    directory_button.pack(side=TOP)
    button_frame.pack(side=BOTTOM)

    output_frame = Frame(bindata_win)
    olabel = Label(output_frame, text="Choose the required output directory", font='Arial 16 bold')
    olabel.pack(side=TOP, pady=10)

    b1.pack(side=LEFT)
    b2.pack(side=LEFT)
    b3.pack(side=LEFT)
    confirm_button.pack(side=BOTTOM)
    input_frame.pack(side=LEFT)
    output_frame.pack(side=RIGHT, padx=40)

    bindata_win.mainloop()

'''--------------dpigen-------------------------------------------------------------------------------------------'''
def dpigengui():
    global dpigen_win
    global evt_label
    global fits_label
    global evt
    global badpix
    global dir
    dpigen_win = Tk() #creating the window
    dpigen_win.title('cztdpigen')
    dpigen_win.geometry("1200x400")
    input_frame = Frame(dpigen_win) #creating the frame
    ilabel = Label(input_frame, text="Choose the required input files", font='Arial 16 bold')
    ilabel.pack(side=TOP)
    b1 = Button(input_frame, text='Browse event file', command=open_evt_file)
    b2 = Button(input_frame, text='Browse badpix file', command=open_badpix_file)
    dir = StringVar()
    evt_label = Label(input_frame, text="Event file: "+evt)
    fits_label = Label(input_frame, text="Badpix file: "+fits)
    evt_label.pack(pady=10)
    fits_label.pack(pady=10)

    confirm_button = Button(dpigen_win, text='Confirm', command=confirm_dpigen)
    button_frame = Frame(dpigen_win)
    directory_button = Button(button_frame, text="Open", command=output)
    directory_label = Label(dpigen_win,  width=200, bg='white', fg='black', textvariable=dir)

    directory_label.pack(side=TOP, padx=10, pady=20)
    directory_button.pack(side=TOP)
    button_frame.pack(side=BOTTOM)

    output_frame = Frame(dpigen_win)
    olabel = Label(output_frame, text="Choose the required output directory", font='Arial 16 bold')
    olabel.pack(side=TOP, pady=10)

    b1.pack(side=LEFT)
    b2.pack(side=LEFT)
    confirm_button.pack(side=BOTTOM)
    input_frame.pack(side=LEFT)
    output_frame.pack(side=RIGHT, padx=40)

    dpigen_win.mainloop()

def imagegui():
    global image_win
    global evt_label
    global fits_label
    global evt
    global badpix
    global dir
    image_win = Tk() #creating the window
    image_win.title('cztimage')
    image_win.geometry("1200x400")
    input_frame = Frame(image_win) #creating the frame
    ilabel = Label(input_frame, text="Choose the required input files", font='Arial 16 bold')
    ilabel.pack(side=TOP)
    b1 = Button(input_frame, text='DPI/DPH', command=choose_dpiordph)
    b2 = Button(input_frame, text='Browse event file', command=open_evt_file)
    b3 = Button(input_frame, text='Browse aspect file for Q0', command=open_aspectQ0_file)
    b4 = Button(input_frame, text='Browse aspect file for Q1', command=open_aspectQ1_file)
    b5 = Button(input_frame, text='Browse aspect file for Q2', command=open_aspectQ2_file)
    b6 = Button(input_frame, text='Browse aspect file for Q3', command=open_aspectQ3_file)
    dir = StringVar()
    evt_label = Label(input_frame, text="Event file: "+evt)
    fits_label = Label(input_frame, text="Badpix file: "+fits)
    evt_label.pack(pady=10)
    fits_label.pack(pady=10)

    confirm_button = Button(image_win, text='Confirm', command=confirm_image)
    button_frame = Frame(image_win)
    directory_button = Button(button_frame, text="Open", command=output)
    directory_label = Label(image_win,  width=200, bg='white', fg='black', textvariable=dir)

    directory_label.pack(side=TOP, padx=10, pady=20)
    directory_button.pack(side=TOP)
    button_frame.pack(side=BOTTOM)

    output_frame = Frame(image_win)
    olabel = Label(output_frame, text="Choose the required output directory", font='Arial 16 bold')
    olabel.pack(side=TOP, pady=10)

    b1.pack(side=LEFT)
    b2.pack(side=LEFT)
    b3.pack(side=LEFT)
    b4.pack(side=LEFT)
    b5.pack(side=LEFT)
    b6.pack(side=LEFT)
    confirm_button.pack(side=BOTTOM)
    input_frame

def rspgengui():
    global rspgen_win
    global evt_label
    global fits_label
    global evt
    global badpix
    global dir
    rspgen_win = Tk() #creating the window
    rspgen_win.title('cztrspgen')
    rspgen_win.geometry("1200x400")
    input_frame = Frame(rspgen_win) #creating the frame
    ilabel = Label(input_frame, text="Choose the required input files", font='Arial 16 bold')
    ilabel.pack(side=TOP)

    confirm_button = Button(rspgen_win, text='Confirm', command=confirm_rspgen)
    button_frame = Frame(rspgen_win)
    directory_button = Button(button_frame, text="Open", command=output)
    directory_label = Label(rspgen_win,  width=200, bg='white', fg='black', textvariable=dir)

    directory_label.pack(side=TOP, padx=10, pady=20)
    directory_button.pack(side=TOP)
    button_frame.pack(side=BOTTOM)

    output_frame = Frame(rspgen_win)
    olabel = Label(output_frame, text="Choose the required output directory", font='Arial 16 bold')
    olabel.pack(side=TOP, pady=10)


    confirm_button.pack(side=BOTTOM)
    input_frame.pack(side=LEFT)
    output_frame.pack(side=RIGHT, padx=40)

    rspgen_win.mainloop()


'''------functions for gui-----------------------------------------------------------------------------------------'''

'''--------------opening and storing location of evt files----------------------------------------------------------------------------------------'''''
def open_evt_file():
    file = filedialog.askopenfile(mode='r', title="Open File")
    global evt
    evt = os.path.abspath(file.name)
    evt_label.config(text="Event file: "+ evt)

'''--------------opening and storing location of mkf files----------------------------------------------------------------------------------------'''
def open_mkf_file():
    file = filedialog.askopenfile(mode='r', title="Open File")
    global mkf
    mkf = os.path.abspath(file.name)
    mkf_label.config(text='Mkf file: '+mkf)

'''--------------opening and storing location of mkf threshold files----------------------------------------------------------------------------------------'''
def open_mkf_thres_file():
    file = filedialog.askopenfile(mode='r', title="Open File")
    global mkf_threshold
    mkf_threshold = os.path.abspath(file.name)
    mkf_thres_label.config(text='Mkf Threshold file: '+mkf_threshold)

'''--------------opening and storing location of gti files----------------------------------------------------------------------------------------'''
def open_gti_file():
    file = filedialog.askopenfile(mode='r', title="Open File")
    global gti
    gti = os.path.abspath(file.name)
    gti_label.config(text='GTI file: '+gti)

'''--------------asking for the gtitype------------------------------------------------------------------------------------------------------------'''
def gti_type_common():
    global gtitype
    global gtitype_label
    gtitype = 'COMMON'
    gtitype_label.config(text='GTI Type: '+gtitype)

'''--------------opening and storing location of livetime files----------------------------------------------------------------------------------------'''
def open_livetime_file():
    file = filedialog.askopenfile(mode='r', title="Open File")
    global livetime
    livetime = os.path.abspath(file.name)
    livetime_label.config(text='Livetime file: '+livetime)

'''--------------opening and storing location of badpix files----------------------------------------------------------------------------------------'''
def open_badpix_file():
    file = filedialog.askopenfile(mode='r', title="Open File")
    global badpix
    badpix = os.path.abspath(file.name)
    badpix_label.config(text='Badpix file: '+badpix)

'''--------------opening and storing location of dpi files----------------------------------------------------------------------------------------'''
def open_dpi_file():
    file = filedialog.askopenfile(mode='r', title="Open File")
    global dpi
    dpi = os.path.abspath(file.name)
    dpi_label.config(text='DPI file: '+dpi)

'''--------------opening and storing location of dph files----------------------------------------------------------------------------------------'''
def open_dph_file():
    file = filedialog.askopenfile(mode='r', title="Open File")
    global dph
    dph = os.path.abspath(file.name)
    dph_label.config(text='DPH file: '+dph)



'''--------------opening and storing location of output files----------------------------------------------------------------------------------------'''
def output():
    global out
    global dir
    out = filedialog.askdirectory() +'/'
    dir.set(out)
    
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
    man_dpigen(dph, badpix)

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
    global dir
    auto_win = Tk() #creating the window
    auto_win.title('Automatic Pipeline')
    auto_win.geometry("1200x600")
    input_frame = Frame(auto_win) #creating the frame
    ilabel = Label(input_frame, text="Choose the required input files", font='Arial 16 bold')
    ilabel.pack(side=TOP)
    b1 = Button(input_frame, text='Browse event file', command=open_evt_file)
    b2 = Button(input_frame, text='Browse mkf file', command=open_mkf_file)
    b3 = Button(input_frame, text='Browse mkf threshold file', command=open_mkf_thres_file)
    b4 = Button(input_frame, text='Browse livetime file', command=open_livetime_file)
    dir = StringVar()
    evt_label = Label(input_frame, text="Event file: "+evt)
    mkf_label = Label(input_frame, text="Mkf file: "+mkf)
    mkf_thres_label = Label(input_frame, text="Mkf Threshold file: "+mkf_threshold)
    livetime_label = Label(input_frame, text="Livetime file: "+livetime)
    evt_label.pack(pady=10)
    mkf_label.pack(pady=10)
    mkf_thres_label.pack(pady=10)
    livetime_label.pack(pady=10)

    confirm_button = Button(auto_win, text='Confirm', command=confirm_auto)
    button_frame = Frame(auto_win)

    b1.pack(side=LEFT)
    b2.pack(side=LEFT)
    b3.pack(side=LEFT)
    b4.pack(side=LEFT)
    button_frame.pack(side=BOTTOM)
    confirm_button.pack(side=BOTTOM)
    input_frame.pack(side=LEFT)
    
    auto_win.mainloop()
automaticgui()