#! /usr/bin/env python
import sys
import logging
import subprocess
import os
import shutil
import hashlib
from pathlib2 import Path
import click
import numpy as np
from astropy.io import fits
from astropy.table import Table, QTable
from astropy.table import join
from astropy.time import Time
from astropy.coordinates import Angle
import os
import argparse
import FrenchMcBands
from astropy.units import Quantity
from random import *
from astropy.io.fits import table_to_hdu
from astropy.io import fits
import math
from astropy.table import Column


#argument: config, prog,nombre de simus, loi spectrale: PWL ou EXP, si PWL le gamma voulu, si EXP le gamma et le beta
#./make_false_data_nsimus.py  "ash_stereo_thsq64" "Prod15_4_stereo" 200 "PWL" 2 1000
#./make_false_data_nsimus.py  "ash_stereo_thsq64" "Prod15_4_stereo" 200 "PWL" 2.2 1000
#./make_false_data_nsimus.py "ash_stereo_thsq64" "Prod15_4_stereo" 200 "EXP" 2.3 0.2 10000
#./make_false_data_nsimus.py "ash_stereo_thsq64" "Prod15_4_stereo" 200 "EXP" 2.1 0.2 10000
#./make_false_data_nsimus.py "ash_stereo_thsq64" "Prod15_4_stereo" 200 "EXP" 2 0.2 10000


class MC:
    """Helper functions to compute file and folder names.
    """

    # filetypes = ['events', 'aeff', 'edisp', 'psf_3gauss']
    filetypes = ['events']

    def __init__(self):
        self.zenMC = [0]
        self.effMC = [80,90]
        self.offMC = [0.5]
        self.gammaMC = [2]
        self.azMC = [180]
        self.MCband=FrenchMcBands.FrenchMcBands(i=4,energy_type="spectrum")

    def folder(self,eff,run_number):
        return Path(str(eff)+"/run"+str(run_number))

    def out_filename(self, filetype,eff,run_number):
        """Name of FITS file in out folder"""
        filename = self.filename(filetype=filetype,eff=eff,run_number=run_number)
        # return Path(dir) / filename
        return filename

    def filename(self, filetype,eff,run_number):

        TAGS = dict(
            events='events',
            aeff='aeff_2d',
            edisp='edisp_2d',
            psf_3gauss='psf_3gauss',
            psf_king='psf_king',
            psf_table='psf_table',
            background='bkg_offruns',
        )

        tag = TAGS[filetype]

        if (filetype == "events"):
            filename = Path("run_"+str(run_number)+"_eventlist.fits")
        else:
            filename = Path(tag+'_'+str(run_number)+'.fits')

        return self.folder(eff=eff,run_number=run_number) / filename



if __name__ == '__main__':
    Observation=MC()
    parser = argparse.ArgumentParser(description='Make the index and obs table')
    parser.add_argument('config', action="store", help='Prod Configuration, cut we apply')
    parser.add_argument('prod', action="store", help='Prod')
    parser.add_argument('Nsimus', action="store", help='Number of MC data we want to perform')
    parser.add_argument('loi_spectrale', action="store", help='Loi spectrale PWL or EXP')
    parser.add_argument('gamma_input', action="store", help='Indice spectrale de la power law')
    #import IPython; IPython.embed()
    if parser.parse_known_args()[0].loi_spectrale== "EXP":
        parser.add_argument('beta_input', action="store", help='coupure de la loi de puissance avec exponentielle cutoff')
    parser.add_argument('n_events_to_draw', action="store", help='Number of events to draw from the trigger list')
    arg = parser.parse_args()
    config = arg.config
    prod = arg.prod
    Nsimus = int(arg.Nsimus)
    loi_spectrale= arg.loi_spectrale
    gamma_input=float(arg.gamma_input)
    if loi_spectrale=="EXP":
         beta_input=float(arg.beta_input)
    indir = os.path.expandvars('$CALDB')+"/"+prod+"/MC_data/"+config
    for (ieff, eff) in enumerate(Observation.effMC):
        for (ioff, off) in enumerate(Observation.offMC):
            for (izen, zen) in enumerate(Observation.zenMC):
                for (igam, gamma) in enumerate(Observation.gammaMC):
                    for (iaz, az) in enumerate(Observation.azMC):
                        run_number=int(Observation.MCband.run_number(az,zen, off, gamma))
                        events_filename = Path(indir) / Observation.filename('events',eff,run_number)
                        aeff_filename = Path(indir) / Observation.filename('aeff',eff,run_number)
                        edisp_filename = Path(indir) / Observation.filename('edisp',eff,run_number)
                        psf_filename = Path(indir) / Observation.filename('psf_table',eff,run_number)
                        try:
                            table = QTable.read(str(events_filename), hdu=1)
                        except Exception:
                            print "fits corrupted for file " + str(events_filename)
                            continue

                        #Comme le premier bin ne commence pas exactement au lower edge du bin pour etre sur d avoir la bonne erngie min pour trouver le flux integre on prend le nombre d evenment simules allant du deuxieme bin non nulle au dernier bin non nul
                        table_MC_simulated=QTable.read(str(events_filename), hdu=5)
                        i_MC=np.where(table_MC_simulated["N_EVT"]>0)[0]
                        Nsimulated_firstbin_nonzeroevents=table_MC_simulated["N_EVT"][i_MC[0]].value
                        Ntot=table_MC_simulated["N_EVT"].sum().value
                        Nsimulated=Ntot-Nsimulated_firstbin_nonzeroevents
                        #Ici correspond du coup au lower edge du bin 13 = deuxieme bin non nul
                        emin=table_MC_simulated["Ebin_MIN"][i_MC[1]]
                        #Ici correspond du coup au upper edge du bin 64 = dernier bin non nul
                        emax=table_MC_simulated["Ebin_MAX"][i_MC[-1]]
                        ikept=np.where(table["MC_ENERGY"]>emin)
                        table=table[ikept]
                        Ntrigered=len(table)
                        #Pour simuler un run de 28min
                        livetime=Quantity(10,"hr")
                        #livetime=Quantity(100,"hr")
                        #livetime=table.meta["LIVETIME"]
                        table_MC=QTable.read(str(events_filename), hdu=4)
                        #Quand Bruno aura changer le bug unite sera deja bien en m2 dans la table astropy...
                        Amc=math.pi*(table_MC["CORE_MAX"][0]**2-table_MC["CORE_MIN"][0]**2)
                        gamma_simulated=-table_MC["INDEX"][0].value
                        n_events_to_draw=int(arg.n_events_to_draw)
                        #n_events_to_draw=19000
                        flux_int_simulated=(Nsimulated/(Amc*livetime)).to("cm^-2 s-1")
                        eref=Quantity(1,"TeV")
                        if gamma_simulated!=1:
                            inter=eref*((emax/eref)**(1-gamma_simulated)-(emin/eref)**(1-gamma_simulated))/((1-gamma_simulated))
                        else:
                            inter=np.log(emax/emin)
                        #import IPython; IPython.embed()
                        if loi_spectrale=="PWL":
                            #proba=(table["MC_ENERGY"]/Emin)**(-float(gamma_input))/(table["MC_ENERGY"]/Emin)**(-gamma_simulated)
                            proba=(table["MC_ENERGY"]/eref)**(-float(gamma_input))/(table["MC_ENERGY"]/eref)**(-gamma_simulated)
                        elif loi_spectrale=="EXP":
                            #proba=((table["MC_ENERGY"]/Emin)**(-float(gamma_input))*np.exp(-Quantity(beta_input,"TeV^-1")*table["MC_ENERGY"]))/(table["MC_ENERGY"]/Emin)**(-gamma_simulated)
                            proba=((table["MC_ENERGY"]/eref)**(-float(gamma_input))*np.exp(-Quantity(beta_input,"TeV^-1")*table["MC_ENERGY"]))/(table["MC_ENERGY"]/eref)**(-gamma_simulated)
                        else:
                            break
                        flux_1TeV_simulated=(flux_int_simulated/inter) *(n_events_to_draw/proba.sum())
                        flux_1Crab_1TeV=Quantity(3.45e-11,"cm^-2 s^-1")
                        percent_flux_Crab_simulated=flux_1TeV_simulated/flux_1Crab_1TeV

                        #import IPython; IPython.embed()
                        for ns in range(Nsimus):
                            index_events=np.random.choice(np.arange(0,Ntrigered,1),n_events_to_draw,replace=False,p=proba/proba.sum())
                            new_events_table=Table(table[index_events])
                            RA_PNT=Angle(new_events_table.meta["RA_PNT"],"deg")
                            DEC_PNT=Angle(new_events_table.meta["DEC_PNT"],"deg")
                            ra=RA_PNT+new_events_table["SKYX_RADEC"]
                            dec=DEC_PNT+new_events_table["SKYY_RADEC"]
                            new_events_table["RA"]=ra.value
                            new_events_table["DEC"]=dec.value
                            new_events_table.meta["LIVETIME"]=livetime.to("s").value
                            new_events_table.meta["ONTIME"]=livetime.to("s").value
                            new_events_table.meta["DEC_OBJ"]=np.median(dec).value
                            new_events_table.meta["RA_OBJ"]=np.median(ra).value
                            primary_hdu=fits.open(str(events_filename))[0]
                            hdu1=table_to_hdu(new_events_table)
                            hdu2=table_to_hdu(Table.read(str(events_filename), hdu=2))
                            hdu3=table_to_hdu(Table.read(str(events_filename), hdu=3))
                            hdu4=table_to_hdu(Table.read(str(events_filename), hdu=4))
                            hdu5=table_to_hdu(Table.read(str(events_filename), hdu=5))
                            hdulist=fits.HDUList([primary_hdu,hdu1,hdu2,hdu3,hdu4,hdu5])
                            if loi_spectrale=="PWL":
                                outdir=indir+"/gamma_"+str(gamma_input)+"_nevents_"+str(n_events_to_draw)+"_time"+str(livetime.value)+"_Nsimus_"+str(ns)
                                outdir2=indir+"/gamma_"+str(gamma_input)+"_nevents_"+str(n_events_to_draw)+"_time"+str(livetime.value)+"_Nsimus_"+str(ns)+"/"+str(Observation.folder(eff,run_number))
                            elif loi_spectrale=="EXP":
                                outdir=indir+"/gamma_"+str(gamma_input)+"_beta_"+str(beta_input)+"_nevents_"+str(n_events_to_draw)+"_time"+str(livetime.value)+"_Nsimus_"+str(ns)
                                outdir2=indir+"/gamma_"+str(gamma_input)+"_beta_"+str(beta_input)+"_nevents_"+str(n_events_to_draw)+"_time"+str(livetime.value)+"_Nsimus_"+str(ns)+"/"+str(Observation.folder(eff,run_number))
                            try:
                                os.makedirs(outdir2)
                            except Exception:
                                pass
                            #import IPython; IPython.embed()
                            obstable=Table.read(indir+"/obs-index.fits.gz")
                            obstable["LIVETIME"]=livetime.to("s").value
                            obstable["ONTIME"]=livetime.to("s").value
                            obstable["RA_OBJ"]= new_events_table.meta["RA_OBJ"]
                            obstable["DEC_OBJ"]= new_events_table.meta["DEC_OBJ"]
                            obstable.write(outdir+"/obs-index.fits.gz", overwrite=True)
                            os.system("cp "+indir+"/hdu-index.fits.gz "+outdir+"/hdu-index.fits.gz")
                            hdulist.writeto(outdir+"/"+str(Observation.filename('events',eff,run_number)),clobber=True)
                            a=os.system("ln -s "+str(aeff_filename)+" "+outdir+"/"+str(Observation.filename('aeff',eff,run_number)))
                            if a==256:
                                os.system("unlink "+ outdir+"/"+str(Observation.filename('aeff',eff,run_number)))
                                os.system("ln -s "+str(aeff_filename)+" "+outdir+"/"+str(Observation.filename('aeff',eff,run_number)))
                            a=os.system("ln -s "+str(edisp_filename)+" "+outdir+"/"+str(Observation.filename('edisp',eff,run_number)))
                            if a==256:
                                os.system("unlink "+ outdir+"/"+str(Observation.filename('edisp',eff,run_number)))
                                os.system("ln -s "+str(edisp_filename)+" "+outdir+"/"+str(Observation.filename('edisp',eff,run_number)))
                            a=os.system("ln -s "+str(psf_filename)+" "+outdir+"/"+str(Observation.filename('psf_table',eff,run_number)))
                            if a==256:
                                 os.system("unlink "+ outdir+"/"+str(Observation.filename('psf_table',eff,run_number)))
                                 os.system("ln -s "+str(psf_filename)+" "+outdir+"/"+str(Observation.filename('psf_table',eff,run_number)))
