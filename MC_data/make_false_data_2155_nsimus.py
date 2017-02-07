#! /usr/bin/env python
import sys
from gammapy.data import DataStore
import logging
import subprocess
import os
from astropy.coordinates import SkyCoord,Angle
import shutil
import hashlib
from pathlib2 import Path
import click
import numpy as np
from astropy.io import fits
from astropy.table import Table,QTable, vstack
from astropy.table import join
from astropy.time import Time
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
#./make_false_data_2155_nsimus.py "ash_stereo_thsq64" "Prod15_4_stereo" 200 "PWL" 2
#./make_false_data_2155_nsimus.py "ash_stereo_thsq64" "Prod15_4_stereo" 200 "PWL" 2.2
#./make_false_data_2155_nsimus.py "ash_stereo_thsq64" "Prod15_4_stereo" 200 "EXP" 2.3 0.2
#./make_false_data_2155_nsimus.py "ash_stereo_thsq64" "Prod15_4_stereo" 200 "EXP" 2 0.2
#./make_false_data_2155_nsimus.py "ash_stereo_thsq64" "Prod15_4_stereo" 200 "EXP" 2.1 0.2

class ObservationData:
    """Helper functions to compute file and folder names.
    """

    # filetypes = ['events', 'aeff', 'edisp', 'psf_3gauss']
    filetypes = ['events']

    def __init__(self, obs_id, hap_config=None, telpattern=None):
        self.obs_id = obs_id
        self.hap_config = hap_config
        self.telpattern = telpattern

    @property
    def obs_group(self):
        obs_id_min = self.obs_id - (self.obs_id % 200)
        obs_id_max = obs_id_min + 199
        return obs_id_min, obs_id_max

    @property
    def _obs_group_folder(self):
        return Path('run{:06d}-{:06d}'.format(self.obs_group[0], self.obs_group[1]))

    @property
    def _obs_folder(self):
        return Path('run{:06d}'.format(self.obs_id))

    def folder(self, step=None):
        """Create folder for a given step.
        """
        if step is None:
            return self._obs_group_folder / self._obs_folder
        else:
            return Path(step) / self._obs_group_folder / self._obs_folder

    def hap_filename(self, filetype):
        """Name of FITS file generated by HAP"""
        if filetype == 'events':
            return self.folder('events') / 'run_{:07d}_{}_eventlist.fits'.format(self.obs_id, self.hap_config)
            # return self.folder('events') / 'events_{:06d}.fits.gz'.format(self.obs_id)
        elif filetype == 'aeff':
            return self.folder('irfs') / 'aeff_{:06d}.fits'.format(self.obs_id)
        elif filetype == 'edisp':
            return self.folder('irfs') / 'edisp_{:06d}.fits'.format(self.obs_id)
        elif filetype == 'psf_3gauss':
            return self.folder('irfs') / 'psf_3gauss_{:06d}.fits'.format(self.obs_id)
        else:
            raise ValueError('Invalid {} {}'.format(filetype))

    def out_filename(self, filetype, dir, format='old'):
        """Name of FITS file in out folder"""
        filename = self.filename(filetype=filetype, format=format)
        return Path(dir) / filename

    def filename(self, filetype, format='old'):
        if format == 'old':
            TAGS = dict(
                events='events',
                aeff='aeff_2d',
                edisp='edisp_2d',
                psf_3gauss='psf_3gauss',
                psf_king='psf_king',
                psf_table='psf_table',
                background='bkg_offruns',
            )
        elif format == 'new':
            TAGS = dict(
                events='events',
                aeff='aeff',
                edisp='edisp',
                psf_3gauss='psf_3gauss',
                psf_king='psf_king',
                psf_table='psf_table',
                background='background',
            )

        tag = TAGS[filetype]
        if (filetype == "events"):
               filename = '{}_{:06d}.fits.gz'.format(tag, self.obs_id)
        else:
            if(self.obs_id>99999):
                filename = '{}_0{:06d}.fits'.format(tag, self.obs_id)
            else:
                filename = '{}_{:06d}.fits'.format(tag, self.obs_id)
        return self.folder() / filename

    def mkdir(self, step):
        """Make directory (parts=True, exists_ok=True)"""
        path = self.folder(step)
        if not path.exists():
            path.mkdir(parents=True)

        return path

    def check_out_files_exist(self):
        """Check if all out files exist"""
        for filetype in self.filetypes:
            filename = self.out_filename(filetype)
            if not filename.is_file():
                log.error('MISSING: {}'.format(filename))
                return False

        return True

class MC:
    """Helper functions to compute file and folder names.
    """

    # filetypes = ['events', 'aeff', 'edisp', 'psf_3gauss']
    filetypes = ['events']

    def __init__(self):
        self.zenMC = [0, 18, 26, 32, 37, 41, 46, 50, 53, 57, 60, 63, 67, 70]
        self.effMC = [50, 60, 70, 80, 90, 100]
        self.offMC = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5]
        self.gammaMC = [2]
        self.azMC = [180,0]
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
    arg = parser.parse_args()
    #import IPython; IPython.embed()
    config = arg.config
    prod = arg.prod
    Nsimus = int(arg.Nsimus)
    loi_spectrale= arg.loi_spectrale
    gamma_input=float(arg.gamma_input)
    if loi_spectrale=="EXP":
         beta_input=float(arg.beta_input)
    indir_MC = os.path.expandvars('$CALDB')+"/"+prod+"/MC_data/"+config
    indir_data= os.path.expandvars('$CALDB')+"/"+prod+"/new_data/"+config
    ds_real_data=DataStore.from_dir(indir_data)
    iobs=945
    table_obs_realdata=ds_real_data.obs_table[iobs]
    obs_id=table_obs_realdata["OBS_ID"]
    RealObservation=ObservationData(obs_id)
    effobs=table_obs_realdata["MUONEFF"]*100
    zenobs=table_obs_realdata["ZEN_PNT"]
    azobs=table_obs_realdata["AZ_PNT"]
    offobs=np.sqrt((table_obs_realdata["RA_PNT"]-table_obs_realdata["RA_OBJ"])**2+(table_obs_realdata["DEC_PNT"]-table_obs_realdata["DEC_OBJ"])**2)
    effMC=Observation.effMC[np.abs(Observation.effMC-effobs).argmin()]
    zenMC=Observation.zenMC[np.abs(Observation.zenMC-zenobs).argmin()]
    offMC=Observation.offMC[np.abs(Observation.offMC-offobs).argmin()]
    azMC=Observation.azMC[np.abs(Observation.azMC-azobs).argmin()]
    gammaMC=Observation.gammaMC[0]
    run_number=int(Observation.MCband.run_number(azMC,zenMC, offMC, gammaMC))
    events_filename_MC = Path(indir_MC) / Observation.filename('events',effMC,run_number)
    events_filename_data = Path(indir_data) / RealObservation.filename('events', format="old")
    aeff_filename = Path(indir_data) / RealObservation.filename('aeff', format="old")
    edisp_filename = Path(indir_data) / RealObservation.filename('edisp', format="old")
    psf_filename = Path(indir_data) / RealObservation.filename('psf_table', format="old")
    try:
        table_data = QTable.read(str(events_filename_data), hdu=1)
    except Exception:
        print "fits corrupted for file " + str(events_filename_data)
        exit()
    try:
        table_data_MC = QTable.read(str(events_filename_MC), hdu=1)
    except Exception:
        print "fits corrupted for file " + str(events_filename_data)
        exit()
    #Comme le premier bin ne commence pas exactement au lower edge du bin pour etre sur d avoir la bonne erngie min pour trouver le flux integre on prend le nombre d evenment simules allant du deuxieme bin non nulle au dernier bin non nul
    table_MC_simulated=QTable.read(str(events_filename_MC), hdu=5)
    i_MC=np.where(table_MC_simulated["N_EVT"]>0)[0]
    Nsimulated_firstbin_nonzeroevents=table_MC_simulated["N_EVT"][i_MC[0]].value
    Ntot=table_MC_simulated["N_EVT"].sum().value
    Nsimulated=Ntot-Nsimulated_firstbin_nonzeroevents
    #Ici correspond du coup au lower edge du bin 13 = deuxieme bin non nul
    emin=table_MC_simulated["Ebin_MIN"][i_MC[1]]
    #Ici correspond du coup au upper edge du bin 64 = dernier bin non nul
    emax=table_MC_simulated["Ebin_MAX"][i_MC[-1]]
    ikept=np.where(table_data_MC["MC_ENERGY"]>emin)
    table=table_data_MC[ikept]
    Ntrigered=len(table_data_MC)
    #Pour simuler un run de 28min
    livetime=Quantity(table_obs_realdata["LIVETIME"],"s")
    table_MC=QTable.read(str(events_filename_MC), hdu=4)
    #Quand Bruno aura changer le bug unite sera deja bien en m2 dans la table astropy...
    Amc=math.pi*(table_MC["CORE_MAX"][0]**2-table_MC["CORE_MIN"][0]**2)
    gamma_simulated=-table_MC["INDEX"][0].value

    flux_int_simulated=(Nsimulated/(Amc*livetime)).to("cm^-2 s-1")
    eref=Quantity(1,"TeV")
    if gamma_simulated!=1:
        inter=eref*((emax/eref)**(1-gamma_simulated)-(emin/eref)**(1-gamma_simulated))/((1-gamma_simulated))
    else:
        inter=np.log(emax/emin)
    if loi_spectrale=="PWL":
        #proba=(table["MC_ENERGY"]/Emin)**(-float(gamma_input))/(table["MC_ENERGY"]/Emin)**(-gamma_simulated)
        proba=(table["MC_ENERGY"]/eref)**(-float(gamma_input))/(table["MC_ENERGY"]/eref)**(-gamma_simulated)
    elif loi_spectrale=="EXP":
        #proba=((table["MC_ENERGY"]/Emin)**(-float(gamma_input))*np.exp(-Quantity(beta_input,"TeV^-1")*table["MC_ENERGY"]))/(table["MC_ENERGY"]/Emin)**(-gamma_simulated)
        proba=((table["MC_ENERGY"]/eref)**(-float(gamma_input))*np.exp(-Quantity(beta_input,"TeV^-1")*table["MC_ENERGY"]))/(table["MC_ENERGY"]/eref)**(-gamma_simulated)
    else:
        print("you didn't give a valid spectral law")
    flux_1TeV_simulated=flux_int_simulated/inter
    flux_Crab_1TeV=Quantity(3.45e-11,"cm^-2 s^-1 TeV^-1")
    flux_wanted_1TeV=Quantity(3.45e-10,"cm^-2 s^-1 TeV^-1")
    n_events_to_draw=(flux_wanted_1TeV/flux_1TeV_simulated)*proba.sum()
    percent_flux_Crab_simulated=flux_wanted_1TeV/flux_Crab_1TeV

    source_pos=SkyCoord(table_data.meta["RA_OBJ"],table_data.meta["DEC_OBJ"], unit="deg")
    events_pos=SkyCoord(table_data["RA"],table_data["DEC"], unit="deg")
    i_events_remove=np.where(events_pos.separation(source_pos)<Angle(0.2,"deg"))[0]
    table_data.remove_rows(i_events_remove)
    for ns in range(Nsimus):
        index_events_MC=np.random.choice(np.arange(0,Ntrigered,1),n_events_to_draw,replace=False,p=proba/proba.sum())
        new_events_MC_table=Table(table_data_MC[index_events_MC])
        RA_PNT=table_data.meta["RA_PNT"]
        DEC_PNT=table_data.meta["DEC_PNT"]
        ra=RA_PNT-new_events_MC_table["SKYX_RADEC"]
        dec=DEC_PNT-new_events_MC_table["SKYY_RADEC"]
        new_events_MC_table["RA"]=ra
        new_events_MC_table["DEC"]=dec
        new_data_table=vstack([new_events_MC_table,Table(table_data)])
        new_data_table.meta["DEC_OBJ"]=np.median(dec)
        new_data_table.meta["RA_OBJ"]=np.median(ra)
        if loi_spectrale=="PWL":
            outdir=indir_MC+"/MC_simus/gamma_"+str(gamma_input)+"_flux_"+str(percent_flux_Crab_simulated.value)+"Crab/Nsimus_"+str(ns)
            outdir2=indir_MC+"/MC_simus/gamma_"+str(gamma_input)+"_flux_"+str(percent_flux_Crab_simulated.value)+"Crab/Nsimus_"+str(ns)+"/"+str(RealObservation.folder())
        elif loi_spectrale=="EXP":
            outdir=indir_MC+"/MC_simus/gamma_"+str(gamma_input)+"_beta_"+str(beta_input)+"_flux_"+str(percent_flux_Crab_simulated.value)+"Crab/Nsimus_"+str(ns)
            outdir2=indir_MC+"/MC_simus/gamma_"+str(gamma_input)+"_beta_"+str(beta_input)+"_flux_"+str(percent_flux_Crab_simulated.value)+"Crab/Nsimus_"+str(ns)+"/"+str(RealObservation.folder())
        try:
            os.makedirs(outdir2)
        except Exception:
            pass
        obstable=Table.read(indir_data+"/obs-index.fits.gz")
        #import IPython; IPython.embed()
        obstable[iobs]["RA_OBJ"]= new_data_table.meta["RA_OBJ"]
        obstable[iobs]["DEC_OBJ"]= new_data_table.meta["DEC_OBJ"]
        obstable.write(outdir+"/obs-index.fits.gz", overwrite=True)
        os.system("cp "+indir_data+"/hdu-index.fits.gz "+outdir+"/hdu-index.fits.gz")
        primary_hdu=fits.open(str(events_filename_data))[0]
        hdu1=table_to_hdu(new_data_table)
        hdu2=table_to_hdu(Table.read(str(events_filename_data), hdu=2))
        hdu3=table_to_hdu(Table.read(str(events_filename_data), hdu=3))
        hdu4=table_to_hdu(Table.read(str(events_filename_data), hdu=4))
        hdulist=fits.HDUList([primary_hdu,hdu1,hdu2,hdu3,hdu4])
        hdulist.writeto(outdir+"/"+str(RealObservation.filename('events', format="old")),clobber=True)
        a=os.system("ln -s "+str(aeff_filename)+" "+outdir+"/"+str(RealObservation.filename('aeff', format="old")))
        if a==256:
            os.system("unlink "+ outdir+"/"+str(RealObservation.filename('aeff', format="old")))
            os.system("ln -s "+str(aeff_filename)+" "+outdir+"/"+str(RealObservation.filename('aeff', format="old")))
        a=os.system("ln -s "+str(edisp_filename)+" "+outdir+"/"+str(RealObservation.filename('edisp', format="old")))
        if a==256:
            os.system("unlink "+ outdir+"/"+str(RealObservation.filename('edisp', format="old")))
            os.system("ln -s "+str(edisp_filename)+" "+outdir+"/"+str(RealObservation.filename('edisp', format="old")))
        a=os.system("ln -s "+str(psf_filename)+" "+outdir+"/"+str(RealObservation.filename('psf_table', format="old")))
        if a==256:
            os.system("unlink "+ outdir+"/"+str(RealObservation.filename('psf_table', format="old")))
            os.system("ln -s "+str(psf_filename)+" "+outdir+"/"+str(RealObservation.filename('psf_table', format="old")))
