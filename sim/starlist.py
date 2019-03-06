from __future__ import print_function
if hasattr(__builtins__, 'raw_input'):
    input = raw_input # adjustment for python < 3.0
import numpy as np
import pandas as pd
import argparse
# coding: utf-8

### Purpose: Query Simbad to obtain ICRS J=2000 coordinates of a target and make a Keck starlist

from cpsutils.hires.exposure import exposure_time, exposure_counts

def simbad_query(query_id):
    """This is a basic simbad query that returns a pandas dataframe"""
    from astroquery.simbad import Simbad
    try: Simbad.add_votable_fields('flux(g)')
    except: pass
    try: Simbad.add_votable_fields('flux(V)')
    except: pass
    try: Simbad.add_votable_fields('flux(R)')
    except: pass
    try: Simbad.add_votable_fields('flux(J)')
    except: pass
    try: Simbad.add_votable_fields('flux(H)')
    except: pass
    try: Simbad.add_votable_fields('flux(K)')
    except: pass
    try: Simbad.add_votable_fields('ids')
    except: pass
    try: Simbad.add_votable_fields('sptype')
    except: pass
    try:
        Simbad.add_votable_fields('PMRA')
        Simbad.add_votable_fields('PMDEC')
    except:
        try:
            Simbad.add_votable_fields('pmra')
            Simbad.add_votable_fields('pmdec')
        except:
            pass
    import pandas as pd
    result_table = Simbad.query_object(query_id)
    result_table.pprint(show_unit=True)
    df = result_table.to_pandas()
    return df

def write_starlist(df):
    """Returns an ascii line that can be pasted or exported directly into a Keck starlist.
    Inputs:
        a pandas database (df) with the following columns:
            cps_id - what we will call the object
            RA_now - current RA (already advanced to include pmra)
            DEC_now - current Declination (already advanced to include pmdec)
            vmag - magnitude in V band
            counts - desired Keck HIRES counts (kcounts)
            err - desired RV error (m/s)
            exptime - expected exposure time (sec)
            maxtime - desired maximum exposure time (sec)
            decker - HIRES decker (C2, B5, B3, B1, etc.)
            i2 - iodine 'in' or 'out'
            shots - the number of consecutive shots
            prio - priority (p1, p2, or p3)
            program - the two-letter program code, or the target owner's initials
            comments - any additional comments"""

    print("\n\nKeck starlist lines:\n")
    starlist_string=""
    for i,row in df.iterrows():
        ra, dec = row.RA_now, row.DEC_now
        ra0, ra1 = ra.split('.')
        if int(ra1[1]) >= 5:
            ra_print = ra0+'.'+str(int(ra1[0]) + 1) # only leading digit in RA decimal
        else:
            ra_print = ra0+'.'+ra1[0]
        dec0, dec1 = dec.split('.')
        dec_print = dec0 # no decimal in DEC
        magstring = "vmag=%0.1f" % df.vmag[i]
        row_string = "%15s %8s %8s 2000 %9s %4.f/%4.f %4s %2s %ix %3s %2s %2s %s\n" % (row.cps_id, ra_print, dec_print, magstring, row.exptime, row.maxtime, row.counts, row.decker, row.shots, row.i2, row.prio, row.program, row.comment)
        print(row_string)
        starlist_string+=row_string
    return starlist_string

def keck_starlist_writer(cps_id=None,
                        query_id=None,
                        vmag=None,
                        err=None,
                        counts=None,
                        maxtime=None,
                        decker=None,
                        i2=None,
                        shots=None,
                        prio=None,
                        program=None,
                        comment=None
                        ):
    """This is the user-input Keck starlist writer.  Run this function if you don't know what to do!"""
#    print("You are running the Keck-HIRES starlist writer.  \nAn example fast input is:\ncps_id='K02732',query_id='KOI-2732',decker='C2',err=6,counts=30,maxtime=1800,i2='y',shots=1,program='LW',prio='p2',comment='KGPS')")
    import numpy as np
    import pandas as pd
    import argparse

    # support for command line arguments
    parser = argparse.ArgumentParser(description='Generate one line of a Keck HIRES starlist with [user-prompted] command line input.  To generate many lines, see write_starlist.py')
    parser.add_argument("--cps_id", type=str,
                    help="CPS object ID")
    parser.add_argument("--simbad_id", type=str,
                    help="a SIMBAD-resolvable object ID")
    parser.add_argument("--vmag", type=float,
                    help="apparent V magnitude")
    parser.add_argument("--err", type=float,
                    help="desired RV error (m/s; only used if no counts are supplied)")
    parser.add_argument("--counts", type=int,
                    help="desired HIRES exposure meter counts (in units of kilocounts)")
    parser.add_argument("--maxtime", type=float,
                    help="maximum exposure time (seconds)")
    parser.add_argument("--decker", type=str, choices=['C2','B5','B3','B1','C1'],
                    help="HIRES decker")
    parser.add_argument("--i2", type=bool,choices=[True,False],
                    help="HIRES i2 cell position (True = in / False = out)")
    parser.add_argument("--shots", type=int,
                    help="number of consecutive shots (max 100)")
    parser.add_argument("--prio", type=int, choices=[1,2,3],
                    help="object priority")
    parser.add_argument("--program", type=str,
                    help="object program code (Two letters, usually human initials")
    parser.add_argument("--comment", type=str,
                    help="additional comments go here")
    parser.add_argument("--advance_coords", type=bool, choices=[True,False],
                    help="argument to advance coordinates to current epoch (default=True)")
    args = parser.parse_args()

    #### WELCOME MESSAGE
    print("""Welcome to the interactive CPS starlist writer!
             This tool is customized for Keck-HIRES starlists.
             Written by L. M. Weiss, 2018""")
    parser.print_usage() # for just the usage line
    #### CPS ID ####
    if args.cps_id:
        cps_id = args.cps_id
    if cps_id == None: cps_id = str(input("What is the CPS target name?")).replace('"','').replace("'","")
    print(("cps_id: ", cps_id))

    #### SIMBAD query ID & DO QUERY ####
    if args.simbad_id:
        query_id = args.simbad_id
    if query_id == None: query_id = str(input("What is the Simbad query id?")).replace('"','').replace("'","")
    print(("query_id: ", query_id)) 
    try:
        df = simbad_query(query_id)
    except:
        raise ValueError('The argument you passed is not a resolvable SIMBAD object.  Please try a different name.  (Note: HD names like HD 189733 are SIMBAD-resolvable; but names like HD-189733 are not.  Hyphens may be resolvable in other catalogs of object names.)')
    df['cps_id'] = cps_id

    #### advance coords if possible ####
    if (args.advance_coords!=False):
        try:
            from astropy.time import Time
            ra_sex, dec_sex = advance_coords(df,Time.now().jd) # computed for J2000.0
            ra_string, dec_string = stringify_radec(ra_sex.value, dec_sex.value)
            df['RA_now'] = ra_string
            df['DEC_now'] = dec_string
            print("advanced coordinates using: \n\tpmra_cos_dec = %s mas/yr\n\tpmdec = %s mas/yr\n\tfor date %s" % (df.PMRA[0], df.PMDEC[0], Time.now()))
        except:
            print("No proper motion in Simbad or something else went wrong...not advancing coordinates")
            df['RA_now'] = df.RA
            df['DEC_now'] = df.DEC
    
    #### vmag ####
    if args.vmag:
        vmag = args.vmag
    if vmag==None:
        vmag = df.FLUX_V#.tolist()
        print("This star has magnitude V=%0.1f " % df.FLUX_V)
        if df['FLUX_V'].isnull().values.all():
            vmag = df.FLUX_g#.tolist()
            print("This star has magnitude g=%0.1f " % df.FLUX_g)
            if df['FLUX_g'].isnull().values.all():
                vmag = np.array(df.FLUX_J) + 1.5 # last ditch hack to estimate V mag
                print("This star has magnitude J=%0.1f and V approx. %0.1f" % (df.FLUX_J, df.FLUX_J + 1.5))
    df['vmag'] = vmag
    
    #### counts ####
    if args.counts:
        counts = args.counts
    while counts == None:
        counts = int(input("What is the desired number of kilo-counts (k) on the exposure meter? [Enter an int; 0=calculate based on RV err]"))   
    df['counts'] = str(counts)+'k'
    print("counts: ", counts)

    #### decker ####
    if args.decker:
        decker = args.decker
    if decker==None:
        decker = ''
    valid_deckers = ['C2','B5','B3','B1','C1']
    while decker not in valid_deckers:
        decker = str(input("Which decker would you like? [Enter C2, B5, B3, B1, or C1]")).replace('"','').replace("'","")
    print(("decker: ", decker))
    df['decker'] = decker

    #### I2 ####
    if args.i2:
        i2_bool = args.i2
    else:
        valid_i2 = ['yes','no','y','n','in','out',True,False]
        while i2 not in valid_i2:
            i2 = str(input("Would you like I2 in? [Enter y or n]")).replace('"','').replace("'","")
            if i2 ==('y' or 'in' or 'yes' or True):
                i2_bool = True
            elif i2 ==('n' or 'out' or 'no' or False):
                i2_bool = False
    i2_dict = {True:'in',False:'out'}
    i2_string = i2_dict[i2_bool]
    df['i2'] = i2_string

    #### RV error ####
    if args.err:
        err = args.err
    while ((err == None) and (counts == 0)):
        err = float(input("What RV error would you like (m/s)? [Enter a float]")) # only necessary if there is no count supplied
        try:
            counts = int(err_to_counts(err))
            df['counts'] = str(counts)+'k'
            print("counts: ", counts)
        except:
            print("Could not set counts, please try again")
            counts = 0
    if ((err == None) and (i2_bool==True)):
        print("Iodine-in RV error (scaling from 2.5 m/s at 250k for *quiet, well-behaved* stars):")
        err = counts_to_err(counts)
    if i2_bool==False:
        err = 1000.
    print("RV error (m/s): ", err)
    df['err'] = err

    #### calculate expected exposure time ####
    exptime=exposure_time(vmag, counts, iod=i2_bool, t1=110.0, v1=8.0, exp1=250.0) # compute the expected exposure time scaling from Vmag8 at 250k takes 110 seconds
    print(("Expected exposure time: ", exptime))
    df['exptime'] = exptime
    
    #### maximum time ####
    if args.maxtime:
        maxtime = args.maxtime
    if maxtime == None: maxtime = float(input("What is the maximum exposure time for this target (sec)? [Enter a float]"))
    print(("maxtime: ", maxtime))
    df['maxtime'] = maxtime

    #### shots ####
    if args.shots:
        shots = args.shots
    if shots in range(1,101):
        pass
    else:
        try:
            shots = int(input("How many (consecutive) shots would you like? [Default = 1; enter a natural number between 1 and 100]"))
        except:
            shots = 1
    print(("shots:", shots))
    df['shots'] = shots

    #### priority ####
    if args.prio:
        prio = args.prio
    valid_prio = ['p1','p2','p3','1','2','3',1,2,3]
    while prio not in valid_prio:
        prio = input("What is the program priority? [p1, p2, or p3]").replace('"','').replace("'","")
    if len(str(prio)) < 2:
        prio = "p"+str(prio)
    print("prio:", prio)
    df['prio'] = prio

    #### program ####
    if args.program:
        program = args.program
    if program==None: program = input("What is the program code? [Two-character string or human initials]").replace('"','').replace("'","")
    df['program'] = program

    #### comments ####
    if args.comment:
        comment = args.comment
    if comment==None: comment = input("Add any comments here (close companions, specific observation time, etc.):").replace('"','').replace("'","")
    df['comment'] = comment


    print("\nComparing counts and desired RV error:")
    import numpy as np
    counts_iod_in = err_to_counts(err)
    print(('%ik counts for an RV precision of %0.1f m/s' % (round(counts_iod_in,-1), err)))
    counts = round(counts, -1)
    starlist = write_starlist(df)
    return df, starlist

def err_to_counts(err):
    '''Compute exposure meter setting (in kcounts), scaling from 2.5 m/s for 250k counts'''
    from numpy import sqrt
    return (2.5/err)**2. * 250.

def counts_to_err(counts):
    '''Compute the expected RV error for an iodine-in observation, scaling from 2.5 m/s at 250k counts'''
    from numpy import sqrt
    return 2.5 / sqrt(counts/250.)

def advance_coords(df, new_time, old_time = 2451544.5):
    """A homemade coordinate advancer, since astropy code to advance coordinates does not work for python <= 2.7.
    Note: pmra given on Simbad is already multiplied by cosdec to give a space velocity, rather than an angular velocity.  we want to convert back to predict changes in RA.
    RA - in degrees
    Dec - in degrees"""
    from numpy import cos, pi, array
    from astropy.coordinates import SkyCoord, ICRS, Angle
    import astropy.units as u
    from astropy.time import Time
    rah,ram,ras = [float(s) for s in df.RA[0].split()]
    decd, decm,decs = [float(s) for s in df.DEC[0].split()]
    ra_sex = Angle([rah,ram,ras], unit=u.deg)
    dec_sex = np.array([float(decd),float(decm),float(decs)]) * u.deg
    pm_ra_cosdec = float(df.PMRA[0])*u.mas/u.yr
    pm_dec= float(df.PMDEC[0]) *u.mas/u.yr
    ra_deg, dec_deg = hms2deg(ra_sex, dec_sex)
    cosdec = cos(dec_deg * pi/180.)
    mura = pm_ra_cosdec / cosdec # in degrees
    mudec = pm_dec # in degrees
    mura = mura.to('deg/yr')
    mudec = mudec.to('deg/yr')
    delta_t = (new_time - old_time) / 365.25 # 2451544.5 for J2000.0
    ra_new_deg = ra_deg + mura * delta_t * u.year
    dec_new_deg = dec_deg + mudec * delta_t * u.year
    ra_new, dec_new = deg2hms(ra_new_deg, dec_new_deg)
    return ra_new, dec_new

def hms2deg(ra_sex, dec_sex):
    rah, ram, ras = ra_sex # unpack
    decd, decm, decs = dec_sex
    ra = rah * 15. # hours to degrees
    ra = ra + ram/60. # add minutes
    ra = ra + ras/3600.

    dec = decd
    dec = dec + decm/60.
    dec = dec + decs/3600.
    return ra, dec

def deg2hms(ra_deg, dec_deg):
    from numpy import cos, pi, array, floor, round
    from astropy.coordinates import SkyCoord, ICRS, Angle
    import astropy.units as u
    rah = int(floor(ra_deg.value) / 15) # hours not degrees
    ra_deg.value % 1
    ram = int(floor((ra_deg.value % 1) * 60))
    ras = round((((ra_deg.value % 1) * 60) % 1) * 60,2)

    decd = int(floor(dec_deg.value)) #  degrees
    ra_deg.value % 1
    decm = int(floor((dec_deg.value % 1) * 60))
    decs = round((((dec_deg.value % 1) * 60) % 1) * 60,1)
        
    ra_sex = Angle([rah,ram,ras], unit=u.deg)
    dec_sex = np.array([decd, decm, decs]) * u.deg
    return ra_sex, dec_sex

def stringify_radec(ra_sex,dec_sex):
    rah, ram, ras = ra_sex
    decd, decm, decs = dec_sex
    ra_string = "%02i %02i %05.2f" % (rah, ram, ras)
    dec_string = "%02i %02i %04.1f" % (abs(decd), decm, decs)
    if decd >= 0:
        dec_string = '+'+dec_string
    else:
        dec_string = '-'+dec_string
    return ra_string, dec_string
#    from astropy import ICRS
#    from astropy.coordinates import SkyCoord
#    from astropy.time import Time
#    coord = SkyCoord(ra,dec,pmra=pmra,pmdec=pmdec)
#    coord.apply_space_motion(new_time)

def find_hd_names(df):
    """Find the Henry Draper names from an object queried in Simbad.
       Input:
          df - a pandas dataframe containing the result of a simbad_query
       Output: the HD name of the object (if found), or the MAIN_ID otherwise"""
    hd_names = []
    for idx, row in df.iterrows():
        ids = row.IDS.decode("utf-8").split("|")
        for s in ids:
            if 'HD' in s:
                hd_name = s
            else:
                pass
        try:
            hd_names.append(hd_name)
        except:
            print("No HD name resolvable")
            pass
    return hd_names

def test_write_starlist():
    my_string1 = write_starlist('209458', 'HD 209458', 2.5, 250, 500, 'B5', 'in', 1, 'p2','LW', comment='transiting planet')
    my_string2 = write_starlist('K00277', 'Kepler-36', 4, 60, 1200, 'C2', 'in', 1, 'p2','LW', comment='Kepler-36')
    my_string3 = write_starlist('K00142', 'KOI-142', 5, 30, 1200, 'C2', 'in', 1, 'p2','LW', comment='King of TTVs')
    return


if __name__ == '__main__':
    keck_starlist_writer()
#else:
#    print("example starlists:")
#    test_write_starlist()
#    print("example command:")
#    print("cps_id='K02732',query_id='KOI-2732',decker='C2',err=6,counts=30,maxtime=1800,i2='y',shots=1,program='LW',prio='p2',comment='KGPS')")





