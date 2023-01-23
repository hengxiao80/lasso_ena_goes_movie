from datetime import datetime, timedelta
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import metpy  
import numpy as np
import s3fs
import glob
import os
import imageio

def download_daily_data(year, month, day, 
                        dir='data', sat='goes16', domain='F', channels=['C01', 'C02', 'C03']):
    hours = np.arange(0, 24, 1)
    for hour in hours:
        t = datetime(year, month, day, hour)
        download_hourly_data(t, dir, sat, domain, channels)

def download_hourly_data(t, dir='data', sat='goes16', domain='F', channels=['C01', 'C02', 'C03']):
    '''
    Download 1 hour worth of GOES-R Radiance data from AWS.
    '''
    fs = s3fs.S3FileSystem(anon=True)
    year = t.timetuple().tm_year
    doy = t.timetuple().tm_yday
    hour = t.timetuple().tm_hour

    # List contents of GOES-16 bucket.
    files = fs.ls(f's3://noaa-{sat}/ABI-L1b-Rad{domain}/{year}/{doy:03d}/{hour:02d}/')

    files2down = []
    for file in files:
        for c in channels:
            if c in str(file):
                files2down.append(file)
                print(file)
                fs.get(file, dir+'/'+file.split('/')[-1])

# Define the rebin function that will be used to resample the band resolution
# Rebin function from https://stackoverflow.com/questions/8090229/resize-with-averaging-or-rebin-a-numpy-2d-array
def rebin(a, shape):
    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
    return a.reshape(sh).mean(-1).mean(1)

def make_goes16_ena_daily_plots(year, month, day, datadir='data', figdir='figs'):
    start = datetime(year,month,day,0,0)
    end = datetime(year,month,day,23,59)
    dt = timedelta(minutes=15) 
    while start < end:
        plot_goes16_ena(start, datadir=datadir, figdir=figdir)
        start += dt

def make_goes16_ena_daily_movie(year, month, day, dir='figs'):
    filenames = glob.glob(f'{dir}/GOES16_truecolor_ENA_10deg_{year}{month:02d}{day:02d}*.png')
    print(sorted(filenames))
    with imageio.get_writer(f'{year}{month:02d}{day:02d}.mp4', mode='I', fps=5, macro_block_size=4, format='ffmpeg') as writer:
        for filename in sorted(filenames):
            image = imageio.imread(filename)
            writer.append_data(image)

def plot_goes16_ena(t, datadir='data', figdir='figs', show_island=False):
    '''
    plot a GOES16 true color image around the ENA site
    Current extent is set to [-33, -23, 34, 44].
    Code based on the following:

    '''

    extent = [-33, -23, 34, 44]

    year = t.timetuple().tm_year
    doy = t.timetuple().tm_yday
    hour = t.timetuple().tm_hour
    minute = t.timetuple().tm_min
    tstamp = f'{year:04d}{doy:03d}{hour:02d}{minute:02d}'
    tstamp2 = t.strftime("%m/%d/%Y, %H:%M")
    tstamp3 = t.strftime("%Y%m%d%H%M")
    figfile = f'{figdir}/GOES16_truecolor_ENA_10deg_{tstamp3}.png'
    if os.path.exists(figfile):
        print(f'{figfile} already there!')
        return
    if len(glob.glob(f'{datadir}/OR*-M3C01_G16_s{tstamp}*.nc')) == 1:
        c01 = xr.open_dataset(glob.glob(f'{datadir}/OR*-M3C01_G16_s{tstamp}*.nc')[0])
    else:
        print(f'{tstamp} C01 file not found!')
        return
    if len(glob.glob(f'{datadir}/OR*-M3C02_G16_s{tstamp}*.nc')) == 1:
        c02 = xr.open_dataset(glob.glob(f'{datadir}/OR*-M3C02_G16_s{tstamp}*.nc')[0])
    else:
        print(f'{tstamp} C02 file not found!')
        return
    if len(glob.glob(f'{datadir}/OR*-M3C03_G16_s{tstamp}*.nc')) == 1:
        c03 = xr.open_dataset(glob.glob(f'{datadir}/OR*-M3C03_G16_s{tstamp}*.nc')[0])
    else:
        print(f'{tstamp} C02 file not found!')
        return
    # subsetting the input data manually to make plotting easier.
    R = c02['Rad'].data[1128*2:2537*2, 8000*2:10000*2]
    G = c03['Rad'].data[1128:2537, 8000:10000]
    B = c01['Rad'].data[1128:2537, 8000:10000]
    dat = c01.metpy.parse_cf('Rad')
    goes = dat.metpy.cartopy_crs
    x = dat.x[8000:10000]
    y = dat.y[1128:2537]

    kappa_B = c01['kappa0'].data
    kappa_R = c02['kappa0'].data
    kappa_G = c03['kappa0'].data
    # To convert radiance to reflectance, use formula:
    # reflectance (ρf(υ)) = kappa factor(κ) * radiance (L(ν))
    # Source: GOES-R Series Product Definition and User Guide (PUG) Volume 3, Revision 2.2, pages 27-28
    R_ref = kappa_R * R
    G_ref = kappa_G * G  
    B_ref = kappa_B * B 
    # Apply range limits for each channel. Reflectance values must be between 0 and 1.
    R_ref = np.clip(R_ref, 0, 1)
    G_ref = np.clip(G_ref, 0, 1)
    B_ref = np.clip(B_ref, 0, 1)
    # Apply a gamma correction to the image to correct ABI detector brightness
    gamma = 2.2
    R = np.power(R_ref, 1/gamma)
    G = np.power(G_ref, 1/gamma)
    B = np.power(B_ref, 1/gamma)
    # Resample the Red Band resolution
    R_rescaled = rebin(R, G.shape) 
    # GOES-R Series satellites do not have a channel in the visible green range. Band 3 is a NIR channel typically used to monitor vegetation.
    # Calculate the "True" Green Band to serve as a green proxy for the RGB True Color image, using a fractional combination.
    # Source: "Generation of GOES‐16 True Color Imagery without a Green Band" - https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2018EA000379
    G_true = 0.45 * R_rescaled + 0.1 * G + 0.45 * B
    G_true = np.clip(G_true, 0, 1)  # Apply band limits again, just in case.
    # The RGB array for the true color image
    RGB = np.dstack([R_rescaled, G_true, B])

    fig = plt.figure(figsize=(10, 8))
    pc = ccrs.PlateCarree()
    # Generate an Cartopy projection
    lc = ccrs.LambertConformal(central_longitude=-28.0, standard_parallels=(39, 39))
    ax = fig.add_axes((0.05, 0.05, 0.9, 0.9), projection=lc)
    ax.set_extent(extent, crs=pc)
    ax.imshow(RGB, origin='upper', extent=(x.min(), x.max(), y.min(), y.max()), transform=goes)
    # Add Coastlines and States
    # ax.coastlines(resolution='50m', color='green', linewidth=0.25)
    if show_island:
        ax.add_feature(ccrs.cartopy.feature.STATES, linewidth=0.25, color='green')
    # ENA location
    ax.plot(-28.025, 39.092, marker='+', color='green', markersize=10) 

    ax.set_title(f'{extent[0]}~{extent[1]}E, {extent[2]}~{extent[3]}N', loc='left', fontsize=12)
    ax.set_title(f'{tstamp2}', loc='right')
    fig.savefig(f'{figdir}/GOES16_truecolor_ENA_10deg_{tstamp3}.png')
    plt.close()

