import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
import re
import os
import sys
import pandas as pd
import flopy.utils.binaryfile as bf
from matplotlib.colors import BoundaryNorm, Normalize, LinearSegmentedColormap
import matplotlib.ticker as ticker

'''
Readme
    - Use Python Version 3
    - To run script: use "python main.py input_file"
      for example: python main.py input/input.csv
    - Check output in the output folder.     


'''


def create_new_dir(directory):
    # directory = os.path.dirname(file_path)
    try:
        os.stat(directory)
    except:
        os.mkdir(directory)
        print(f'Created a new directory {directory}\n')


def create_output_folders():
    create_new_dir('output')
    create_new_dir('output/png')
    create_new_dir('output/shp')
    create_new_dir(f'output/png/conc_{var}')
    create_new_dir(f'output/shp/conc_{var}')


def read_ucn(ifile_ucn):
    ucnobj = bf.UcnFile(ifile_ucn, precision='double')
    times = ucnobj.get_times()
    data = ucnobj.get_alldata(mflay=None, nodata=-1)/1000
    ntimes, nlay, nr, nc = data.shape
    # times = ucnobj.get_times()
    # for t in times:
    #    conc = ucnobj.get_data(totim=t)
    return data, ntimes, nlay, nr, nc, times, ucnobj


def generate_map1(arr, ofile, ptitle, levels, colors, xy):
    '''
        - Generate 2D plume maps
        - Last updated on 3/15/2022 by hpham
    '''

    # Mapping using GeoPandas
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 8))

    #
    dx = pd.read_csv(f'input/delx.csv')
    dy = pd.read_csv(f'input/dely.csv')
    dx_mesh_edge = dx.X-dx.delx/2
    dy_mesh_edge = dy.Y+dy.dely/2

    cmap = LinearSegmentedColormap.from_list("", colors)

    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    im = ax.pcolormesh(dx_mesh_edge, dy_mesh_edge, arr,
                       cmap=cmap, norm=norm, alpha=1)

    # Plot 2 - Show Basalt Above Water Table
    show_Basalt = True
    if show_Basalt:
        ifile = f'input/shp_files/Basalt_Above_Water_Table.shp'
        layer = gpd.read_file(ifile)
        layer.plot(ax=ax, alpha=1,  color='#bdbdbd',  # linewidth=0.25,
                   edgecolor='none', zorder=1, legend=True, label='Basalt Above Water Table')

    # Plot 3 - GWIA
    show_gwia = False
    if show_gwia:
        ifile_zone = f'input/shp_files/GWIA_2017.shp'
        # ifile_zone = f'input/GIS/hss/test_domain.shp'
        zones = gpd.read_file(ifile_zone)
        zones.plot(ax=ax, alpha=0.1, linewidth=0.75, color='none',
                   edgecolor='k', zorder=1, legend=True, label='GWIA')
        # zones.apply(lambda x: ax.annotate(s=x.GWIA_NAME,
        #                                  xy=x.geometry.centroid.coords[0], ha='center'), axis=1)
    # Plot 4 - Inner_and_Outer_Boundary
    show_Inner_and_Outer_Boundary = False
    if show_Inner_and_Outer_Boundary:
        ifile_zone = f'input/shp_files/Inner_and_Outer_Boundary.shp'
        zones = gpd.read_file(ifile_zone)
        # zones.plot(ax=ax, alpha=0.25, linewidth=0.5, color='none',
        #           edgecolor='#636363', zorder=1, legend=True, label='Waste Site')  # darkred

    # Plot 5 - show former operational areas (e.g. 200W, 200E, etc.)
    show_OU = True
    if show_OU:
        ifile_zone = f'input/shp_files/bdjurdsv.shp'
        zones = gpd.read_file(ifile_zone)
        zones.plot(ax=ax, alpha=0.25, linewidth=0.5, color='none',
                   edgecolor='#636363', zorder=1, legend=True, label='Waste Site')  # darkred
        # zones.apply(lambda x: ax.annotate(s=x.NAME,
        #                                  xy=x.geometry.centroid.coords[0], ha='center'), axis=1)

    # Plot 6 - show River
    show_river = True
    if show_river:
        ifile_zone = f'input/shp_files/River.shp'
        zones = gpd.read_file(ifile_zone)
        # print(zones.head())
        zones.plot(ax=ax, alpha=0.3, linewidth=0.75, color='#2b8cbe',
                   edgecolor='#2b8cbe', zorder=2, legend=True, label='River')
    # Plot 7
    show_mcali_zone = False
    if show_mcali_zone:
        ifile_zone = f'input/shp_files/P2Rv831_focus_calibration_area.shp'
        zones = gpd.read_file(ifile_zone)
        # print(zones.head())
        zones.plot(ax=ax, alpha=0.75, linewidth=1.25, color='none',
                   edgecolor='blue', zorder=2, legend=True, label='Focused_MCali_Zone')

    # Plot 8 - show AWLN
    show_AWLN = False
    if show_AWLN:
        ifile_zone = f'input/shp_files/AWLN.shp'
        zones = gpd.read_file(ifile_zone)
        # print(zones.head())
        zones.plot(ax=ax, alpha=0.75, linewidth=0.75, color='#eff3ff',
                   edgecolor='#3182bd', zorder=3, legend=True, label='AWLN')
        zones.apply(lambda x: ax.annotate(s=x.Name,
                                          xy=x.geometry.centroid.coords[0], ha='center'), axis=1)
    # Plot 9 - Show CA Compliance boundary
    show_CA_Boundary = False
    if show_CA_Boundary:
        ifile = f'input/shp_files/CompositeAnalysis_1998.shp'
        zones = gpd.read_file(ifile)
        zones.plot(ax=ax, alpha=0.25, linewidth=1.5, color='none',
                   edgecolor='k', zorder=1, legend=True, label='CA Compliance Boundary')

    #divider = make_axes_locatable(ax)
    #cax = divider.append_axes("right", size="5%", pad=0.05)
    # fig.colorbar(im, ax=ax, format=ticker.FuncFormatter(
    #    fmt), fraction=0.046, pad=0.04)
    fig.colorbar(im, fraction=0.02, pad=0.04, format=ticker.FuncFormatter(fmt))

    #
    ax.set_title(ptitle)
    ax.set_xlim([dx.X.min(), dx.X.max()])
    ax.set_ylim([dy.Y.min(), dy.Y.max()])
    #
    ax.set_xlim([xy[0], xy[1]])  # xmin, xmax
    ax.set_ylim([xy[2], xy[3]])  # ymin, ymax

    # plt.gca().set_aspect('equal', adjustable='box')
    fig.savefig(ofile, dpi=300, transparent=False, bbox_inches='tight')
    print(f'Saved {ofile}\n')
    # plt.show()
    plt.close('all')


def conv_str2num(str):
    '''
    Convert a list of strings to list of numbers
    '''
    clevels = []
    for i in str:
        clevels.append(float(i))
    return clevels


def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)


def arr2shp(arr, ofile):
    '''
    Export an arry to a shapefile
    '''
    nr, nc = arr.shape
    print(nr, nc)
    val = np.reshape(arr, nr*nc)

    # use model grid shapefile geometry
    # model grid shapefile
    try:
        gridShp = os.path.join('input', 'shp_files', 'grid_274_geo_rc.shp')
    except:
        print('ERROR: Specify path to grid_274_geo_rc.shp')
    gdf = gpd.read_file(gridShp)
    df = pd.DataFrame()
    df['row'] = gdf['row']
    df['column'] = gdf['column']
    df['val'] = val

    # export shapefile
    gdf1 = gpd.GeoDataFrame(df, crs='EPSG:4326', geometry=gdf.geometry)
    # ofile_shp = os.path.join(
    #    work_dir, 'scripts', 'output', 'shp', 'Cmax_trit_ts76.shp')
    gdf1.to_file(driver='ESRI Shapefile', filename=ofile)
    print(f'Saved {ofile} ! ! !')


if __name__ == "__main__":

    # [1] Load input file -----------------------------------------------------
    #ifile = f'input/input.csv'
    ifile = sys.argv[1]

    # read input file ---------------------------------------------------------
    dfin = pd.read_csv(ifile)
    dfin = dfin.set_index('var')

    # Read lines in the input file --------------------------------------------
    sce = dfin['name'].loc['sce']
    var = dfin['name'].loc['var']
    ucn_file = dfin['name'].loc['ucn_file']  # full path to ucn file
    conc_cutoff = float(dfin['name'].loc['conc_cutoff'])
    contour_levels = re.split(',', dfin['name'].loc['contour_levels'])
    contour_levels = conv_str2num(contour_levels)  # convert to list of numbers
    colors = re.split(',', dfin['name'].loc['color_levels'])
    list_sp = re.split(',', dfin['name'].loc['list_sp'])
    list_sp = conv_str2num(list_sp)

    list_layer = re.split(',', dfin['name'].loc['list_layer'])
    list_layer = conv_str2num(list_layer)
    map_dim = re.split(',', dfin['name'].loc['map_dim'])
    map_dim = conv_str2num(map_dim)

    ucn2png = dfin['name'].loc['ucn2png']
    ucn2shp = dfin['name'].loc['ucn2shp']

    # [1] Map spatial distribution of total mass/activity arriving at water table
    if ucn2png == 'yes':
        '''
        This generates plume maps (in png files): 
            + for a given layers and stress periods, or 
            + for maximum plume footprint (max over all model layers)
        '''

        # Create some output folders to write outputs
        create_output_folders()

        # Read ucn file using flopy -------------------------------------------
        data, ntimes, nlay, nr, nc, times, ucnobj = read_ucn(ucn_file)
        print(f'nrow={nr}, ncol={nc}, nlay={nlay}, nsp={ntimes}\n')

        data = np.ma.masked_less_equal(data, conc_cutoff)
        Cmax_over_layer = np.nanmax(data, 1)

        vmin, vmax = np.nanmin(data), np.nanmax(data)
        print(f'Cmin={vmin}, Cmax={vmax}\n')

        for ilay in list_layer:
            for isp in list_sp:
                if ilay == 999:
                    arr = Cmax_over_layer[int(isp)-1, :, :]
                    ptitle = f'COC: {var}, Layer: Max Footprint, SP: {isp}'
                    # output png file
                    ofile_png = f'output/png/conc_{var}/Conc_{var}_Lay_max_SP_{int(isp)}.png'
                    if ucn2shp == 'yes':
                        ofile_shp = f'output/shp/conc_{var}/Conc_{var}_Lay_max_SP_{int(isp)}.shp'
                        arr2shp(arr, ofile_shp)
                        print(f'Saved {ofile_shp}\n')
                else:
                    arr = data[int(isp)-1, int(ilay)-1, :, :]
                    ptitle = f'COC: {var}, Layer: {ilay}, SP: {isp}'
                    # output png file
                    ofile_png = f'output/png/conc_{var}/Conc_{var}_Lay_{int(ilay)}_SP_{int(isp)}.png'
                    if ucn2shp == 'yes':
                        ofile_shp = f'output/shp/conc_{var}/Conc_{var}_Lay_{int(ilay)}_SP_{int(isp)}.shp'
                        arr2shp(arr, ofile_shp)
                        print(f'Saved {ofile_shp}\n')

                # Map array arr -----------------------------------------------
                generate_map1(arr, ofile_png, ptitle, contour_levels,
                              colors, map_dim)
                print(f'Saved {ofile_png}\n')
