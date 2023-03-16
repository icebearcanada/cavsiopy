#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 11 19:04:18 2021
Last change: 21 OCT 2021
This module contains miscellanous programs for:
1. legend insertion to figures (put_legend)
2. combining images (combine_horizontal, combine_vertical)
3. determine background colors for line plots
4. marking specific times on map(mark_onmap) and altitude plots(mark_altitudeplots)
5. annotating selected points: annotate_onmap, annotate_3D
6. finding index of desired values in parameters: find_index
7. finding where magnetic dip equals elevation angle
@author: ceren
"""
# =============================================================================
# insert figure legends
# =============================================================================
def put_legend(ax, n, location, labelspace, anchorx, anchory,
                labels={}, linestyles={}, markers={}, colors={}, edgecolors={}):

    """
    Inserts legend to figure for cases when plt.legend() cannot insert all
    legends or there are customized symbols.

    Parameters
    ----------
    ax: axis
    n: how many columns should the legend have
    location: location of the legend (upper right, lower left, exc.)
    labelspace: space between the labels
    anchorx, anchory: coordinate of the lower left corner of legend box

    The following are for the custamization of the symbols we want to display:
    labels: labels to display in the legend
    linestyles: linestyles used in plotting the parameters
    markers: markers used in plotting the parameters
    colors: colors used in plotting the parameters
    edgecolors: edgecolors used in plotting the parameters

    Returns
    -------
    None
    Inserts legend on plot
    """
    from matplotlib.legend_handler import HandlerBase
    class MarkerHandler(HandlerBase):
        def create_artists(self, legend, tup, xdescent, ydescent,
                            width, height, fontsize, trans):
            return [ax.scatter([width/2], [height/2],
                           marker=tup[2], color=tup[1], edgecolor=tup[0],
                           s=80, transform=trans)]

    # just use the default location if there is no anchorx and anchory
    if anchorx == None and anchory == None:
        ax.legend(list(zip(edgecolors, colors, markers)), labels,
                  handler_map={tuple:MarkerHandler()}, fontsize=12,
                  labelspacing=labelspace, handletextpad=0.02, handlelength=2,
                  ncol=n, loc=location, frameon = False)

    else:
        ax.legend(list(zip(edgecolors, colors, markers)), labels,
                  handler_map={tuple:MarkerHandler()}, fontsize=12,
                  labelspacing=labelspace, handlelength=2,
                  handletextpad=0.02, ncol=n, loc=location,
                  bbox_to_anchor=(anchorx,anchory), frameon = False)

    return
# =============================================================================
# insert figure legends
# =============================================================================
def put_legend_fnt(ax, n, location, labelspace, anchorx, anchory, fontsize,
                labels={}, linestyles={}, markers={}, colors={}, edgecolors={}):

    """
    Inserts legend to figure for cases when plt.legend() cannot insert all
    legends or there are customized symbols.

    Parameters
    ----------
    ax: axis
    n: how many columns should the legend have
    location: location of the legend location of the legend (upper right, lower left, exc.)
    labelspace: space between the labels
    anchorx, anchory: coordinate of the lower left corner of legend box
    fontsize: legend font size

    The following are for the custamization of the symbols we want to display:
    labels: labels to display in the legend
    linestyles: linestyles used in plotting the parameters
    markers: markers used in plotting the parameters
    colors: colors used in plotting the parameters
    edgecolors: edgecolors used in plotting the parameters

    Returns
    -------
    None
    Inserts legend on plot
    """
    from matplotlib.legend_handler import HandlerBase
    class MarkerHandler(HandlerBase):
        def create_artists(self, legend, tup, xdescent, ydescent,
                            width, height, fontsize, trans):
            return [ax.scatter([width/2], [height/2],
                           marker=tup[2], color=tup[1], edgecolor=tup[0],
                           s=80, transform=trans)]

    # just use the default location if there is no anchorx and anchory
    if anchorx == None and anchory == None:
        ax.legend(list(zip(edgecolors, colors, markers)), labels,
                  handler_map={tuple:MarkerHandler()}, fontsize=fontsize,
                  labelspacing=labelspace, handletextpad=0, handlelength=2,
                  ncol=n, loc=location, frameon = False)

    else:
        ax.legend(list(zip(edgecolors, colors, markers)), labels,
                  handler_map={tuple:MarkerHandler()}, fontsize=fontsize,
                  labelspacing=labelspace, handlelength=2,
                  handletextpad=0, ncol=n, loc=location,
                  bbox_to_anchor=(anchorx,anchory),
                  frameon = False)

    return
# =============================================================================
# insert figure legends
# =============================================================================
def put_legend_vectors(ax, n, location, labelspace, anchorx, anchory, fontsize,
                labels={}, colors={}, markers={}):

    """
    Inserts legend to figure for cases when plt.legend() cannot insert all
    legends or there are customized symbols. This function does not have the
    linestyles and edgecolors options.

    Parameters
    ----------
    ax: axis
    n: how many columns should the legend have
    location: location of the legend
    labelspace: space between the labels
    anchorx, anchory: coordinate of the lower left corner of legend box

    The following are for the custamization of the symbols we want to display:
    labels: labels to display in the legend
    markers: markers used in plotting the parameters
    colors: colors used in plotting the parameters

    Returns
    -------
    None
    Inserts legend on plot
    """
    from matplotlib.legend_handler import HandlerBase
    class MarkerHandler(HandlerBase):
        def create_artists(self, legend, tup, xdescent, ydescent,
                            width, height, fontsize, trans):
            return [ax.scatter([width/2], [height/2],
                           marker=tup[1], color=tup[0],
                           s=80, transform=trans)]

    # just use the default location if there is no anchorx and anchory
    if anchorx == None and anchory == None:
        ax.legend(list(zip(colors, markers)), labels,
                  handler_map={tuple:MarkerHandler()}, fontsize=fontsize,
                  labelspacing=labelspace, handletextpad=0, handlelength=2,
                  ncol=n, loc=location, frameon = False)

    else:
        ax.legend(list(zip(colors, markers)), labels,
                  handler_map={tuple:MarkerHandler()}, fontsize=fontsize,
                  labelspacing=labelspace, handlelength=2,
                  handletextpad=0, ncol=n, loc=location,
                  bbox_to_anchor=(anchorx,anchory),
                  frameon = False)

    return

# =============================================================================
# combine images
# =============================================================================
import PIL
def combine_horizontal(f, fnew):
    """
    Horizontally combines images with arbitrary image size as a new image and
    saves the new image with the given file name.

    Parameters
    ----------
    f: Name of image files, list
    fnew: Output file name

    Returns
    -------
    Saves horizontally combined image as fnew to the directory
    Reference: https://stackoverflow.com/a/30228789/15165141"""

    list_im = f
    # open all images in the list
    imgs    = [ PIL.Image.open(i) for i in list_im ]
    # find the minimum size in given images
    min_shape = sorted( [(np.sum(i.size), i.size ) for i in imgs])[0][1]
    # resize all images to match the image with minimum size
    # and stack horizontally
    imgs_comb = np.hstack([np.asarray(i.resize(min_shape)) for i in imgs])

    # Create an image memory from an object
    imgs_comb = PIL.Image.fromarray(imgs_comb)
    # save new image
    imgs_comb.save(fnew)

def combine_vertical(f, fnew):
    """
    Vertically combines images with arbitrary image size as a new image and
    saves the new image with the given file name.

    Parameters
    ----------
    f: Name of image files, list
    fnew: Output file name

    Returns
    -------
    Saves vertically combined image as fnew to the directory
    Reference: https://stackoverflow.com/a/30228789/15165141 """

    list_im = f
    # open all images in the list
    imgs    = [ PIL.Image.open(i) for i in list_im ]
    # find the minimum size in given images
    min_shape = sorted( [(np.sum(i.size), i.size ) for i in imgs])[0][1]
    # resize all images according to the minimum size
    imgs_comb = np.vstack([np.asarray(i.resize(min_shape)) for i in imgs])

    # Create an image memory from an object
    imgs_comb = PIL.Image.fromarray(imgs_comb)
    # save new image
    imgs_comb.save(fnew)

# =============================================================================
# find the indices of close values for specific angle values
# =============================================================================
def find_index(data, closeto, tol):
    """
    Finds the index of numbers that are close to the value we are searching for
    with a tolerance 'tol'.

    Parameters
    ----------
    data: input data for inspection, array
    closeto: the value we are searching for
    tol: amount of tolerance

    Returns
    -------
    Index of close values in provided data
    """

    # initialize the index array of close points
    ind_close=[]

    # are there any close numbers within a tolerance of tol?
    boolean = np.isclose(closeto, data, atol=tol)
    # no close values
    if np.all(boolean==0):
        print("".join(['No values found close to', str(closeto)]))
        ind_close=np.empty(0)

    # if there are close values
    else:
        # find how many close values
        how_many_close_vals=np.count_nonzero(boolean)
        # or np.sum(boolean!=0)
        # append the index of close values in ind_close
        ind_close=np.where(boolean==1)
        print(' '.join(['found', str(how_many_close_vals),\
                        'values close to', str(closeto), 'within',\
                        str(tol), 'tolerance.\n']))
        # col gives the number of close values
        row, col = np.shape(ind_close)
        # print which values
        for i in range (0, col):
            ind=ind_close[row-1][i]
            print('Those are\n',''.join([str(i+1),': data[',\
                                          str(ind),'] = ', str(data[ind]),'\n']))

    return ind_close

# =============================================================================
# find the indice when magnetic dip angle is very close to elevation angle
# =============================================================================
def find_dip_equal_elev(data, closeto, tol, anglename1, anglename2):
    """ index(data, closeto, tol)

    Parameters
    ----------
    data: input data for inspection, array
    closeto: the value we are searching for
    tol: amount of tolerance
    anglename1, anglename2: name of angles you want to check

    Returns
    -------
    Index of close values in provided data
    """

    ind_close = []
    # are there any close numbers within a tolerance of tol?
    ind_close = np.where(np.isclose(closeto, data, atol=tol))

    if np.size(ind_close)>0:
        how_many_close_vals=np.size(ind_close)
        print(' '.join(['found', str(how_many_close_vals),\
                        anglename1,'angle values close to',
                        anglename2, 'angle within',\
                        str(tol), 'tolerance.\n']))

    else:
        print(" ".join(['No close values between',
                        anglename1,'and', anglename2,'angle\n',
                        'within', str(tol), 'degree tolerance']))
        ind_close=np.empty(0)

    return ind_close

# =============================================================================
# mark the beginning of the pass and insert the location of the selected city
# to the plot
# =============================================================================
import numpy as np

def mark_onmap(ax, px, py, color1, z, y, x, color2):
    """ marks the beginning of the pass and location of Ottawa on the plots.

    Parameters
    ----------
    ax: axes
    px, py: Coordinates of city in lon, lat
    x: in the form: lon[0], y: lat[0]
        = first coordinate data points of pass

    Returns
    -------
    None

    Example
    -------
    mark_init(ax, px, py, x, y)

    """

    import cartopy.crs as ccrs

    # add desired ground location to plot
    ax.scatter(px, py, color='cyan', edgecolor = color1, marker='*', s=150,\
      transform=ccrs.PlateCarree(), zorder= 20)

    # mark the beginning of the pass
    arrow_plotter(ax, z, y, x)
    # ax.scatter(x, y, color=color2, marker='x', linewidth=2,
    #             alpha=1, s=100, transform=ccrs.PlateCarree(), zorder=30)

    return

def mark_altitude_plots(ax, px, py, color1, x, y, color2):
    """ Marks the beginning of the pass and location of Ottawa on the plots.

    Parameters
    ----------
    ax: axes
    px, py: Coordinates of the city
    x, y:  can be any of the lat[0], lon[0], alt[0]
        = first coordinate data points of pass

    Returns
    -------
    None

    """

    # add desired ground location to plot
    ax.scatter(px, py,  color='cyan', edgecolor = color1,
               marker='*', s=150, zorder= 20)

    # mark the beginning of the pass
    # ap(ax, z, y, x)
    ax.scatter(x, y, color=color2, marker='x', linewidth=2,
                alpha=1, s=100, zorder=30)

    return

# =============================================================================
# mark the coordinates of the points in lat-long for specific angle values
# =============================================================================

def annotate_onmap(ax, x, y, ind, val, which_angle, symbol, color, s):
    """
    Marks specific values of angles over map plots.

    Parameters
    ----------
    ax: axes to use
    x, y = coordinates to mark: float, array
    ind: index values of close points: float, array
    val: the specific value we are searching for in the data

    Returns
    -------
    None
    """
    import cartopy.crs as ccrs

    if np.size(ind)>0:
        row, col = np.shape(ind)
        # if there is more than one close value
        if col>1:
            diff_ind= np.diff(ind)
            if (diff_ind[row-1][0]==1):
                ix=ind[row-1][0]
                ax.scatter(x[ix], y[ix], marker = symbol, color= color,\
                        alpha=1, s=s, transform=ccrs.PlateCarree(), zorder= 20)
            if (diff_ind[row-1][-1]==1):
                ix=ind[row-1][-1]
                ax.scatter(x[ix], y[ix], marker = symbol, color= color,\
                        alpha=1, s=s, transform=ccrs.PlateCarree(), zorder= 20)
            else:
                for i in range(1,len(diff_ind)-1):
                    if (diff_ind[row-1][i]==1 and diff_ind[row-1][i+1]!=1):
                        ix=ind[row-1][i+1]
                        ax.scatter(x[ix], y[ix], marker = symbol, color= color,\
                                  alpha=1, s=s, transform=ccrs.PlateCarree(), zorder= 20)
                    elif (diff_ind[row-1][i+1]==1 and diff_ind[row-1][i]!=1):
                        ix=ind[row-1][i+2]
                        ax.scatter(x[ix], y[ix], marker = symbol, color= color,\
                                 alpha=1, s=s, transform=ccrs.PlateCarree(), zorder= 20)

        # if there is only one close value
        else:
            ix=ind[row-1][0]
            ax.scatter(x[ix], y[ix], marker =  symbol, color= color,\
              alpha=1, s=s, transform=ccrs.PlateCarree(), zorder= 20)

    else:
        print(" ".join(['no annotations for \u2220',which_angle, '=', str(val),\
                        u'\N{DEGREE SIGN}']))
    return

# =============================================================================
# mark the coordinates of the points which have specific angle values
# =============================================================================
def annotate_altitude_plots(ax, x, y, ind, val, which_angle, symbol, color, s):
    """
    Marks specific values of angles over plots without map. Choose this routine
    for annotation if x coordinate is not lon. and y coordinate is not lat.

    Parameters
    ----------
    ax: axes to use
    x, y = coordinates to mark: float, array
    ind: index values of close points: float, array
    val: the specific value we are searching for in the data

    Returns
    -------
    None
    """

    if np.size(ind)!=0:
        row, col = np.shape(ind)

        # if there is more than one close value
        if col>1:
            diff_ind= np.diff(ind)
            if (diff_ind[row-1][0]==1):
                ix=ind[row-1][0]
                ax.scatter(x[ix], y[ix], marker = symbol, color= color,\
                        alpha=1, s=s, zorder= 20)
            if (diff_ind[row-1][-1]==1):
                ix=ind[row-1][-1]
                ax.scatter(x[ix], y[ix], marker = symbol, color= color,\
                        alpha=1, s=s, zorder= 20)
            else:
                for i in range(1,len(diff_ind)-1):
                    if (diff_ind[row-1][i]==1 and diff_ind[row-1][i+1]!=1):
                        ix=ind[row-1][i+1]
                        ax.scatter(x[ix], y[ix], marker = symbol, color= color,\
                                  alpha=1, s=s, zorder= 20)
                    elif (diff_ind[row-1][i+1]==1 and diff_ind[row-1][i]!=1):
                        ix=ind[row-1][i+2]
                        ax.scatter(x[ix], y[ix], marker = symbol, color= color,\
                                 alpha=1, s=s, zorder= 20)

        # if there is only one close value
        else:
            ix=ind[row-1][0]
            ax.scatter(x[ix], y[ix], marker =  symbol, color= color,\
              alpha=1, s=s, zorder= 20)

    else:
        print(" ".join(['no annotations for \u2220',which_angle, '=', str(val),\
                        u'\N{DEGREE SIGN}']))

    return
# =============================================================================
# mark the coordinates of the points which have specific angle values in 3D plots
# =============================================================================
def annotate_index_3D(ax, x, y, z, ind, val, which_angle, symbol, color, s):
    """
    Annotations on 3D plots.

    Parameters
    ----------
    ax: axes to use
    x, y = coordinates to mark: float, array
    ind: index values of close points: float, array
    val: the specific value we are searching for in the data

    Returns
    -------
    None
    """
    if np.size(ind)!=0:
        row, col = np.shape(ind)

        # if there is more than one close value
        if col>1:
            diff_ind= np.diff(ind)
            if (diff_ind[row-1][0]==1):
                ix=ind[row-1][0]
                ax.scatter3D(x[ix], y[ix], z[ix], marker = symbol, color= color,\
                        alpha=1, s=s, zorder= 20)
            if (diff_ind[row-1][-1]==1):
                ix=ind[row-1][-1]
                ax.scatter3D(x[ix], y[ix], z[ix], marker = symbol, color= color,\
                        alpha=1, s=s, zorder= 20)
            else:
                for i in range(1,len(diff_ind)-1):
                    if (diff_ind[row-1][i]==1 and diff_ind[row-1][i+1]!=1):
                        ix=ind[row-1][i+1]
                        ax.scatter3D(x[ix], y[ix], z[ix], marker = symbol, color= color,\
                                  alpha=1, s=s, zorder= 20)
                    elif (diff_ind[row-1][i+1]==1 and diff_ind[row-1][i]!=1):
                        ix=ind[row-1][i+2]
                        ax.scatter3D(x[ix], y[ix], z[ix], marker = symbol, color= color,\
                                 alpha=1, s=s, zorder= 20)

        # if there is only one close value
        else:
            ix=ind[row-1][0]
            ax.scatter3D(x[ix], y[ix], z[ix], marker =  symbol, color= color,\
              alpha=1, s=s, zorder= 20)

    else:
        print(" ".join(['no annotations for \u2220',which_angle, '=', str(val),\
                        u'\N{DEGREE SIGN}']))

    return


# =============================================================================
# assign background colors for line plots of angle values
# =============================================================================
def determine_background_colors(angle, ticklimits, colors_for_mapping):

    no_int = len(angle)
    background_color=[[] for i in range(no_int)]

    for i in range(0, no_int):
        if ticklimits[0] <= angle[i] < ticklimits[1]:
            background_color[i]=colors_for_mapping[0]
        elif ticklimits[1]<= angle[i] < ticklimits[2]:
            background_color[i]=colors_for_mapping[1]
        elif ticklimits[2]<= angle[i] < ticklimits[3]:
            background_color[i]=colors_for_mapping[2]
        elif ticklimits[3]<= angle[i] < ticklimits[4]:
            background_color[i]=colors_for_mapping[3]
        elif ticklimits[4]<= angle[i] < ticklimits[5]:
            background_color[i]=colors_for_mapping[4]
        elif ticklimits[5]<= angle[i] < ticklimits[6]:
            background_color[i]=colors_for_mapping[5]
        elif ticklimits[6]<= angle[i] < ticklimits[7]:
            background_color[i]=colors_for_mapping[6]
        elif ticklimits[7]<= angle[i] < ticklimits[8]:
            background_color[i]=colors_for_mapping[7]
        elif ticklimits[8]<= angle[i] < ticklimits[9]:
            background_color[i]=colors_for_mapping[8]
        elif ticklimits[9]<= angle[i] < ticklimits[10]:
            background_color[i]=colors_for_mapping[9]
        elif angle[i] == ticklimits[10]:
             background_color[i]=colors_for_mapping[10]
        elif np.isnan(angle[i])==True:
            background_color[i]='None'
        else:
            background_color[i]='red'

    return background_color


import cartopy.crs as ccrs


def arrow_plotter(ax, Alt, Lat, Lon):
    """
    Plots the arrows to depict the type of trajectory: Ascending /Descending .
    """
    L=np.size(Alt)
    # calculate the differences to determine where the spacecraft is heading to
    # diff_alt  > 0 = descending in altitude
    # diff_alt  < 0 = ascending in altitude
    diff_alt= Alt[0]-Alt[5]
    # diff_lon  > 0 = heading west (as we are westward from Greenwhich)
    # diff_lon  < 0 = heading east
    diff_lon= Lon[0]-Lon[5]
    # diff_lat  > 0 = heading south (as we are northward from equator)
    # diff_lat  < 0 = heading north
    diff_lat= Lat[0]-Lat[5]
    # determine lat_min and max to put the arrows to the start of the passage
    Lat_min, Lat_max=np.min(Lat), np.max(Lat)
    # specify arrrow length and head size using constants
    # c1: arrow length to denote the altitudinal direction
    # c2: arrow length to denote the lat-long plane direction
    # c3: head size for altitudinal direction
    # c4: head size for lat-long plane direction
    c1=1 ; c2= -1.5; c3=0.3; c4=0.3;

    # arrows to denote the changes in altitude, , color: blue
    if (diff_alt>0 and diff_lon>0 and diff_lat>0):
        # ax.arrow(Lon[0],Lat_max, 0, -c1,
        # head_width=c3, head_length=c3, fc='b', ec='b', zorder=25)
        print("Descending in altitude-Going towards: SW")
    elif (diff_alt>0 and diff_lon<0 and diff_lat>0):
        # ax.arrow(Lon[0],Lat_max, 0, -c1,
        # head_width=c3, head_length=c3, fc='b', ec='b', zorder=25)
        print("Descending in altitude-Going towards: SE")
    elif (diff_alt>0 and diff_lon<0 and diff_lat<0):
        # ax.arrow(Lon[L-1],Lat_max, 0, c1,
        # head_width=c3, head_length=c3, fc='b', ec='b', zorder=25)
        print("Descending in altitude-Going towards: NE")
    elif (diff_alt>0 and diff_lon>0 and diff_lat<0):
        # ax.arrow(Lon[L-1],Lat_max, 0, c1,
        # head_width=c3, head_length=c3, fc='b', ec='b', zorder=25)
        print("Descending in altitude-Going towards: NW")
    elif (diff_alt<0 and diff_lon>0 and diff_lat>0):
        # ax.arrow(Lon[0],Lat_max, 0, -c1,
        #      head_width=c3, head_length=c3, fc='b', ec='b', zorder=25)
        print("Ascending in altitude-Going towards: SW")
    elif (diff_alt<0 and diff_lon<0 and diff_lat>0):
        # ax.arrow(Lon[0],Lat_max, 0, -c1,
        # head_width=c3, head_length=c3, fc='b', ec='b', zorder=25)
        print("Ascending in altitude-Going towards: SE")
    elif (diff_alt<0 and diff_lon<0 and diff_lat<0):
        # ax.arrow(Lon[L-1],Lat_max, 0, c1,
        # head_width=c3, head_length=c3, fc='b', ec='b', zorder=25)
        print("Ascending in altitude-Going towards: NE")
    elif (diff_alt<0 and diff_lon>0 and diff_lat<0):
        # ax.arrow(Lon[L-1],Lat_max, 0, c1,
        # head_width=c3, head_length=c3, fc='b', ec='b', zorder=25)
        print("Ascending in altitude-Going towards: NW")

    # arrows for denoting the direction, color: black

    if (diff_alt>0 and diff_lon>0 and diff_lat>0):
        ax.arrow(Lon[0],Lat[0]-.2, diff_lon*c2, diff_lat*c2,
             head_width=c4, head_length=c4, fc='k', ec='k',
             transform=ccrs.PlateCarree(), zorder=25)

    elif (diff_alt>0 and diff_lon<0 and diff_lat>0):
        ax.arrow(Lon[0],Lat[0]-.2, diff_lon*c2, diff_lat*c2,
        head_width=c4, head_length=c4, fc='k', ec='k',
        transform=ccrs.PlateCarree(), zorder=25)

    elif (diff_alt>0 and diff_lon<0 and diff_lat<0):
        ax.arrow(Lon[0],Lat[0]-.2, diff_lon*c2, diff_lat*c2,
        head_width=c4, head_length=c4, fc='k', ec='k',
        transform=ccrs.PlateCarree(), zorder=25)

    elif (diff_alt>0 and diff_lon>0 and diff_lat<0):
        ax.arrow(Lon[0],Lat[0]-.2, diff_lon*c2, diff_lat*c2,
        head_width=c4, head_length=c4, fc='k', ec='k',
        transform=ccrs.PlateCarree(), zorder=25)

    elif (diff_alt<0 and diff_lon>0 and diff_lat>0):
        ax.arrow(Lon[0],Lat[0]-.2, diff_lon*c2, diff_lat*c2,
             head_width=c4, head_length=c4, fc='k', ec='k',
             transform=ccrs.PlateCarree(), zorder=25)

    elif (diff_alt<0 and diff_lon<0 and diff_lat>0):
        ax.arrow(Lon[0],Lat[0]-.2, diff_lon*c2, diff_lat*c2,
        head_width=c4, head_length=c4, fc='k', ec='k',
        transform=ccrs.PlateCarree(), zorder=25)

    elif (diff_alt<0 and diff_lon<0 and diff_lat<0):
        ax.arrow(Lon[0],Lat[0]-.2, diff_lon*c2, diff_lat*c2,
        head_width=c4, head_length=c4, fc='k', ec='k',
        transform=ccrs.PlateCarree(), zorder=25)

    elif (diff_alt<0 and diff_lon>0 and diff_lat<0):
        ax.arrow(Lon[0],Lat[0]+.1, diff_lon*c2, diff_lat*c2,
        head_width=c4, head_length=c4, fc='k', ec='k',
        transform=ccrs.PlateCarree(), zorder=25)


# # =============================================================================
# # multiple x axis
# # =============================================================================
def multiple_x_axes(ax, start_date, end_date, time_array_1sec, step, \
                    parameter, parameter_name, position, pos, n):

    import matplotlib as plt
    # x = -0.065
    # y = -0.365

    # x and y values for implication plot
    x = -0.1
    y = -0.2

    # # x and y values for summary plot
    # x = -0.065
    # y = -0.35
    # x and y values for fft
    # x = -0.065
    # y = -0.12
    offset = plt.rcParams['xtick.major.pad']
    new_tick_locations = time_array_1sec[1::step]
    new_tick_labels=[];
    for i in range(1,len(parameter),step):
        new_tick_labels.append("{0:.1f}".format(float(parameter[i])))

    twin1 = ax.twiny()
    twin1.set_xlim(ax.get_xlim())
    twin1.yaxis.set_visible(False) # hide the yaxis
    twin1.set_xticks(new_tick_locations)
    twin1.set_xticklabels(new_tick_labels, fontsize = 9)
    # set the position of the second x-axis to bottom
    twin1.xaxis.set_ticks_position('bottom')
    twin1.spines['bottom'].set_position(('outward', position+pos*n))
    twin1.spines["bottom"].set_color("white")
    twin1.tick_params(axis='both', which='both', length=2.5)
    twin1.annotate(parameter_name, xy=(x, y), xytext=(x, -offset*n*4.5), \
                ha='left', va='center', xycoords='axes fraction', \
                    textcoords='offset points', fontsize = 9)

    return

def find_ylimits(*var):
    minimums=[]
    maximums=[]
    for v in var:
        minimums.append(np.nanmin(v))
        maximums.append(np.nanmax(v))

    ymin = np.nanmin(minimums)-10
    ymax = np.nanmax(maximums)+10

    return ymin, ymax

def arrange_axes(ax, ax_min, ax_max):
    from matplotlib.ticker import MultipleLocator

    diff = ax_max-ax_min

    if 0 < diff <= 25:
        ax.yaxis.set_major_locator(MultipleLocator(5))
        ax.yaxis.set_minor_locator(MultipleLocator(2.5))

    elif 25 < diff <= 50:
        ax.yaxis.set_major_locator(MultipleLocator(10))
        ax.yaxis.set_minor_locator(MultipleLocator(2.5))
    elif 50 < diff < 90:
        ax.yaxis.set_major_locator(MultipleLocator(20))
        ax.yaxis.set_minor_locator(MultipleLocator(5))
    elif 90 < diff < 180:
        ax.yaxis.set_major_locator(MultipleLocator(45))
        ax.yaxis.set_minor_locator(MultipleLocator(15))

    return

def insert_rectangle_to_background(ax, time1, time2, ydif):
    import matplotlib.patches as patches
    width=time2-time1
    rect = patches.Rectangle((time1, 0), width, ydif,
                           edgecolor='None',
                           facecolor = 'grey', alpha=0.3,
                           fill=True, lw=0)
    ax.add_patch(rect)

    return
