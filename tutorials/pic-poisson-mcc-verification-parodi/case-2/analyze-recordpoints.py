#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# System Libaries
import os
import sys
import h5py
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import numpy as np
from progress.bar import IncrementalBar
from termcolor import colored
import matplotlib

def plotExOverTime():
    RPSet  = 'Line_RPSet.h5'
    ProjectName = 'cvwm'

    # Extract the RP position
    print('┌'+78*'─')
    print('│ Extracting RP positions')
    with h5py.File(RPSet,'r', driver='core', backing_store=False) as h5:
        FileType = h5.attrs['File_Type']

        if FileType[0].decode('UTF-8') != 'RecordPoints':
            print(' │ File type is not RecordPoints. Exiting...')
            sys.exit()

        # Read in the group names
        GroupNames = h5['GroupNames'][:]
        group = np.where(GroupNames == str.encode('Line'))

        # Read in the line coordinates
        # Get x,y and z for each point: x is not sorted in ascending order
        coords = h5['xF_RP'][:]
        print('│ shape of xF_RP', np.shape(coords))
        # print('\n'+'coords:')
        # print(coords)

        # Sort x,y and z coords along the dimension 0, which is the x-coordinate (in descending order)
        arg = coords[:, 0].argsort()
        # print('\n'+'arg:')
        # print(arg)
        # print('\n'+'arg[::-1]:')
        # print(arg[::-1])
        xF_RP_sorted = coords[arg[::-1]]  # when you do a[::-1], it starts from the end towards the first taking each element. So it reverses a
        # print('\n'+'xF_RP_sorted:')
        # print(xF_RP_sorted)

        # Get only the x-coordinate. It is sorted in descending order
        xF_RP = xF_RP_sorted[:, 0]
        # print('\n'+'xF_RP_sorted[:, 0]:')
        # print(xF_RP)

    filenames = [f for f in os.listdir('.') if os.path.isfile(f)]
    RPFiles = [kk for kk in filenames if '.h5' in kk and '{}_RP_'.format(ProjectName) in kk]
    RPFiles.sort()

    print('├'+78*'─'+'\n'+'│ READING RPDATA...')

    # Load the data
    print('│ Loading data from {} files'.format(len(RPFiles)))

    if os.path.exists('RP_Data.npz'):
        print(colored('├─ Loading numpy data from compressed file "{}"'.format('RP_Data.npz')))
        pbar = IncrementalBar('├─ Loading numpy data        ',max=1,suffix='%(percent)d%% [%(elapsed)ds / %(eta)ds]',check_tty=False)
        pbar.start()
        RP_Data = np.load('RP_Data.npz')['RP_Data']
        pbar.next()
        pbar.finish()

    else:
        # Just a progress bar
        pbar = IncrementalBar('├─ Loading RP files          ',max=len(RPFiles),suffix='%(percent)d%% [%(elapsed)ds / %(eta)ds]',check_tty=False)
        pbar.start()

        # print('\n'+'├'+78*'─')
        FirstFile = RPFiles.pop(0)  # git stash pop
        # print(FirstFile)
        with h5py.File(FirstFile, 'r', driver='core', backing_store=False) as h5:
            RP_Data_h5 = h5['RP_Data'][:]
        pbar.next()

        for RPnum, RPfile in enumerate(RPFiles):
            # print(RPfile)
            with h5py.File(RPfile, 'r', driver='core', backing_store=False) as h5:
                current_data = h5['RP_Data'][:][1:, :, :]
                RP_Data_h5 = np.vstack((RP_Data_h5, current_data))
            pbar.next()

        pbar.finish()
        plt.rcParams.update({'font.size': 18})
        # plt.rcParams['errorbar.capsize'] = 3
        # plt.rcParams['font.size'] = '12'
        # plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
        plt.rcParams.update({
            "text.usetex": True,
            "font.family": "Helvetica"
        })

        vmin = -10.
        vmax = 10.
        print('│ shape of RP_Data', np.shape(RP_Data_h5))
        cutoff=2.25
        # cutoff=5.00
        if cutoff < 5:
            fig, ax = plt.subplots(1, 1, figsize=(4, 5))
            scale=0.85
            fig, ax = plt.subplots(1, 1, figsize=(4*scale, 5*scale))
        else:
            fig, ax = plt.subplots(1, 1, figsize=(4, 11))

        # Apply a fancy colormap to the figure
        cmap = plt.get_cmap('seismic')

        xF_RP_sorted   = coords[arg[::-1], 0]
        Ex_sorted = RP_Data_h5[:, arg[::-1], 2]
        print('│ shape of Ex_sorted', np.shape(Ex_sorted))
        time = RP_Data_h5[:, 0, 0]*1e6
        print('│ shape of RP_Data_h5[:, 0, 0]', np.shape(RP_Data_h5[:, 0, 0]))
        print('│ shape of time', np.shape(time))
        print('│ shape of time[time <= %s]' % cutoff, np.shape(time[time <= cutoff]))
        print('│ shape of Ex_sorted[time <= %s, :]' % cutoff, np.shape(Ex_sorted[time <= cutoff, :]))

        # contourf_ = ax.pcolormesh(xF_RP[:], time[:], Ex_sorted[:, :], vmin=vmin, vmax=vmax, cmap=cmap)
        contourf_ = ax.pcolormesh(xF_RP[:], time[time <= cutoff], Ex_sorted[time <= cutoff, :], vmin=vmin, vmax=vmax, cmap=cmap)

        # create an axes on the right side of ax. The width of cax will be 5%
        # of ax and the padding between cax and ax will be fixed at 0.05 inch.
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(contourf_, cax=cax, label='test')
        cbar.set_label(r'$E_x~\mathrm{[V/m]}$')

        # cbar.set_label('Mach number', rotation=270)
        ax.set_xlabel(r'$x~\mathrm{[m]}$')
        ax.set_ylabel(r'$t~\mathrm{[µs]}$')
        ax.set_aspect('auto')  # or 'equal'
        if cutoff < 5:
            plt.savefig('Ex-over-t-cut-off-%se-6.png' % cutoff, bbox_inches='tight', dpi=300)
        else:
            plt.savefig('Ex-over-t.png', bbox_inches='tight', dpi=300)
    print('└'+78*'─')


# Program starts here
plotExOverTime()
