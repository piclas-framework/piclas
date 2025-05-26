#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# System Libaries
import os
import sys
import h5py
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib
matplotlib.rcParams.update({'font.size': 18})
matplotlib.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})
import matplotlib.pyplot as plt
import numpy as np
from progress.bar import IncrementalBar
from termcolor import colored
import skimage


def plotExOverTime():
    RPSet  = 'Line_RPSet.h5'
    ProjectName = 'cvwm'
    # ProjectName = 'turner2013'

    # Extract the RP position
    print('┌'+78*'─')
    print('│ Extracting RP positions')
    with h5py.File(RPSet,'r', driver='core', backing_store=False) as h5:
        FileType = h5.attrs['File_Type']

        if FileType[0].decode('UTF-8') != 'RecordPoints':
            print(' │ File type is not RecordPoints. Exiting...')
            sys.exit()

        # Read in the line coordinates: Get x,y and z for each point: x is not sorted in ascending order
        coords = h5['xF_RP'][:]
        print('│ shape of xF_RP', np.shape(coords))

        # Sort x,y and z coords along the dimension 0, which is the x-coordinate (in descending order)
        arg = coords[:, 0].argsort()
        xF_RP_sorted = coords[arg[::-1]]  # when you do a[::-1], it starts from the end towards the first taking each element. So it reverses a

        # Get only the x-coordinate. It is sorted in descending order
        xF_RP = xF_RP_sorted[:, 0]

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

        FirstFile = RPFiles.pop(0)  # git stash pop
        with h5py.File(FirstFile, 'r', driver='core', backing_store=False) as h5:
            RP_Data_h5 = h5['RP_Data'][:]
        pbar.next()

        for RPnum, RPfile in enumerate(RPFiles):
            with h5py.File(RPfile, 'r', driver='core', backing_store=False) as h5:
                current_data = h5['RP_Data'][:][1:, :, :]
                RP_Data_h5 = np.vstack((RP_Data_h5, current_data))
            pbar.next()

        pbar.finish()

        xF_RP_sorted = coords[arg[::-1], 0]
        Ex_sorted = RP_Data_h5[:, arg[::-1], 2]
        print('│ shape of Ex_sorted', np.shape(Ex_sorted))
        time = RP_Data_h5[:, 0, 0]*1e6
        print('│ shape of RP_Data_h5[:, 0, 0]', np.shape(RP_Data_h5[:, 0, 0]))
        print('│ shape of time', np.shape(time))

        scale=10.0
        # scale=1e4
        vmin = -scale
        vmax = scale
        print('│ shape of RP_Data', np.shape(RP_Data_h5))
        for cutoff in [2.25, 5.00]:
            if cutoff < 5:
                fig, ax = plt.subplots(1, 1, figsize=(4, 5))
                scale=0.85
                fig, ax = plt.subplots(1, 1, figsize=(4*scale, 5*scale))
            else:
                fig, ax = plt.subplots(1, 1, figsize=(4, 11))

            # Apply a fancy colormap to the figure
            cmap = plt.get_cmap('seismic')
            print('│ shape of time[time <= %s]' % cutoff, np.shape(time[time <= cutoff]))
            print('│ shape of Ex_sorted[time <= %s, :]' % cutoff, np.shape(Ex_sorted[time <= cutoff, :]))

            # contourf_ = ax.pcolormesh(xF_RP[:], time[:], Ex_sorted[:, :], vmin=vmin, vmax=vmax, cmap=cmap)
            contourf_ = ax.pcolormesh(xF_RP[:], time[time <= cutoff], Ex_sorted[time <= cutoff, :], vmin=vmin, vmax=vmax, cmap=cmap)

            # create an axes on the right side of ax. The width of cax will be 5% of ax and the padding between cax and ax will be fixed at 0.05 inch.
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            cbar = plt.colorbar(contourf_, cax=cax, label='test')
            cbar.set_label(r'$E_x~\mathrm{[V/m]}$')

            # cbar.set_label('Mach number', rotation=270)
            ax.set_xlabel(r'$x~\mathrm{[m]}$')
            ax.set_ylabel(r'$t~\mathrm{[\mu s]}$')
            ax.set_aspect('auto')  # or 'equal'
            if cutoff < 5:
                plt.savefig('Ex-over-t-cut-off-%se-6.png' % cutoff, bbox_inches='tight', dpi=300)
            else:
                plt.savefig('Ex-over-t.png', bbox_inches='tight', dpi=300)

    print('└'+78*'─')
    return xF_RP, time/1e6, Ex_sorted


def plotFFT(x, t, Ex):
    # Set parameters
    omega_pe  = 3.98911503E07
    lambda_De = 1.051113E-02
    dt        = t[1] - t[0]
    dx        = x[1] - x[0]

    Ex = np.transpose(Ex)
    print('│ shape of Ex', np.shape(Ex))

    # Hann filter in 2D
    win = skimage.filters.window('hann', (len(x), len(t)))

    # Create FFTs
    F_x_omega = np.fft.fft(win*Ex, axis=1)
    F_x_omega = np.fft.fftshift(F_x_omega, axes=1)  # centre frequencies

    F_k_omega = np.fft.fft(F_x_omega, axis=0)
    F_k_omega = np.fft.fftshift(F_k_omega, axes=0)  # centre frequencies

    # Frequenzachsen
    Nt    = len(t)
    Nx    = len(x)
    omega = np.fft.fftshift(np.fft.fftfreq(Nt, dt)) * 2 * np.pi /omega_pe # Umrechnung in rad/s
    k     = np.fft.fftshift(np.fft.fftfreq(Nx, dx)) * 2 * np.pi *lambda_De # Umrechnung in rad/s

    scale = 10.0
    vmin  = 0.0
    vmax  = np.max(np.abs(F_k_omega.T))/scale

    # Plot
    fig, ax = plt.subplots(1, 1, figsize=(5, 4.25))
    plt.pcolormesh(k, omega, np.abs(F_k_omega.T), vmin=vmin, vmax=vmax, cmap=plt.get_cmap('magma'))
    plt.xlabel(r'$k\lambda_{De}$ [-]')
    plt.ylabel(r'$\omega/\omega_{pe}$ [-]')
    ax.set_xlim(-5, 5)
    ax.set_ylim(-5, 5)

    # Plot analytical function
    x = np.arange(-5.0, 5.0, 0.1)
    y = 50*np.sin(x)
    y = np.sqrt(3.0*x**2+1.0)
    plt.plot(x,y, 'w--', lw=1)

    plt.savefig('FFT-omega-over-lambda.png', bbox_inches='tight', dpi=300)
    plt.close('all')

    print('└'+78*'─')


# Program starts here

# Create Ex plot
x, t, Ex = plotExOverTime()

# Create FFT plot
plotFFT(x, t, Ex)
