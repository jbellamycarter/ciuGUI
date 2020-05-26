#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ciuGUI is a simple GUI for comparing two CIU datasets.

Two .csv files with the standard CIU format are imported--these must
contain the same number of rows and columns. The data can then be
smoothed with either 1D area or intensity normalisation.

The data are displayed as a four-pane figure. The top two panes show the
heatmaps for each data set in the selected colormap. The bottom-left
pane shows the RMSD heatmap and the rolling RMSD (across the collision
energies). The bottom-right pane shows the overlaid mobiligrams for a
selected energy, the corresponding energy is also displayed as a line on
the top panes.

This is a free and open-source software for the specific task of
comparing CIU datasets. Other, more comprehensive, software is available
such as CIUSuite, ORIGAMI and PULSAR.

Copyright (c) 2016-2020 Jedd Bellamy-Carter, University of Nottingham
Released under a GNU Public License

Created on Thu Dec  8 23:15:32 2016
"""

#    Copyright 2016 Jedd Bellamy-Carter, University of Nottingham

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.

import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator
import numpy as np
import tkinter as tk
import tkinter.filedialog as filedialog
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backend_bases import key_press_handler

mpl.rcParams['font.size'] = 9


class appGui(tk.Frame):
    def __init__(self, master):
        tk.Frame.__init__(self, master)

        ## Define main figure with gridspec
        self.fig = plt.figure(frameon=False, dpi=120)
        self.gs = gridspec.GridSpec(2, 2)

        ## First CIU Profile
        self.ciu1 = plt.subplot(self.gs[0, 0])
        self.ciu1ax = self.ciu1.imshow([[0, 0], [0, 0]],
                                       aspect='auto',
                                       cmap='jet',
                                       origin='lower',
                                       interpolation='gaussian',
                                       extent=[0, 1, 0, 1],
                                       vmin=0,
                                       vmax=1)
        self.ciu1.yaxis.set_major_locator(MaxNLocator(integer=True))
        self.ciu1text = self.ciu1.text(0.95,
                                       0.9,
                                       '',
                                       color='w',
                                       horizontalalignment='right',
                                       transform=self.ciu1.transAxes)
        self.ciu1v = self.ciu1.axvline(0, c='w', ls=':', lw=1)
        self.ciu1.set_ylabel('Arrival Time (ms)')
        self.ciu1.set_xlabel('Collision Energy (eV)')

        ## Second CIU Profile
        self.ciu2 = plt.subplot(self.gs[0, 1])
        self.ciu2ax = self.ciu2.imshow([[0, 0], [0, 0]],
                                       aspect='auto',
                                       cmap='jet',
                                       origin='lower',
                                       interpolation='gaussian',
                                       extent=[0, 1, 0, 1],
                                       vmin=0,
                                       vmax=1)
        self.ciu2.yaxis.set_major_locator(MaxNLocator(integer=True))
        self.ciu2text = self.ciu2.text(0.95,
                                       0.9,
                                       '',
                                       color='w',
                                       horizontalalignment='right',
                                       transform=self.ciu2.transAxes)
        self.ciu2v = self.ciu2.axvline(0, c='w', ls=':', lw=1)
        self.ciu2.set_ylabel('Arrival Time (ms)')
        self.ciu2.set_xlabel('Collision Energy (eV)')

        self.gs1 = gridspec.GridSpecFromSubplotSpec(3,
                                                    2,
                                                    subplot_spec=self.gs[2])

        ## Diff Plot
        self.diffFig = plt.subplot(self.gs1[:2, :2])
        self.diffax = self.diffFig.imshow([[0, 0], [0, 0]],
                                          aspect='auto',
                                          cmap='RdBu',
                                          origin='lower',
                                          interpolation='gaussian',
                                          extent=[0, 1, 0, 1],
                                          vmin=-1,
                                          vmax=1)
        self.difftext = self.diffFig.text(0.95,
                                          0.85,
                                          '',
                                          horizontalalignment='right',
                                          transform=self.diffFig.transAxes)
        self.diffFig.set_ylabel('Arrival Time (ms)')
        plt.setp(self.diffFig.get_xticklabels(), visible=False)

        ## RMSD
        self.rmsdFig = plt.subplot(self.gs1[2, :2], sharex=self.diffFig)
        self.rmsdax = self.rmsdFig.plot([0, 1], [0, 0], 'k.')
        self.rmsdhln = self.rmsdFig.axhline(0, c='k', ls=':')
        self.rmsdFig.set_ylim(0)
        self.rmsdFig.locator_params(axis='y', nbins=3)
        self.rmsdFig.set_ylabel('RMSD (%)')
        self.rmsdFig.set_xlabel('Collision Energy (eV)')

        ## Slice
        self.slice = plt.subplot(self.gs[1, 1])
        self.ln1 = self.slice.plot([], [], '-', c='#0571b0', label='')
        self.ln2 = self.slice.plot([], [],
                                   '--',
                                   dashes=(2, 2),
                                   c='#ca0020',
                                   label='')
        self.slicetext = self.slice.text(0.95,
                                         0.1,
                                         '',
                                         horizontalalignment='right',
                                         transform=self.slice.transAxes)
        self.slice.set_ylabel('Relative Intensity (%)')
        self.slice.set_xlabel('Arrival Time (ms)')
        self.slice.set_yticks([0, .2, .4, .6, .8, 1])
        self.slice.set_yticklabels([0, 20, 40, 60, 80, 100])
        self.slice.set_ylim(0, 1)

        self.gs.tight_layout(self.fig, h_pad=1)

        self.master = master
        self.initUI()

    def initUI(self):
        self.master.title("CIU Visualiser")

        ## Parameters
        tk.Button(self.master,
                  text="Open CIU Data",
                  command=self.askopenfilename).grid(row=0,
                                                     column=0,
                                                     sticky="ew",
                                                     padx=5)
        tk.Button(self.master, text="Save CIU Data",
                  command=self.saveCIU).grid(row=0,
                                             column=1,
                                             sticky="ew",
                                             padx=5)

        params1 = tk.LabelFrame(self.master, text="First Species")
        params1.grid(row=2, column=0, columnspan=2, sticky="ew", padx=5)
        params1.grid_columnconfigure(0, weight=1)
        tk.Label(self.master, text="File").grid(in_=params1,
                                                row=0,
                                                column=0,
                                                sticky="w")
        tk.Label(self.master, textvariable=file1).grid(in_=params1,
                                                       row=0,
                                                       column=1,
                                                       sticky="ew")
        tk.Label(self.master, text="Name").grid(in_=params1,
                                                row=1,
                                                column=0,
                                                sticky="w")
        tk.Entry(self.master, textvariable=name1).grid(in_=params1,
                                                       row=1,
                                                       column=1,
                                                       sticky="news")
        tk.Label(self.master, text="Charge State").grid(in_=params1,
                                                        row=2,
                                                        column=0,
                                                        sticky="w")
        chg1Ent = tk.Entry(self.master, textvariable=chg1,
                           width=5).grid(in_=params1,
                                         row=2,
                                         column=1,
                                         sticky="ew")
        self.master.bind('<Return>', self.update)
        self.master.bind('<Control_L>o', self.askopenfilename)
        self.master.bind('<Control_L>s', self.saveCIU)
        self.master.bind('<Control_L>e', self.exportImage)

        params2 = tk.LabelFrame(self.master, text="Second Species")
        params2.grid(row=4, column=0, columnspan=2, sticky="ew", padx=5)
        params2.grid_columnconfigure(0, weight=1)
        tk.Label(self.master, text="File").grid(in_=params2,
                                                row=0,
                                                column=0,
                                                sticky="w")
        tk.Label(self.master, textvariable=file2).grid(in_=params2,
                                                       row=0,
                                                       column=1,
                                                       sticky="ew")
        tk.Label(self.master, text="Name").grid(in_=params2,
                                                row=1,
                                                column=0,
                                                sticky="w")
        tk.Entry(self.master, textvariable=name2).grid(in_=params2,
                                                       row=1,
                                                       column=1,
                                                       sticky="news")
        tk.Label(self.master, text="Charge State").grid(in_=params2,
                                                        row=2,
                                                        column=0,
                                                        sticky="w")
        tk.Entry(self.master, textvariable=chg2, width=5).grid(in_=params2,
                                                               row=2,
                                                               column=1,
                                                               sticky="ew")

        displayParams = tk.LabelFrame(self.master, text="Parameters")
        displayParams.grid(row=7, column=0, columnspan=2, sticky="ew", padx=5)
        displayParams.grid_columnconfigure(0, weight=1)
        displayParams.grid_columnconfigure(1, weight=1)
        tk.Label(self.master, text="Main Colormap").grid(in_=displayParams,
                                                         row=0,
                                                         column=0,
                                                         sticky="w")
        tk.Entry(self.master, textvariable=colorMap,
                 width=5).grid(in_=displayParams, row=0, column=1, sticky="ew")
        tk.Label(self.master, text="Smooth Window (±)").grid(in_=displayParams,
                                                             row=1,
                                                             column=0,
                                                             sticky="w")
        tk.Entry(self.master, textvariable=smoothWin,
                 width=5).grid(in_=displayParams, row=1, column=1, sticky="ew")
        tk.Label(self.master, text="Number of Smooths").grid(in_=displayParams,
                                                             row=2,
                                                             column=0,
                                                             sticky="w")
        tk.Entry(self.master, textvariable=smoothNum,
                 width=5).grid(in_=displayParams, row=2, column=1, sticky="ew")
        tk.Label(self.master, text="Normalisation").grid(in_=displayParams,
                                                         row=3,
                                                         column=0,
                                                         rowspan=2,
                                                         sticky="w")
        tk.Radiobutton(self.master, text="Area", variable=norm,
                       value=0).grid(in_=displayParams,
                                     row=3,
                                     column=1,
                                     sticky="w")
        tk.Radiobutton(self.master, text="Intensity", variable=norm,
                       value=1).grid(in_=displayParams,
                                     row=4,
                                     column=1,
                                     sticky="w")
        tk.Label(self.master, text="Collision Energy").grid(in_=displayParams,
                                                            row=5,
                                                            column=0,
                                                            sticky="w")
        self.colMenu = tk.OptionMenu(self.master, colEnergy, *colEnergies)
        self.colMenu.grid(in_=displayParams, row=5, column=1, sticky="ew")
        tk.Label(self.master).grid(row=9, column=0, sticky="news")

        tk.Button(self.master,
                  text="Update Graphs",
                  command=lambda: self.update('dump')).grid(row=8,
                                                            column=0,
                                                            sticky="ew",
                                                            padx=5)
        tk.Button(self.master, text="Export Image",
                  command=self.exportImage).grid(row=8,
                                                 column=1,
                                                 sticky="ew",
                                                 padx=5)

        figFrame = tk.Frame(self.master)
        figFrame.grid(row=0, column=2, rowspan=9, sticky="news")
        figFrame.grid_columnconfigure(0, weight=1)
        figFrame.grid_rowconfigure(0, weight=1)
        canvas = FigureCanvasTkAgg(self.fig, self.master)
        #canvas.show()
        self.fig.canvas.draw()
        canvas.get_tk_widget().grid(in_=figFrame,
                                    row=0,
                                    column=0,
                                    sticky="news")

    def askopenfilename(self,dump):
        global drift1
        global energy1
        global data1
        global drift2
        global energy2
        global data2
        global filename
        filename = filedialog.askopenfilenames(
            title='Choose CSV containing CIU data',
            filetypes=[
                ('CSV file', '.csv'),
                ('All files', '.*'),])
        if len(filename) == 1:
            drift1, energy1, data1 = self.readCIU(filename[0])
            drift2, energy2, data2 = drift1, energy1, np.zeros_like(data1)
            file1.set(os.path.basename(filename[0]))
            file2.set('None')
        elif len(filename) >= 2:
            drift1, energy1, data1 = self.readCIU(filename[0])
            drift2, energy2, data2 = self.readCIU(filename[1])
            file1.set(os.path.basename(filename[0]))
            file2.set(os.path.basename(filename[1]))

        self.update('dump')

    def exportImage(self, dump):
        filename = filedialog.asksaveasfilename(title='Export image',
                                                filetypes=[
                                                    ('PNG file', '.png'),
                                                    ('All files', '.*'),])
        self.fig.savefig(filename, transparent=True, dpi=600)

    def update(self, dump):
        colEnergies = energy1 * chg1.get()
        colEnergies2 = energy2 * chg2.get()
        if colEnergy.get() not in colEnergies: colEnergy.set(colEnergies[0])
        ## Update data
        normData1 = self.normData(data1, norm.get())
        normData2 = self.normData(data2, norm.get())
        diffData = normData1 - normData2
        rmsdData = self.rollingRMSD(diffData)
        self.ciu1ax.set_cmap(colorMap.get())
        self.ciu2ax.set_cmap(colorMap.get())
        self.ciu1ax.set_data(normData1)
        self.ciu2ax.set_data(normData2)
        self.ciu1ax.set_extent(
            [colEnergies[0], colEnergies[-1], drift1[0], drift1[-1]])
        self.ciu2ax.set_extent(
            [colEnergies2[0], colEnergies2[-1], drift1[0], drift1[-1]])
        self.ciu1text.set_text(r'%s $^{%d{+}}$' % (name1.get(), chg1.get()))
        self.ciu2text.set_text(r'%s $^{%d{+}}$' % (name2.get(), chg2.get()))
        self.slice.legend(loc=1,
                          fontsize=9,
                          frameon=False,
                          labels=[name1.get(), name2.get()])
        self.diffax.set_extent(
            [colEnergies[0], colEnergies[-1], drift1[0], drift1[-1]])
        self.diffax.set_data(diffData)
        self.difftext.set_text("RMSD = %.2f%%" % rmsdData.mean())
        self.rmsdax[0].set_data(colEnergies, rmsdData)
        self.rmsdFig.set_xlim(colEnergies[0], colEnergies[-1])
        self.rmsdFig.set_ylim(0, np.ceil(rmsdData.max()))
        self.rmsdhln.set_ydata([rmsdData.mean()])

        tarIdx1 = np.argwhere(colEnergies == colEnergy.get())
        tarIdx2 = np.argwhere(colEnergies2 == colEnergy.get())
        if tarIdx1 >= 0 and tarIdx2 >= 0:
            datMax = np.max(
                [normData1[:, tarIdx1].max(), normData2[:, tarIdx2].max()])
            self.ln1[0].set_data(drift1, normData1[:, tarIdx1] / datMax)
            self.ciu1v.set_xdata(colEnergies[tarIdx1])
            self.slicetext.set_text("%.1f eV" % colEnergy.get())
            self.ln2[0].set_data(drift2, normData2[:, tarIdx2] / datMax)
            self.ciu2v.set_xdata(colEnergies2[tarIdx2])
            self.slice.set_xlim(drift1[0], drift1[-1])

        self.fig.canvas.draw()

        ## Reset colMenu
        menu = self.colMenu['menu']
        menu.delete(0, 'end')
        for val in colEnergies:
            menu.add_command(label=val, command=lambda v=val: colEnergy.set(v))

    def _quit(self):
        root.quit()  # stops mainloop
        root.destroy()

    def readCIU(self, _file):
        raw = np.genfromtxt(_file, delimiter=',')
        return raw[1:, 0], raw[0, 1:], raw[1:, 1:]  ## Drifts, energies, data

    def saveCIU(self, dump):
        name1, ext = os.path.splitext(filename[0])
        raw = np.genfromtxt(filename[0], delimiter=',')
        raw[1:, 1:] = self.normData(data1, norm.get())
        np.savetxt(name1 + '_norm.csv', raw, '%.5f', delimiter=',')
        name2, ext = os.path.splitext(filename[1])
        raw = np.genfromtxt(filename[1], delimiter=',')
        raw[1:, 1:] = self.normData(data2, norm.get())
        np.savetxt(name2 + '_norm.csv', raw, '%.5f', delimiter=',')

    def meanSmooth(self, _data, _window, _smooths):
        for smooth in range(_smooths):
            _data = np.convolve(_data,
                                np.ones((_window, )) / _window,
                                mode='same')
        return _data

    def smoothData(self, _tic):
        smoothTic = np.zeros_like(_tic.T, dtype='float64')
        for i in range(len(_tic.T)):
            smoothTic[i] = self.meanSmooth(_tic.T[i],
                                           (smoothWin.get() * 2) + 1,
                                           smoothNum.get())
        return smoothTic

    def normData(self, _tic, _type=0):
        smoothTic = self.smoothData(_tic)
        if _type == 1:
            return smoothTic.T / smoothTic.max(axis=1)
        elif _type == 0:
            area = smoothTic.T / np.trapz(smoothTic, axis=1)
            areaNorm = area / area.sum(axis=0)
            return areaNorm / areaNorm.max()

    def RMSD(self, _diff):
        return np.sqrt(np.sum(np.power(_diff, 2.)) / _diff.size) * 100

    def rollingRMSD(self, _diff):
        _rmsd = np.zeros(len(_diff[0]), dtype='float64')
        for i in range(len(_diff[0])):
            _rmsd[i] = self.RMSD(_diff[:, i])
        return _rmsd


if __name__ == '__main__':
    root = tk.Tk()
    root.option_add("*Font", "Arial, 11")
    root.geometry("1200x800+0+0")
    root.columnconfigure(2, weight=1)
    root.rowconfigure(0, weight=1)
    root.rowconfigure(1, weight=1)
    root.rowconfigure(2, weight=1)
    root.rowconfigure(3, weight=1)
    root.rowconfigure(4, weight=1)
    root.rowconfigure(5, weight=1)
    root.rowconfigure(6, weight=1)
    root.rowconfigure(7, weight=1)
    root.rowconfigure(8, weight=1)
    root.rowconfigure(9, weight=5)
    file1 = tk.StringVar()
    file2 = tk.StringVar()
    colorMap = tk.StringVar()
    colorMap.set("jet")
    smoothWin = tk.IntVar()
    smoothWin.set(3)
    smoothNum = tk.IntVar()
    smoothNum.set(2)
    name1 = tk.StringVar()
    chg1 = tk.IntVar()
    chg1.set(1)
    name2 = tk.StringVar()
    chg2 = tk.IntVar()
    chg2.set(1)
    norm = tk.IntVar()
    norm.set(1)
    colEnergies = [0.]
    colEnergy = tk.DoubleVar()
    exFig = tk.IntVar()
    app = appGui(root)
    root.mainloop()
