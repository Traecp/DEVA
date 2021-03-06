#!/usr/bin/env python
# -*- coding: utf-8 -*-
###############################################################################
## Filename:          libXpad.py
## Description:       Xpad detector library
## Author:            Clement Buton <clement.buton@esrf.fr>
## Author:            $Author: cbuton $
## Created on:        $Date: 2013/02/12 17:20:56 $
## Modified on:       2013/09/18 09:59:38
## Copyright:         2013, Clément Buton (CeCILL-C)
## $Id: libXpad.py, 2013/02/12 17:20:56 cbuton Exp $
###############################################################################

"""
.. _libXpad:

Calibration.libCalibration - Xpad detector library
==================================================

.. todo::
    * Add widgets (sliders, buttons, radio) to plots

.. note::
    * Chips per module: 7
    * Modules number: {D1:8, D5:8, imXpad2:2, S70:1}
    * Pixel size: 0.13 x 0.13 mm²
    * Chip's rows: 120
    * Chip's columns: 80
    * Magnification factor for the pixels at the chip sides: 32/13.
    * Magnification factor for the pixels at the module's edges: 1.2
    * Module periodicity: {D5:19.17} mm
"""

__author__ = 'Clement Buton <clement.buton@esrf.fr>'
__date__ = '2013/02/12 17:20:56'
__adv__ = 'libXpad.py'

import os
import numpy as N

import DEVA.utilities.ColoredText as CT
from DEVA.utilities.Statistics import median_stats

from matplotlib.cm import jet
from collections import namedtuple

# Classes =====================================================================

class Chip(object):
    """Elementary brick of a detector."""

    def __init__(self,
                 chipShape=(120, 80),
                 location='core',
                 name=u'C_#',
                 pixelSize=(0.13, 0.13),
                 edgeFactor=1.2,
                 sideFactor=32/13.):
        """Class 'Chip' initialization.

        :param chipShape: Chip shape (rows,columns) [pixels]
        :param location: Physical location in the module (single, 'core', 'left' or 'right')
        :param name: Nomenclature of the chip (i.e. 'C_4')
        :param pixelSize: Pixels size (y,x) [mm]
        :param edgeFactor: Magnification factor for the pixels at the module's edges (y,x)
        :param sideFactor: Magnification factor for the wider pixels at the chip sides (y)
        """

        # Inputs
        self.shape = chipShape     # [pixels]
        self.name = name

        # Initialize the config file header for a single chip
        self.dacG = OrderedDict([('CMOS_DSBL', 0),
                                 ('AMP_TP', 0),
                                 ('ITHH', 0),
                                 ('VADJ', 0),
                                 ('VREF', 0),
                                 ('IMFP', 52),
                                 ('IOTA', 40),
                                 ('IPRE', 60),
                                 ('ITHL', 25),
                                 ('ITUNE', 120),
                                 ('IBUFFER', 0)])

        # Data initialization
        self.data = N.zeros(self.shape)

        # Thresholds initialization
        self.ith = self.dacG['ITHL']         # One per chip
        self.dacl = N.ones(self.shape) * 32  # One per pixel

        # Number of pixels in the chip
        self.nPixels = self.shape[0] * self.shape[1]

        # Initialize the detector physical description
        self.physical = namedtuple('physical', 'location pixelSize edgeFactor \
                                   sideFactor size shape data dacl')

        self.physical.location = location
        self.physical.pixelSize = pixelSize
        self.physical.edgeFactor = edgeFactor
        self.physical.sideFactor = sideFactor

        # Pysical chip size
        if self.physical.location == 'single':
            size = (chipShape[0] * self.physical.pixelSize[0],
                    chipShape[1] * self.physical.pixelSize[1])

        elif self.physical.location == 'core':
            size = (chipShape[0] * self.physical.pixelSize[0],
                    (chipShape[1] - 2) * self.physical.pixelSize[1] + \
                    (2 * self.physical.pixelSize[1] * self.physical.sideFactor))

        else:
            size = (chipShape[0] * self.physical.pixelSize[0],
                    (chipShape[1] - 2) * self.physical.pixelSize[1] + \
                    (self.physical.pixelSize[1] * self.physical.sideFactor) + \
                    (self.physical.pixelSize[1] * self.physical.edgeFactor))

        self.physical.size = tuple(N.round(size, 2)) # physical size [mm x mm]

        nRows = N.round(self.physical.size[0] / self.physical.pixelSize[0])
        nCols = N.round(self.physical.size[1] / self.physical.pixelSize[1])

        self.physical.shape = (nRows, nCols)
        self.physical.data = N.zeros(self.physical.shape)
        self.physical.dacl = N.zeros(self.physical.shape)

    def __str__(self):
        """Print out Class information.

        :return: return Chip information
        """

        s = "Chip '%s': (%i,%i) pixels" % \
            (self.name, self.shape[0], self.shape[1])

        return s

    def update_dacG(self, keyword, value):
        """Update the the dacG values of the chip.

        :param keyword: Configuration Keyword
        :param value: Corresponding keyword value
        """

        value = int(value)
        setattr(self, keyword.lower(), value)
        self.dacG[keyword] = value  # update a given dacG key

    def set_data(self, array, dacl=False):
        """Set the data array (or the dacl array) of the chip.

        :param array: data array (chip shape)
        :param dacl: if True, set the dacl array of the chip
        """

        if isinstance(array, N.ndarray):
            assert array.shape == self.shape, \
                "The input array should have the following shape: (%i, %i) !" \
                % (self.shape[0], self.shape[1])

            if dacl:
                self.dacl = array

            else:
                self.data = array

        else:
            raise IOError('The input array should be a 2D-numpy.ndarray !')

    def adjust_data(self, log=False):
        """Adjust the intensity of the 320 µm wide pixels (i.e. divided
        it by 32/13.) and the module edge pixels.

        :param log: data in log
        """

        if self.physical.location == 'single':
            testFactor = 2.
            if log is False:
                self.data[0][1:-1] /= testFactor      # Chip's first row
                self.data[-1][1:-1] /= testFactor     # Chip's last row
                self.data[:, 0] /= testFactor         # Chip's first column
                self.data[:, -1] /= testFactor        # Chip's last column

            else:
                self.data[0][1:-1] -= N.log10(testFactor)      # Chip's first row
                self.data[-1][1:-1] -= N.log10(testFactor)     # Chip's last row
                self.data[:, 0] -= N.log10(testFactor)         # Chip's first column
                self.data[:, -1] -= N.log10(testFactor)        # Chip's last column

        elif self.physical.location == 'core':
            if log is False:
                self.data[0][1:-1] /= self.physical.edgeFactor     # Chip's first row
                self.data[-1][1:-1] /= self.physical.edgeFactor    # Chip's last row
                self.data[:, 0] /= self.physical.sideFactor  # Chip's first column
                self.data[:, -1] /= self.physical.sideFactor # Chip's last column

            else:
                self.data[0][1:-1] -= N.log10(self.physical.edgeFactor)  # Chip's first row
                self.data[-1][1:-1] -= N.log10(self.physical.edgeFactor) # Chip's last row
                self.data[:, 0] -= N.log10(self.physical.sideFactor)     # Chip's first column
                self.data[:, -1] -= N.log10(self.physical.sideFactor)    # Chip's last column

        elif self.physical.location == 'left':
            if log is False:
                self.data[0][1:-1] /= self.physical.edgeFactor  # Chip's first row
                self.data[-1][1:-1] /= self.physical.edgeFactor # Chip's last row
                self.data[:, 0] /= self.physical.edgeFactor     # Chip's first column
                self.data[:, -1] /= self.physical.sideFactor    # Chip's last column

            else:
                self.data[0][1:-1] -= N.log10(self.physical.edgeFactor)  # Chip's first row
                self.data[-1][1:-1] -= N.log10(self.physical.edgeFactor) # Chip's last row
                self.data[:, 0] -= N.log10(self.physical.edgeFactor)     # Chip's first column
                self.data[:, -1] -= N.log10(self.physical.sideFactor)    # Chip's last column

        else:
            if log is False:
                self.data[0][1:-1] /= self.physical.edgeFactor     # Chip's first row
                self.data[-1][1:-1] /= self.physical.edgeFactor    # Chip's last row
                self.data[:, 0] /= self.physical.sideFactor  # Chip's first column
                self.data[:, -1] /= self.physical.edgeFactor # Chip's last column

            else:
                self.data[0][1:-1] -= N.log10(self.physical.edgeFactor)     # Chip's first row
                self.data[-1][1:-1] -= N.log10(self.physical.edgeFactor)    # Chip's last row
                self.data[:, 0] -= N.log10(self.physical.sideFactor)  # Chip's first column
                self.data[:, -1] -= N.log10(self.physical.edgeFactor) # Chip's last column

    def get_data(self, pixel, dacl=False):
        """Get the data value at a given pixel.

        :param pix: pixel coordinates (row, column) or pixel number (n)
        :param dacl: if True, get the dacl value of a given pixel
        """

        array = dacl and self.dacl or self.data

        if isinstance(pixel, int):
            assert len(array.ravel()) > pixel, \
                "The given pixel number should be < %i !" % len(array.ravel())

            return array.ravel()[pixel]

        elif isinstance(pixel, tuple) and len(pixel) == 2:
            assert array.shape[0] > pixel[0], \
                "The given pixel row should be < %i !" \
                % array.shape[0]
            assert array.shape[1] > pixel[1], \
                "The given pixel column should be < %i !" \
                % array.shape[1]

            return array[pixel[0], pixel[1]]

        else:
            raise ValueError('Bad pixel value !')

    def statistics(self, dacl=False, level=None):
        """Compute the statistics for the data.

        :param dacl: if True, get the staistics for the dacl array
        :param level: Level from which the signal becomes noise
        """

        if dacl is True:
            data = self.dacl

        else:
            data = self.data

        # self.stat = type('', (), {})()
        if level is not None:
            self.stat = namedtuple('stat', 'row col valid')

        else:
            self.stat = namedtuple('stat', 'row col')

        # Initialization
        self.stat.row = namedtuple('row', 'mean std med nMAD min max')
        self.stat.col = namedtuple('col', 'mean std med nMAD min max')
        if level is not None:
            self.stat.valid = namedtuple('valid', 'mean std med nMAD min max')

        # Global statistics
        self.stat.mean, self.stat.std = N.mean(data), N.std(data)
        self.stat.median, self.stat.nMAD = median_stats(data)
        self.stat.min, self.stat.max = N.min(data), N.max(data)

        # Row statistics
        self.stat.row.mean = N.mean(data, axis=0)
        self.stat.row.std = N.std(data, axis=0)
        self.stat.row.median, self.stat.row.nMAD = median_stats(data, axis=0)
        self.stat.row.min = N.min(data, axis=0)
        self.stat.row.max = N.max(data, axis=0)

        # Column statistics
        self.stat.col.mean = N.mean(data, axis=1)
        self.stat.col.std = N.std(data, axis=1)
        self.stat.col.median, self.stat.col.nMAD = median_stats(data, axis=1)
        self.stat.col.min = N.min(data, axis=1)
        self.stat.col.max = N.max(data, axis=1)

        # Statistics on valid pixels
        self.no_signal = data[data == 0]

        if level is not None:
            self.noisy = data[data >= level]
            bad = self.no_signal | self.noisy

            self.stat.valid.mean = N.mean(data[~bad])
            self.stat.valid.std = N.std(data[~bad])
            self.stat.valid.median, self.stat.valid.nMAD = median_stats(data[~bad])
            self.stat.valid.min = N.min(data[~bad])
            self.stat.valid.max = N.max(data[~bad])

    def plot(self, ax=None, figsize=(6, 9), cmap=jet,
             title='', xlabel='X [pixels]', ylabel='Y [pixels]',
             clabel='[arbitrary units]', dacl=False,
             vmin=None, vmax=None):
        """Display the image.

        :param ax: `matplotlib.axes.AxesSubplot` instance
        :param figsize: Size of the requiered figure (x,y)
        :param cmap: Matplotlib colormap
        :param title: Figure title
        :param xlabel: X-axis label
        :param ylabel: Y-axis label
        :param clabel: Color bar label
        :param dacl: if True, plot the dacl instead of the data
        :param vmin: Minimum scale value
        :param vmax: Maximum scale value
        """

        import matplotlib.pyplot as P

        # Figure and axis
        if ax is None:
            fig = P.figure(figsize=figsize)
            ax = fig.add_subplot(1, 1, 1,
                                 title=title,
                                 xlabel=xlabel,
                                 ylabel=ylabel)

        else:
            fig = ax.get_figure()

        # 2D image
        array = dacl and self.dacl or self.data
        img = ax.imshow(array, vmin=vmin, vmax=vmax,
                        cmap=cmap, aspect='equal',origin='lower',
                        interpolation='nearest')

        # Ticks and labels
        ax.tick_params(width=0)
        ax.set_xticks([-0.5, self.shape[1] - 0.5])
        ax.set_yticks([-0.5, self.shape[0] - 0.5])

        ax.set_xticklabels(['1', str(self.shape[1])])
        ax.set_yticklabels(['1', str(self.shape[0])])

        # Colorbar
        try:
            from mpl_toolkits.axes_grid1 import make_axes_locatable
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.2)
            cb = fig.colorbar(img, cax=cax)
            cb.set_label(clabel)

            return ax, cb

        except ImportError:
            from mpl_toolkits.axes_grid import make_axes_locatable
            divider = make_axes_locatable(ax)
            tmp_ax = divider.new_horizontal(size="5%", pad=0.2, pack_start=False)
            cax = divider._fig.add_axes(tmp_ax)
            cb = fig.colorbar(img, cax=cax)
            cb.set_label(clabel)

            return ax, cb

class Module(object):
    """Intermediate brick of a detector."""

    def __init__(self,
                 chipShape=(120, 80),
                 nChips=7,
                 append='horizontal',
                 name=u'M_#',
                 pixelSize=(0.13, 0.13),
                 edgeFactor=1.2,
                 sideFactor=32/13.,
                 moduleRing=0.55):
        """Class 'Module' initialization.

        :param chipShape: Chip shape (rows,columns) [pixels]
        :param nChips: Number of chips in a module
        :param append: Chips assembly in a module ['horizontal' | 'vertical']
        :param name: Nomenclature of the module (i.e. 'M_2')
        :param pixelSize: Pixels size (y,x) [mm]
        :param edgeFactor: Magnification factor for the pixels at the module's edges (y,x)
        :param sideFactor: Magnification factor for the wider pixels at the chip sides (y)
        :param moduleRing: module guard ring [mm]
        """

        # Inputs
        self.chipShape = chipShape # [pixels]
        self.nChips = nChips
        self.append = append
        self.name = name

        # Chip types (physical location in the module)
        if nChips > 1:
            location = ['core'] * nChips
            location[0] = 'left'
            location[-1] = 'right'

        else:
            location = ['single']

        # Module's chips (Chip instances)
        self.chips = N.array([Chip(chipShape=chipShape,
                                   location=location[i],
                                   name=u'C_%i' % (i+1))
                              for i in range(nChips)])

        # Module shape [pixels] & physical size [mm x mm]
        assert append in ['horizontal', 'vertical'], \
            "Bad assembly orientation: should be 'horizontal' or 'vertical'."

        if append == 'horizontal':      # Increase the number of columns
            self.shape = (chipShape[0], chipShape[1] * nChips)

        else:                           # Increase the number of rows
            self.shape = (chipShape[0] * nChips, chipShape[1])

        self.nPixels = self.shape[0] * self.shape[1]

        # Data initialization
        self.data = N.zeros(self.shape)

        # DacL initialization
        self.dacl = N.zeros(self.shape)

        # DacG initialization
        self.dacG = OrderedDict([('CMOS_DSBL', []),
                                 ('AMP_TP', []),
                                 ('ITHH', []),
                                 ('VADJ', []),
                                 ('VREF', []),
                                 ('IMFP', []),
                                 ('IOTA', []),
                                 ('IPRE', []),
                                 ('ITHL', []),
                                 ('ITUNE', []),
                                 ('IBUFFER', [])])

        for chip in self.chips:
            for key, value in chip.dacG.iteritems():
                self.dacG[key].append(value)

        # Chips indices
        self.chips_idx = N.zeros((nChips, self.shape[0],
                                  self.shape[1]), dtype='bool')

        for i in xrange(len(self.chips_idx)):
            if append == 'horizontal':
                self.chips_idx[i, :, i * chipShape[1]:
                               (i+1) * chipShape[1]] = True
            else:
                self.chips_idx[i, i * chipShape[0]:
                               (i+1) * chipShape[0], :] = True

        # Initialize the detector physical description
        self.physical = namedtuple('physical', 'pixelSize edgeFactor sideFactor \
                                   moduleGap moduleRing size shape data dacl')

        self.physical.pixelSize = pixelSize
        self.physical.edgeFactor = edgeFactor
        self.physical.sideFactor = sideFactor
        self.physical.moduleRing = moduleRing

        if append == 'horizontal':      # Increase the number of columns
            size = (self.chips[0].physical.size[0],
                    N.sum([chip.physical.size[1] for chip in self.chips]))

        else:                           # Increase the number of rows
            size = NotImplementedError

        self.physical.size = tuple(N.round(size, 2)) # physical size [mm x mm]

        nRows = N.round(self.physical.size[0] / self.physical.pixelSize[0])
        nCols = N.round(self.physical.size[1] / self.physical.pixelSize[1])

        self.physical.shape = (nRows, nCols)
        self.physical.data = N.zeros(self.physical.shape)
        self.physical.dacl = N.zeros(self.physical.shape)

    def __str__(self):
        """Print out Class information.

        :return: return Module information
        """

        s = "Module '%s': (%i,%i) [pixels]" % \
            (self.name, self.shape[0], self.shape[1])

        return s

    def update_dacG(self, keyword, values):
        """Update the dacG values (configuration file header) for
        each chip of the module.

        :param keyword: Configuration Keyword
        :param values: Corresponding keyword value for each chip
        """

        # Check input array consistency
        if isinstance(keyword, str):
            if isinstance(values, N.ndarray):
                assert len(values) == self.nChips, \
                    "There should be as many values as chips !"

                # Update dacG values in each chips
                for value, chip in zip(values, self.chips):
                    chip.update_dacG(keyword, value)

                # Update dacG values in the module
                self.dacG[keyword] = values

            else:
                raise IOError('The input array should be a 1D-numpy.ndarray !')

        else:
            raise IOError('The input keyword should be a string !')

    def config_header(self):
        """Create the module configuration header including thedacG
        values.

        :return: string containing the module configuration header
        """

        items = N.array([chip.dacG.items() for chip in self.chips.flatten()])

        keywords = items.T[0, :, 0]
        values = items.T[1]

        s = ''
        for i, keyword in enumerate(keywords):
            s += '%s ' % keyword
            for value in values[i]:
                s += '%s ' % value
            s += '\n'

        return s

    def write_config(self, filename, mask=False):
        """Create a module configuration file.

        :param filename: Output file name (output.dacs)
        :param mask: If True, append a mask of bad data to the
                     calibration output file
        """

        path, name = os.path.split(filename)
        basename = os.path.splitext(name)[0]

        if path == '':
            path = os.path.curdir

        fullname = os.path.extsep.join(
            (os.path.join(path, basename), 'txt'))
        slave = open(fullname, "w")  # open `.dacs` file

        # Write file header
        s = self.config_header()
        slave.write(s)

        # Write the dacl array
        slave.write("\nDACL\n")
        for i, chip in enumerate(self.chips):
            dacl = chip.dacl.copy()   # Reassign NaN to 0 for the config file
            dacl[~N.isfinite(dacl)] = 0

            for j, row in enumerate(N.array(dacl, dtype='i')):
                slave.write('%i %i ' % (i+1, j+1))
                row.tofile(slave, sep=' ', format='%s')
                slave.write('\n')

        # Write the mask array
        if mask is True:
            slave.write("\nMASK\n")
            for i, chip in enumerate(self.chips):
                bad = ~N.isfinite(chip.dacl)

                for j, row in enumerate(N.array(bad, dtype='i')):
                    slave.write('%i %i ' % (i+1, j+1))
                    row.tofile(slave, sep=' ', format='%s')
                    slave.write('\n')

        slave.close()

    def read_config(self, filename, mask=False):
        """Read a module configuration file.

        :param filename: Output file name (output.dacs)
        :param mask: If True, append a mask of bad data to the
                     calibration output file
        """

        module = open(filename, "r")  # open '.dacs' or '.txt' file

        # DacG dictionary
        header = []
        newline = module.readline()
        while newline != '\n':
            header.append(newline)
            newline = module.readline()

        dacG = OrderedDict()
        for line in header:
            dacG[line.split()[0]] = N.array(line.split()[1:], dtype='i')

        # DacL array
        newline = module.readline()
        while newline != 'DACL\n':
            newline = module.readline()

        lines = module.readlines()
        data = N.array([line.split()[2:] for line in lines], dtype='i')
        chips = N.array_split(data, self.nChips)
        dacL = N.hstack(chips)

        module.close()

        if mask:
            raise NotImplementedError

        else:
            return dacG, dacL

    def set_data(self, array, dacl=False):
        """Set the data array (or the dacl array) of the module.

        :param array: data array (module shape)
        :param dacl: if True, set the dacl array of the module
        """

        # Check input array consistency
        if isinstance(array, N.ndarray):
            assert array.shape == self.shape, \
                "The input array should have the following shape: (%i, %i) !" \
                % (self.shape[0], self.shape[1])

            # Fill the chips
            for i, chip in enumerate(self.chips):
                chip.set_data(N.reshape(array[self.chips_idx[i]],
                                        chip.shape), dacl=dacl)

            # Fill the module
            if dacl:
                self.dacl = array

            else:
                self.data = array

        else:
            raise IOError('The input array should be a 2D-numpy.ndarray !')

    def adjust_data(self, log=False):
        """Adjust the intensity of the 320 µm wide pixels (i.e. divided
        it by 32/13.) and the module edge pixels.

        :param log: data in log
        """

        # Adjust chips
        for chip in self.chips:
            chip.adjust_data(log)

        # Propage the modifications of the chips on the module
        if self.append == 'horizontal':
            self.data = N.hstack([chip.data for chip in self.chips])

        else:
            self.data = N.vstack([chip.data for chip in self.chips])

    def plot(self, ax=None, figsize=(13, 3.5), cmap=jet,
             title='', xlabel='X [pixels]', ylabel='Y [pixels]',
             clabel='[arbitrary units]', dacl=False,
             vmin=None, vmax=None):
        """Display the image.

        :param ax: `matplotlib.axes.AxesSubplot` instance
        :param figsize: Size of the requiered figure (x,y)
        :param cmap: Matplotlib colormap
        :param title: Figure title
        :param xlabel: X-axis label
        :param ylabel: Y-axis label
        :param clabel: Color bar label
        :param dacl: if True, plot the dacl instead of the data
        :param vmin: Minimum scale value
        :param vmax: Maximum scale value
        """

        from matplotlib.ticker import MaxNLocator
        import matplotlib.pyplot as P

        # Figure and axis
        if ax is None:
            fig = P.figure(figsize=figsize)
            ax = fig.add_subplot(1, 1, 1,
                                 title=title,
                                 xlabel=xlabel,
                                 ylabel=ylabel)

        else:
            fig = ax.get_figure()

        # 2D image
        array = dacl and self.dacl or self.data
        img = ax.imshow(array, vmin=vmin, vmax=vmax,
                        cmap=cmap, aspect='equal',origin='lower',
                        interpolation='nearest')

        # Ticks and labels
        ax.tick_params(width=0)
        tmp = N.arange(-0.5, self.shape[1] + self.chipShape[1] - 0.5, self.chipShape[1])
        ax.set_xticks(tmp)
        tmp = N.arange(0, self.shape[1] + self.chipShape[1], self.chipShape[1])
        tmp[0] = 1
        ax.set_xticklabels(tmp)
        ax.set_yticks([-0.5, self.shape[0] - 0.5])
        ax.set_yticklabels(['1', str(self.shape[0])])

        # Grid
        ax.grid(axis='x')

        # Colorbar
        try:
            from mpl_toolkits.axes_grid1 import make_axes_locatable
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.2)
            cb = fig.colorbar(img, cax=cax)
            cb.set_label(clabel)
            cb.ax.yaxis.set_major_locator(MaxNLocator(4))

            return ax, cb

        except ImportError:
            from mpl_toolkits.axes_grid import make_axes_locatable
            divider = make_axes_locatable(ax)
            tmp_ax = divider.new_horizontal(size="5%", pad=0.2, pack_start=False)
            cax = divider._fig.add_axes(tmp_ax)
            cb = fig.colorbar(img, cax=cax)
            cb.set_label(clabel)
            cb.ax.yaxis.set_major_locator(MaxNLocator(4))

            return ax, cb

    def statistics(self, dacl=False, level=None):
        """Compute the statistics for the data.

        :param dacl: if True, get the staistics for the dacl array
        :param level: Level from which the signal becomes noise
        """

        if dacl is True:
            data = self.dacl

        else:
            data = self.data

        # self.stat = type('', (), {})()
        if level is not None:
            self.stat = namedtuple('stat', 'row col valid')

        else:
            self.stat = namedtuple('stat', 'row col')

        # Initialization
        self.stat.row = namedtuple('row', 'mean std med nMAD min max')
        self.stat.col = namedtuple('col', 'mean std med nMAD min max')
        if level is not None:
            self.stat.valid = namedtuple('valid', 'mean std med nMAD min max')

        # Global statistics
        self.stat.mean, self.stat.std = N.mean(data), N.std(data)
        self.stat.median, self.stat.nMAD = median_stats(data)
        self.stat.min, self.stat.max = N.min(data), N.max(data)

        # Row statistics
        self.stat.row.mean = N.mean(data, axis=0)
        self.stat.row.std = N.std(data, axis=0)
        self.stat.row.median, self.stat.row.nMAD = median_stats(data, axis=0)
        self.stat.row.min = N.min(data, axis=0)
        self.stat.row.max = N.max(data, axis=0)

        # Column statistics
        self.stat.col.mean = N.mean(data, axis=1)
        self.stat.col.std = N.std(data, axis=1)
        self.stat.col.median, self.stat.col.nMAD = median_stats(data, axis=1)
        self.stat.col.min = N.min(data, axis=1)
        self.stat.col.max = N.max(data, axis=1)

        # Statistics on valid pixels
        self.no_signal = data[data == 0]

        if level is not None:
            self.noisy = data[data >= level]
            bad = self.no_signal | self.noisy

            self.stat.valid.mean = N.mean(data[~bad])
            self.stat.valid.std = N.std(data[~bad])
            self.stat.valid.median, self.stat.valid.nMAD = median_stats(data[~bad])
            self.stat.valid.min = N.min(data[~bad])
            self.stat.valid.max = N.max(data[~bad])


class Detector(object):
    """Full detector constructor."""

    def __init__(self,
                 chipShape=(120, 80),
                 nChips=7,
                 nModules=8,
                 append='vertical',
                 name='D5',
                 pixelSize=(0.13, 0.13),
                 edgeFactor=1.2,
                 sideFactor=32/13.,
                 moduleGap=3.57,
                 moduleRing=0.55):
        """Class initialization.

        :param chipShape: Chip shape (rows, columns) [pixels]
        :param nChip: Number of chips in a module
        :param nModule: Number of modules in the detector
        :param append: Chips assembly in a module ['horizontal' | 'vertical']
        :param name: Detector name (e.g. 'D5', 'S70', ...)
        :param pixelSize: Pixels size (y,x) [mm]
        :param edgeFactor: Magnification factor for the pixels at the module's edges (y,x)
        :param sideFactor: Magnification factor for the wider pixels at the chip sides (y)
        :param moduleGap: vertical gap between modules [mm]
        :param moduleRing: module guard ring [mm]
        """

        # Inputs
        self.chipShape = chipShape
        self.nChips = nChips
        self.nModules = nModules
        self.append = append
        self.name = name

        # Detector's modules (Module instances)
        self.modules = N.array([Module(chipShape=chipShape,
                                       nChips=nChips,
                                       name=u'M_%i' % (i+1))
                                for i in range(nModules)])

        # Detector's chips (Chip instances)
        self.chips = N.array([module.chips for module in self.modules])

        # Detector shape [pixels] & physical size [mm x mm]
        assert append in ['horizontal', 'vertical'], \
            "Bad assembly orientation: should be 'horizontal' or 'vertical'."


        if append == 'horizontal':      # Increase the number of columns
            self.shape = (self.modules[0].shape[0],
                          self.modules[0].shape[1] * self.nModules)

        else:                           # Increase the number of rows
            self.shape = (self.modules[0].shape[0] * self.nModules,
                          self.modules[0].shape[1])

        self.nPixels = self.shape[0] * self.shape[1]

        # Data initialization
        self.data = N.zeros(self.shape)

        # DacL initialization
        self.dacl = N.zeros(self.shape)

        # DacG initialization
        self.dacG = OrderedDict([('CMOS_DSBL', []),
                                 ('AMP_TP', []),
                                 ('ITHH', []),
                                 ('VADJ', []),
                                 ('VREF', []),
                                 ('IMFP', []),
                                 ('IOTA', []),
                                 ('IPRE', []),
                                 ('ITHL', []),
                                 ('ITUNE', []),
                                 ('IBUFFER', [])])

        for mod in self.modules:
            for key, values in mod.dacG.iteritems():
                self.dacG[key].append(values)

        # Modules indices
        self.modules_idx = N.zeros((self.nModules, self.shape[0],
                                    self.shape[1]), dtype='bool')

        for i in xrange(len(self.modules_idx)):
            if append == 'horizontal':
                self.modules_idx[i, :, i * self.modules[0].shape[1]:
                                 (i+1) * self.modules[0].shape[1]] = True
            else:
                self.modules_idx[i, i * self.modules[0].shape[0]:
                                 (i+1) * self.modules[0].shape[0], :] = True

        # Chips indices
        self.chips_idx = N.zeros((self.nModules, self.nChips,
                                    self.shape[0], self.shape[1]),
                                   dtype='bool')

        for i in xrange(self.nModules):
            for j in xrange(self.nChips):
                self.chips_idx[i, j,
                               i * self.modules[i].shape[0]:
                               (i+1) * self.modules[i].shape[0],
                               j * chipShape[1]:
                               (j+1) * chipShape[1]] = True

        # Initialize the detector physical description
        self.physical = namedtuple('physical', 'pixelSize edgeFactor sideFactor \
                                   moduleGap moduleRing size shape data dacl')

        self.physical.pixelSize = pixelSize
        self.physical.edgeFactor = edgeFactor
        self.physical.sideFactor = sideFactor
        self.physical.moduleGap = moduleGap
        self.physical.moduleRing = moduleRing

        if append == 'horizontal':      # Increase the number of columns
            size = NotImplementedError

        else:                           # Increase the number of rows
            if self.nModules > 1:
                size = (N.sum([mod.physical.size[0] for mod in self.modules]) + \
                        (self.nModules - 1) * self.physical.moduleGap - \
                        self.physical.moduleRing, self.modules[0].physical.size[1])

            else:
                size = (N.sum([mod.physical.size[0] for mod in self.modules]),
                        self.modules[0].physical.size[1])

        self.physical.size = tuple(N.round(size, 2)) # physical size [mm x mm]

        nRows = N.round(self.physical.size[0] / self.physical.pixelSize[0])
        nCols = N.round(self.physical.size[1] / self.physical.pixelSize[1])

        self.physical.shape = (nRows, nCols)
        self.physical.data = N.zeros(self.physical.shape)
        self.physical.dacl = N.zeros(self.physical.shape)

    def __str__(self):
        """Print out Class information.

        :return: return Detector information
        """

        s = '    * name: %s\n' % \
            CT.strc('%s' % self.name, 'blue')
        s += '    * shape: (%s,%s) pixels\n' % \
             (CT.strc('%i' % self.shape[0], 'blue'),
              CT.strc('%i' % self.shape[1], 'blue'))
        s += '    * %s modules\n' % \
             CT.strc('%i' % self.nModules, 'blue')
        s += '    * %s chips per module\n' % \
             CT.strc('%i' % self.nChips, 'blue')
        s += '    * chip shape: (%s,%s) pixels' % \
             (CT.strc('%i' % self.chipShape[0], 'blue'),
              CT.strc('%i' % self.chipShape[1], 'blue'))
        return s

    def update_dacG(self, keyword, values):
        """Update the dacG values (configuration file header) for
        each chip of the detector.

        :param keyword: Configuration Keyword
        :param values: Corresponding keyword values for each module
        """

        # Check input array consistency
        if isinstance(keyword, str):
            if isinstance(values, N.ndarray):
                assert values.shape == (self.nModules, self.nChips), \
                    "The input array shape should be (%i, %i) !" % \
                    (self.nModules, self.nChips)

                # Update dacG values in each module (and therefore in
                # each chip)
                for array, module in zip(values, self.modules):
                    module.update_dacG(keyword, array)

                # Update dacG values in the detector
                self.dacG[keyword] = values

            else:
                raise IOError('The input array should be a 2D-numpy.ndarray !')

        else:
            raise IOError('The input keyword should be a string !')

    def read_cppm(self, filename):
        """Read file containing a single chip data array using the CPPM
        data format.

        :param filename: File to be read containing a single chip data array
        """

        infile = open(filename, "r")
        header = infile.readline()   # First line

        # DacG
        keys, values = N.array([hdr.split(' = ')
                                for hdr in header.split(', ')]).T
        keys[0] = keys[0][20:]
        values = N.array(values, dtype='i')

        dacG = {}
        for key, value in zip(keys, values):
            dacG[key] = value

        # DacL
        dacL =  N.array([line.split()[0]
                         for line in infile.readlines()], dtype='f')

        infile.close()

        return dacG, dacL.reshape(self.shape)

    def write_config(self, filename, mask=False):
        """Create the detector configuration file.

        :param filename: Output file name (output.master)
        :param mask: If True, append a mask of bad data to the
                     calibration output file
        """

        path, name = os.path.split(filename)
        basename = os.path.splitext(name)[0]

        if path == '':
            path = os.path.curdir

        # `.master` file ------------------------------

        fullname = os.path.extsep.join(
            (os.path.join(path, basename), 'master'))
        master = open(fullname, "w")  # open `.master` file

        # Loop over modules (write the full address of the associated
        # `.dacs` files)
        module_filenames = []
        for n in range(self.nModules):
            module_filename = os.path.extsep.join((
                os.path.join(path, basename) + '_%i' % (n + 1), 'txt'))
            module_filenames.append(module_filename)
            master.write("mod_%i %s\n" % (n + 1, module_filename))

        master.close()              # close `.master` file

        # `.dacs` files ------------------------------

        for i, module in enumerate(self.modules):
            module.write_config(module_filenames[i], mask)

    def read_config(self, filename, mask=False, set_dac=False):
        """Read a detector configuration file.

        :param filename: Output file name (output.master)
        :param mask: If True, append a mask of bad data to the
                     calibration output file
        :param set_dac: Automatically set the detector's DAC values
        """

        master = open(filename, "r")  # open `.master` file
        dacs = master.readlines()
        filenames = [dac.split()[-1] for dac in dacs]

        # DacL array
        dacL = []
        headers = []
        for module,name  in zip(self.modules, filenames):
            header, data = module.read_config(name, mask)
            headers.append(header)
            dacL.append(data)

        dacL = N.vstack(dacL)

        # DacG dictionary
        dacG = headers[0]
        for header in headers[1:]:
            for key, array in header.iteritems():
                dacG[key] = N.vstack((dacG[key], array))

        master.close()              # close `.master` file

        if set_dac is True:
            self.set_data(array=dacL, dacl=True)
            for keyword, values in dacG.iteritems():
                self.update_dacG(keyword, values)

        else:
            return dacG, dacL

    def set_data(self, array=None, dacl=False,
                 test=False, init=False, random=False):
        """Set the data array (or the dacl array) of the detector.

        :param array: data array (detector shape)
        :param dacl: if True, set the dacl array of the detector with
                     given values.
        :param test: if True, set the dacl/data array of the detector with a
                     fake pattern of concentric circles.
        :param init: if True, set the dacl/data array of the detector with a
                     flat field.
        :param random: if True, set the dacl/data array of the detector with a
                       random field.
        """

        if init:
            # Flat values (32)
            array = N.ones(self.shape, dtype='int') * 32

        if test:
            # Concentric circles of different values
            lims = (self.shape[0] / 2, self.shape[1] / 2)
            x, y = N.mgrid[-lims[0]:lims[0]:1, -lims[1]:lims[1]:1]
            z = N.sqrt(x**2 + y**2) + N.cos(x**2 + y**2)
            array = N.round(z / 10 + 5)

        if random:
            array = N.random.randint(0, 1e4, size=self.shape)

        # Check input array consistency
        if isinstance(array, N.ndarray):
            assert array.shape == self.shape, \
                "The input array should have the following shape: (%i, %i) !" \
                % (self.shape[0], self.shape[1])

            # Fill the modules
            for i, module in enumerate(self.modules):
                module.set_data(N.reshape(array[self.modules_idx[i]],
                                          module.shape), dacl=dacl)

            # Fill the detector
            if dacl:
                self.dacl = array

            else:
                self.data = array

            # # Fill the physical detector
            # self.physical.set_data(self)

        else:
            raise IOError('The input array should be a 2D-numpy.ndarray !')

    def adjust_data(self, log=False):
        """Adjust the intensity of the 320 µm wide pixels (i.e. divided
        it by 32/13.) and the module edge pixels.

        :param log: data in log
        """

        # Adjust the modules
        for module in self.modules:
            module.adjust_data(log)

        # Propagate the modifications of the modules on the detector
        if self.append == 'vertical':
            self.data = N.vstack([mod.data for mod in self.modules])

        else:
            self.data = N.hstack([mod.data for mod in self.modules])

    def set_physical_data(self):
        """Set the detector data into the physical model."""

        cols = self.chipShape[0]
        rows = self.chipShape[1]

        for i in xrange(self.nModules):
            if self.nModules > 1:
                m_gap = i * (self.physical.shape[0] - self.nModules * cols) /\
                        (self.nModules - 1)

            else:
                m_gap = 0

            for j in xrange(self.nChips):
                if self.nChips > 1:
                    c_gap = j * (self.physical.shape[1] - self.nChips * rows) /\
                            (self.nChips - 1)

                else:
                    c_gap = 0

                self.physical.dacl[i * cols + m_gap:(i+1) * cols + m_gap,
                                       j * rows + c_gap:(j+1) * rows + c_gap] = \
                                       self.modules[i].chips[j].dacl

                self.physical.data[i * cols + m_gap:(i+1) * cols + m_gap,
                                       j * rows + c_gap:(j+1) * rows + c_gap] = \
                                       self.modules[i].chips[j].data

    def reshape_pixels(self):
        """Reshape the 320 µm wide inter-chip pixels.

        Side pixels of each 'core' chip are 320 µm wide. Last and first
        column of side by side chips can then be represented as 5
        'classic' pixel columns: 2 x 320 ~ 5 x 130.

        .. todo::
            * Robustify (i.e. if the geometry is different)

        .. note::
            * The 320 µm wide pixel's intensity has to be adjusted
              (i.e. divided by 32/13.) before it is reshaped into 2.5
              pixels.
       """

        rows = self.chipShape[1]

        # Among the 5 pixel columns: the first two are identical and
        # equal to the first, the last two are identical and equal to
        # the last column, and the middle column corresponds to the mean
        # of the first and last columns.
        for j in xrange(self.nChips - 1):
            c_gap = j * (self.physical.shape[1] - self.nChips * rows) /\
                    (self.nChips - 1)

            # 2nd column is identical to the first one (which is
            # actually the last column of the chip)
            self.physical.data[:, (j+1) * rows + c_gap] = \
                self.physical.data[:, (j+1) * rows + c_gap - 1]

            # Middle column is the average of the first and last column
            self.physical.data[:, (j+1) * rows + c_gap + 1] = \
                N.mean([self.physical.data[:, (j+1) * rows + c_gap - 1],
                        self.physical.data[:, (j+1) * rows + c_gap + 3]], axis=0)

            # second to last column is identical to the last one (which
            # is actually the first column of the chip)
            self.physical.data[:, (j+1) * rows + c_gap + 2] = \
                self.physical.data[:, (j+1) * rows + c_gap + 3]

    def chip2det(self, M=1, C=1, row=1, col=1):
        """Compute the pixel coordinates (Row,Col) in the detector from
        the coordinates of a pixel in a given chip (M,C,row,col).

        :param M: Module number  (1 <= M <= nModules)
        :param C: Chip number (1 <= C <= nChips)
        :param row: Row number in the chip [pixel]
        :param col: Column number in the chip [pixel]

        :return: Row and Column in the detector (Row,Col)
        """

        # Check the options concistency
        assert (M >= 1) and (M <= self.nModules), \
            "Module number should be between 1 and %i!" % self.nModules
        assert (C >= 1) and (C <= self.nChips), \
            "Chip number should be between 1 and %i!" % self.nChips

        assert (row >= 1) and (row <= self.chipShape[0]), \
            "Row number should be between 1 and %i!" % \
            self.chipShape[0]

        assert (col >= 1) and (col <= self.chipShape[1]), \
            "Column number should be between 1 and %i!" % \
            self.chipShape[1]

        assert isinstance(M, int), "Module number should be an integer!"
        assert isinstance(C, int), "Chip number should be an integer!"
        assert isinstance(row, int), "row number should be an integer!"
        assert isinstance(col, int), "column number should be an integer!"

        Row = (M - 1) * self.chipShape[0] + row
        Col = (C - 1) * self.chipShape[1] + col

        return (Row, Col)

    def statistics(self, dacl=False, level=None):
        """Compute the statistics for the data.

        :param dacl: if True, get the staistics for the dacl array
        :param level: Level from which the signal becomes noise
        """

        if dacl is True:
            data = self.dacl

        else:
            data = self.data

        # self.stat = type('', (), {})()
        if level is not None:
            self.stat = namedtuple('stat', 'row col valid')

        else:
            self.stat = namedtuple('stat', 'row col')

        # Initialization
        self.stat.row = namedtuple('row', 'mean std med nMAD min max')
        self.stat.col = namedtuple('col', 'mean std med nMAD min max')
        if level is not None:
            self.stat.valid = namedtuple('valid', 'mean std med nMAD min max')

        # Global statistics
        self.stat.mean, self.stat.std = N.mean(data), N.std(data)
        self.stat.median, self.stat.nMAD = median_stats(data)
        self.stat.min, self.stat.max = N.min(data), N.max(data)

        # Row statistics
        self.stat.row.mean = N.mean(data, axis=0)
        self.stat.row.std = N.std(data, axis=0)
        self.stat.row.median, self.stat.row.nMAD = median_stats(data, axis=0)
        self.stat.row.min = N.min(data, axis=0)
        self.stat.row.max = N.max(data, axis=0)

        # Column statistics
        self.stat.col.mean = N.mean(data, axis=1)
        self.stat.col.std = N.std(data, axis=1)
        self.stat.col.median, self.stat.col.nMAD = median_stats(data, axis=1)
        self.stat.col.min = N.min(data, axis=1)
        self.stat.col.max = N.max(data, axis=1)

        # Statistics on valid pixels
        self.no_signal = data[data == 0]

        if level is not None:
            self.noisy = data[data >= level]
            bad = self.no_signal | self.noisy

            self.stat.valid.mean = N.mean(data[~bad])
            self.stat.valid.std = N.std(data[~bad])
            self.stat.valid.median, self.stat.valid.nMAD = median_stats(data[~bad])
            self.stat.valid.min = N.min(data[~bad])
            self.stat.valid.max = N.max(data[~bad])

    def det2chip(self, Row=1, Col=1):
        """Compute the pixel coordinates (M,C,row,col) in the module M
        and the chip C from the coordinates of a pixel (Row,Col) in the
        detector.

        :param Row: Pixel row number in the detector
        :param Col: Pixel column number in the detector

        :return: Pixel coordinates (row,col) in the module M and the
                 chip C (M,C,row,col)
        """

        assert (Row >= 1) and (Row <= self.shape[0]), \
            "Pixel row number should be between 1 and %i!" % \
            self.shape[0]

        assert (Col >= 1) and (Col <= self.shape[1]), \
            "Pixel column number should be between 1 and %i!" % \
            self.shape[1]

        assert isinstance(Row, int), "Row number should be an integer!"
        assert isinstance(Col, int), "Column number should be an integer!"

        C = (Row + self.chipShape[0] - 1) // self.chipShape[0]
        M = (Col + self.chipShape[1] - 1) // self.chipShape[1]
        row = (Row + self.chipShape[0] - 1) % self.chipShape[0] + 1
        col = (Col + self.chipShape[1] - 1) % self.chipShape[1] + 1

        return (M, C, row, col)

    def get_module_array(self, array, M=1):
        """Retrieve the array corresponding to the indices of a given
        module .

        :param data: Any array with the same shape as the detector
        :param M: Module number (1 <= M <= nModules)

        :return: Array with the same shape as the chosen module
        """

        # Check the options concistency
        assert (M >= 1) and (M <= self.nModules), \
            "Module number should be between 1 and %i!" % self.nModules

        assert isinstance(M, int), \
            "Module number should be an integer!"

        M -= 1 # Module #1 is actually modules_idx[0]

        # Return data corresponding the requiered module
        return N.reshape(array[self.modules_idx[M]], self.modules[0].shape)

    def get_chip_array(self, array, M=1, C=1):
        """Retrieve the array corresponding to the indices of a given
        chip.

        :param data: Any array with the same shape as the detector
        :param M: Module number (1 <= M <= nModules)
        :param C: Chip number in the module (1 <= C <= nChips)

        :return: Array with the same shape as the chosen chip
        """

        # Check the options concistency
        assert (C >= 1) and (C <= self.nChips), \
            "Chip number should be between 1 and %i!" % self.nChips

        assert isinstance(C, int), \
            "Chip number should be an integer!"

        C -= 1 # Chip #1 is actually chips_idx[0]

        # Module
        mod = self.get_module_array(array, M)

        # Return data corresponding the requiered chip
        return N.reshape(mod[self.modules[0].chips_idx[C]],
                         self.chipShape)

    def plot(self, ax=None, figsize=(6, 9), cmap=jet,
             title='', xlabel='X [pixels]', ylabel='Y [pixels]',
             clabel='[arbitrary units]', dacl=False,
             vmin=None, vmax=None, image=None):
        """Display the image.

        :param ax: `matplotlib.axes.AxesSubplot` instance
        :param figsize: Size of the requiered figure (x,y)
        :param cmap: Matplotlib colormap
        :param title: Figure title
        :param xlabel: X-axis label
        :param ylabel: Y-axis label
        :param clabel: Color bar label
        :param dacl: if True, plot the dacl instead of the data
        :param vmin: Minimum scale value
        :param vmax: Maximum scale value
        """

        import matplotlib.pyplot as P
        import mpl_toolkits

        # Figure and axis
        if ax is None:
            fig = P.figure(figsize=figsize)
            ax = fig.add_subplot(1, 1, 1,
                                 title=title,
                                 xlabel=xlabel,
                                 ylabel=ylabel)

        else:
            fig = ax.get_figure()

        # 2D image
        if dacl == True:
            array = self.dacl

        elif image is not None:
            array = image

        else:
            array = self.data

        img = ax.imshow(array, vmin=vmin, vmax=vmax,
                        cmap=cmap, aspect='equal',origin='lower',
                        interpolation='nearest')

        # Ticks and labels
        ax.set_xticks(N.concatenate(([-0.5],
                                     N.arange(self.chipShape[1] - 0.5,
                                              self.shape[1],
                                              self.chipShape[1]))))

        ax.set_xticklabels(range(0, self.chipShape[1] + \
                                 self.shape[1],
                                 self.chipShape[1]))

        ax.set_yticks(N.concatenate(([-0.5],
                                     N.arange(self.chipShape[0] - 0.5,
                                              self.shape[0],
                                              self.chipShape[0]))))
        ax.set_yticklabels(range(0, self.chipShape[0] + \
                                 self.shape[0],
                                 self.chipShape[0]))

        # Grid
        ax.grid()

        # Colorbar
        try:
            from mpl_toolkits.axes_grid1 import make_axes_locatable
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.2)
            cb = fig.colorbar(img, cax=cax)
            cb.set_label(clabel)

            return ax, cb, img

        except ImportError:
            from mpl_toolkits.axes_grid import make_axes_locatable
            divider = make_axes_locatable(ax)
            tmp_ax = divider.new_horizontal(size="5%", pad=0.2, pack_start=False)
            cax = divider._fig.add_axes(tmp_ax)
            cb = fig.colorbar(img, cax=cax)
            cb.set_label(clabel)

            return ax, cb, img

    def physical_plot(self, ax=None, figsize=(6, 9), cmap=jet,
                      title='', xlabel='X [mm]', ylabel='Y [mm]',
                      clabel='Counts [arbitrary units]', dacl=False,
                      vmin=None, vmax=None, image=None):
        """Display the image.

        :param ax: `matplotlib.axes.AxesSubplot` instance
        :param figsize: Size of the requiered figure (x,y)
        :param cmap: Matplotlib colormap
        :param title: Figure title
        :param xlabel: X-axis label
        :param ylabel: Y-axis label
        :param clabel: Color bar label
        :param dacl: if True, plot the dacl instead of the data
        :param vmin: Minimum scale value
        :param vmax: Maximum scale value
        """

        import matplotlib.pyplot as P
        import mpl_toolkits

        # Figure and axis
        if ax is None:
            fig = P.figure(figsize=figsize)
            ax = fig.add_subplot(1, 1, 1,
                                 title=title,
                                 xlabel=xlabel,
                                 ylabel=ylabel)

        else:
            fig = ax.get_figure()

        # 2D image
        # 2D image
        if dacl == True:
            array = self.physical.dacl

        elif image is not None:
            array = image

        else:
            array = self.physical.data

        img = ax.imshow(array, vmin=vmin, vmax=vmax,
                        cmap=cmap, aspect='equal',origin='lower',
                        interpolation='nearest')

        # Ticks and ticks labels
        ax.set_xticks(N.linspace(0, self.shape[1], num=5))
        ax.set_xticklabels(['%.2f' % i for i in N.linspace(0, self.physical.size[1], num=5)])
        ax.set_yticks(N.linspace(0, self.shape[0], num=9))
        ax.set_yticklabels(['%.2f' % i for i in N.linspace(0, self.physical.size[0], num=9)])

        # Colorbar
        try:
            from mpl_toolkits.axes_grid1 import make_axes_locatable
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.2)
            cb = fig.colorbar(img, cax=cax)
            cb.set_label(clabel)

            return ax, cb, img

        except ImportError:
            from mpl_toolkits.axes_grid import make_axes_locatable
            divider = make_axes_locatable(ax)
            tmp_ax = divider.new_horizontal(size="5%", pad=0.2, pack_start=False)
            cax = divider._fig.add_axes(tmp_ax)
            cb = fig.colorbar(img, cax=cax)
            cb.set_label(clabel)

            return ax, cb, img

# Older python/numpy/matplotlib Compatibility ==================================

# Class copied from python 2.7 to ensure compatibility with python 2.6
from _abcoll import KeysView, ValuesView, ItemsView
from thread import get_ident as _get_ident
from collections import MutableMapping

class OrderedDict(dict):
    """Dictionary that remembers insertion order.

    .. warning::

        Class copied from python 2.7 to ensure compatibility with python 2.6
    """

    def __init__(self, *args, **kwds):
        """Initialize an ordered dictionary.  The signature is the same
        as regular dictionaries, but keyword arguments are not
        recommended because their insertion order is arbitrary.
        """

        if len(args) > 1:
            raise TypeError('expected at most 1 arguments, got %d' % len(args))
        try:
            self.__root
        except AttributeError:
            self.__root = root = []                     # sentinel node
            root[:] = [root, root, None]
            self.__map = {}
        self.__update(*args, **kwds)

    def __setitem__(self, key, value, PREV=0, NEXT=1, dict_setitem=dict.__setitem__):
        """od.__setitem__(i, y) <==> od[i]=y"""

        # Setting a new item creates a new link at the end of the linked
        # list, and the inherited dictionary is updated with the new
        # key/value pair.
        if key not in self:
            root = self.__root
            last = root[PREV]
            last[NEXT] = root[PREV] = self.__map[key] = [last, root, key]
        dict_setitem(self, key, value)

    def __delitem__(self, key, PREV=0, NEXT=1, dict_delitem=dict.__delitem__):
        """od.__delitem__(y) <==> del od[y]"""

        # Deleting an existing item uses self.__map to find the link
        # which gets removed by updating the links in the predecessor
        # and successor nodes.
        dict_delitem(self, key)
        link_prev, link_next, key = self.__map.pop(key)
        link_prev[NEXT] = link_next
        link_next[PREV] = link_prev

    def __iter__(self):
        """od.__iter__() <==> iter(od)"""

        # Traverse the linked list in order.
        NEXT, KEY = 1, 2
        root = self.__root
        curr = root[NEXT]
        while curr is not root:
            yield curr[KEY]
            curr = curr[NEXT]

    def __reversed__(self):
        """od.__reversed__() <==> reversed(od)"""

        # Traverse the linked list in reverse order.
        PREV, KEY = 0, 2
        root = self.__root
        curr = root[PREV]
        while curr is not root:
            yield curr[KEY]
            curr = curr[PREV]

    def clear(self):
        'od.clear() -> None.  Remove all items from od.'
        for node in self.__map.itervalues():
            del node[:]
        root = self.__root
        root[:] = [root, root, None]
        self.__map.clear()
        dict.clear(self)

    # -- the following methods do not depend on the internal structure --

    def keys(self):
        """od.keys() -> list of keys in od"""

        return list(self)

    def values(self):
        """od.values() -> list of values in od"""

        return [self[key] for key in self]

    def items(self):
        """od.items() -> list of (key, value) pairs in od"""

        return [(key, self[key]) for key in self]

    def iterkeys(self):
        """od.iterkeys() -> an iterator over the keys in od"""

        return iter(self)

    def itervalues(self):
        'od.itervalues -> an iterator over the values in od'
        for k in self:
            yield self[k]

    def iteritems(self):
        """od.iteritems -> an iterator over the (key, value) pairs in
        od
        """

        for k in self:
            yield (k, self[k])

    update = MutableMapping.update

    __update = update # let subclasses override update without breaking __init__

    __marker = object()

    def pop(self, key, default=__marker):
        """od.pop(k[,d]) -> v, remove specified key and return the
        corresponding value.  If key is not found, d is returned if
        given, otherwise KeyError is raised.
        """

        if key in self:
            result = self[key]
            del self[key]
            return result
        if default is self.__marker:
            raise KeyError(key)
        return default

    def setdefault(self, key, default=None):
        """od.setdefault(k[,d]) -> od.get(k,d), also set od[k]=d if k
        not in od
        """

        if key in self:
            return self[key]
        self[key] = default
        return default

    def popitem(self, last=True):
        """od.popitem() -> (k, v), return and remove a (key, value)
        pair.  Pairs are returned in LIFO order if last is true or FIFO
        order if false.
        """

        if not self:
            raise KeyError('dictionary is empty')
        key = next(reversed(self) if last else iter(self))
        value = self.pop(key)
        return key, value

    def __repr__(self, _repr_running={}):
        'od.__repr__() <==> repr(od)'
        call_key = id(self), _get_ident()
        if call_key in _repr_running:
            return '...'
        _repr_running[call_key] = 1
        try:
            if not self:
                return '%s()' % (self.__class__.__name__,)
            return '%s(%r)' % (self.__class__.__name__, self.items())
        finally:
            del _repr_running[call_key]

    def __reduce__(self):
        """Return state information for pickling"""

        items = [[k, self[k]] for k in self]
        inst_dict = vars(self).copy()
        for k in vars(OrderedDict()):
            inst_dict.pop(k, None)
        if inst_dict:
            return (self.__class__, (items,), inst_dict)

        return self.__class__, (items,)

    def copy(self):
        """od.copy() -> a shallow copy of od"""

        return self.__class__(self)

    @classmethod
    def fromkeys(cls, iterable, value=None):
        """OD.fromkeys(S[, v]) -> New ordered dictionary with keys from S.
        If not specified, the value defaults to None.
        """

        self = cls()
        for key in iterable:
            self[key] = value
        return self

    def __eq__(self, other):
        """od.__eq__(y) <==> od==y.  Comparison to another OD is
        order-sensitive while comparison to a regular mapping is
        order-insensitive.
        """

        if isinstance(other, OrderedDict):
            return len(self)==len(other) and self.items() == other.items()
        return dict.__eq__(self, other)

    def __ne__(self, other):
        """od.__ne__(y) <==> od!=y"""

        return not self == other

    # -- the following methods support python 3.x style dictionary views --

    def viewkeys(self):
        """od.viewkeys() -> a set-like object providing a view on od's
        keys
        """

        return KeysView(self)

    def viewvalues(self):
        """od.viewvalues() -> an object providing a view on od's
        values
        """

        return ValuesView(self)

    def viewitems(self):
        """od.viewitems() -> a set-like object providing a view on od's
        items
        """

        return ItemsView(self)

# End of libXpad.py ===========================================================
