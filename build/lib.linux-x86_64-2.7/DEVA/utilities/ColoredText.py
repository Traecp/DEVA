#!/usr/bin/env python
# -*- coding: utf-8 -*-
###############################################################################
## Filename:          ColoredText.py
## Description:       ANSII Color formatting for output in terminal
## Author:            Clement Buton <clement.buton@esrf.fr>
## Author:            $Author: cbuton $
## Created on:        $Date: 2013/03/06 17:48:52 (UTC)$
## Modified on:       2013/03/22 09:32:25
## Copyright:         2013, Clement Buton
## $Id: ColoredText.py, 2013/03/06 17:48:52 cbuton Exp $
###############################################################################

from __future__ import print_function  # be prepared for python 3.x

"""
Library for Cross-platform ANSII Color formatting for output in
terminal. It allows to print and write some text or draw a progress bar
in a terminal with given ANSII colors and styles.

This simple utility prevent of importing non-standard (although very
nice) libraries such as `colorama`
(https://pypi.python.org/pypi/colorama), `termcolor`
(https://pypi.python.org/pypi/termcolor) and `progressbar`
(http://code.google.com/p/python-progressbar/).
"""

__author__ = 'Clement Buton <clement.buton@esrf.fr>'
__date__ = '2013/03/06 17:48:52'
__adv__ = 'ColoredText.py'

import sys
import logging

colorCodes = {
    'normal': '',
    'black':  '30',
    'red':    '31',
    'green':  '32',
    'yellow': '33',
    'blue':   '34',
    'purple': '35',
    'cyan':   '36',
    'white':  '37'}

styleCodes = {
    'normal':     '0',
    'bold':       '1;',
    'faint':      '2;',  # Not widely supported
    'italic':     '3;',
    'underline':  '4;',
    'blink':      '5;',
    'fast_blink': '6;'}   # Not widely supported

# Class =======================================================================


class ProgressBar:
    """Creates a text-based progress bar."""

    def __init__(self, minValue=0, maxValue=100, totalWidth=80,
                 holds='[]', marker='=', markerEdge='>',
                 color='normal', style='normal'):
        """Initialize the `ProgressBar` class.

        :param minValue: Minimum value of the progress bar
        :param maxValue: Maximum value of the progress bar
        :param totalWidth: With of the progress bar
        :param holds: Sign holding the progress bar (e.g: '[]', '||', '<>', ...)
        :param marker: Progress marker (e.g: '█', '=', '#', ...)
        :param markerEdge: Marker edge (e.g: '', '<', '>', ...)
        :param color: ANSI color among ['normal','black','red','green',
                      'yellow','blue','purple','cyan','white']
        :param style: ANSI style among ['normal','bold','faint','italic',
                      'underline','blink','fast_blink']
        """

        self.min = minValue
        self.max = maxValue
        self.width = totalWidth
        self.progBar = holds  # This holds the progress bar string
        self.marker = marker    # '█', '=', '#', ...
        self.edge = markerEdge  # '', '<', '>', ...
        self.color = color
        self.style = style

        self.span = maxValue - minValue
        self.amount = 0         # When amount == max, we are 100% done
        self.updateAmount(0)    # Build progress bar string

    def updateAmount(self, newAmount=0):
        """Update the progress bar with the new amount (with min and max
        values set at initialization; if it is over or under, it takes
        the min or max value as a default.

        :param newAmount: Amount value for the progress bar update
        """

        if newAmount <= self.min:
            newAmount = self.min
        if newAmount >= self.max:
            newAmount = self.max

        self.amount = newAmount

        # Figure out the new percent done, round to an integer
        diffFromMin = float(self.amount - self.min)
        percentDone = (diffFromMin / float(self.span)) * 100.0
        percentDone = int(round(percentDone))

        # Figure out how many hash bars the percentage should be
        allFull = self.width - 2
        numHashes = (percentDone / 100.0) * allFull
        numHashes = int(round(numHashes))

        # Build a progress bar with an arrow of equal signs; special
        # cases for empty and full
        if numHashes == 0:
            self.progBar = "[%s%s]" % (self.edge, ' ' * (allFull - 1))
        elif numHashes == allFull:
            self.progBar = "[%s]" % (self.marker * allFull)
        else:
            self.progBar = "[%s%s%s]" % (self.marker * (numHashes - 1),
                                         self.edge,
                                         ' ' * (allFull - numHashes))

        # figure out where to put the percentage, roughly centered
        percentPlace = (len(self.progBar) / 2) - len(str(percentDone))
        percentString = str(percentDone) + "%"

        # slice the percentage into the bar
        self.progBar = ''.join([self.progBar[0:percentPlace], percentString,
                                self.progBar[percentPlace+len(percentString):]
                                ])

    def __str__(self):
        """Print out the progress bar."""

        return str(self.progBar)

    def __call__(self, value):
        """Updates the amount, and writes to stdout. Prints a carriage
        return first, so it will overwrite the current line in stdout.

        :param value: Amount value for the progress bar update
        """

        print('\r', end='')
        self.updateAmount(value)
        writec(str(self), self.color, self.style)
        sys.stdout.flush()


class Formatter(logging.Formatter):
    """'logging.Formatter' instance used to convert a LogRecord to
    text."""

    def format(self, record):
        """Format the specified record as text."""

        level_colors = {
            'DEBUG': strc('DEBUG', 'yellow', 'bold'),
            'INFO': strc('INFO', 'blue', 'bold'),
            'WARNING': strc('WARNING', 'yellow', 'bold'),
            'ERROR': strc('ERROR', 'red', 'bold'),
            'CRITICAL': strc('CRITICAL', 'red', 'bold')}

        if record.levelname in level_colors.keys():
            record.levelname = level_colors[record.levelname]
        record.name = strc(record.name, 'black', 'bold')

        return logging.Formatter.format(self, record)

# Definitions =================================================================


def ANSIIcode(color='black', style='normal'):
    """Create an ANSII color code.

    :param color: ANSII color among ['normal','black','red','green',
                  'yellow','blue','purple','cyan','white']
    :param style: ANSII style among ['normal','bold','faint','italic',
                  'underline','blink','fast_blink']

    :return: the ANSII color code corresponding to the chosen color and
             style
    """

    colorCode = colorCodes[color]
    styleCode = styleCodes[style]

    return '\033[' + styleCode + colorCode + 'm'


def strc(text, color='black', style='normal'):
    """Create a string from a given text with corresponding color and
    style.

    :param text: text to be printed
    :param color: ANSI color among ['normal','black','red','green',
                  'yellow','blue','purple','cyan','white']
    :param style: ANSI style among ['normal','bold','faint','italic',
                  'underline','blink','fast_blink']

    :return: Colored string
    """

    ansii = ANSIIcode(color, style)
    back_to_normal = ANSIIcode('normal', 'normal')  # '\033[0m'

    return ansii + text + back_to_normal


def printc(text, color='black', style='normal', **kwargs):
    """Print a given text with corresponding color and style.

    :param text: text to be printed
    :param color: ANSI color among ['normal','black','red','green',
                  'yellow','blue','purple','cyan','white']
    :param style: ANSI style among ['normal','bold','faint','italic',
                  'underline','blink','fast_blink']
    """

    print(strc(text, color, style), **kwargs)
    sys.stdout.flush()


def writec(text, color='black', style='normal'):
    """Write a given text with corresponding color and style.

    :param text: text to be written
    :param color: ANSI color among ['normal','black','red','green',
                  'yellow','blue','purple','cyan','white']
    :param style: ANSI style among ['normal','bold','faint','italic',
                  'underline','blink','fast_blink']
    """

    sys.stdout.write(strc(text, color, style))


def logger(level, format='%(levelname)s %(message)s'):
    """Create a colorized logger (see the logging package for more
    information).

    :param level: verbose level ('DEBUG','INFO','WARNING','ERROR','CRITICAL')
    :type level: string
    :param format: Output string format
    :type format: string
    :return: Single logging channel
    :rtype: 'logging.Logger'
    """

    # Remove previous handlers
    root = logging.getLogger()
    if root.handlers:
        for handler in root.handlers:
            root.removeHandler(handler)

    # Create logger
    logger = logging.getLogger()
    logger.setLevel(getattr(logging, level.upper()))

    # Create console handler and set level to debug
    ch = logging.StreamHandler()
    ch.setLevel(getattr(logging, level.upper()))

    # Create formatter
    formatter = Formatter(format)

    # Add formatter to ch
    ch.setFormatter(formatter)

    # Add console handler to logger
    logger.addHandler(ch)

    return logger

# End of ColoredText.py =======================================================
