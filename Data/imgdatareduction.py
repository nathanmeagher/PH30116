#!/usr/bin/env python

"""
Reading in and manipulating fits files

"""

from astropy.io import fits
import numpy as np
import os

#------------------------------------------------------------------------
class ImageData(object):
    """
    : raw image class
    """
    def __init__(self,filen=None):

        self.headprim = None
        self.raw = None
        self.rawtrim = None
        self.bias = None
        self.chip1_bl = 0
        self.chip2_bl = 0
        self.biassub = None
        self.flat = None
        self.red = None
        if filen is not None:
            self.readfits(filen)

    
    def readfits(self, filen,ext=1):
        """
        : reads in fits file corresponding to the raw image class
        """
        # Get primary header
        self.headprim = fits.getheader("Data/%s"%filen, 0)
        # Read in data
        self.raw = fits.getdata("Data/%s"%filen,ext)
        
        # Use overscan to determine value of biaslevel (to be used in bias_subtract)
        self.chip1_bl = np.mean(self.raw[1:,0:49])
        self.chip2_bl = np.mean(self.raw[1:,2151:2199])
        
        # combine two halves of detector
        self.rawtrim = np.delete(self.raw,slice(2052,2100),0) # remove top border
        self.rawtrim = np.delete(self.rawtrim,slice(0,49),1) # remove lhs border
        self.rawtrim = np.delete(self.rawtrim,slice(2102,2151),1) # remove rhs border
        self.rawtrim = np.delete(self.rawtrim,slice(1023,1077),1) # remove central gap

    def readdata(self, filen,ext=1):
        """
        : reads in a reduced data file
        """
        # Get primary header
        self.redhead = fits.getheader("Data/%s"%filen, 0)
        # Read in data
        self.red = fits.getdata("Data/%s"%filen,ext)

    def readbias(self,bias):
        """
        : reads in bias frame
        """
        biasimg = fits.getdata("Data/%s"%bias,0)
        self.bias = biasimg+self.chip1_bl
        self.bias[:,1025:2048] = self.bias[:,1025:2048]-self.chip1_bl+self.chip2_bl
        
    def bias_subtract(self):
        """
        : subtracts bias frame from raw image correspondong to the raw image class instance
        """
        self.biassub = self.rawtrim - self.bias
        
    def readflat(self,flat):
        """
        : reads in flat-frame
        """
        flatimg = fits.getdata("Data/%s"%flat,0)
        self.flat = flatimg
      
    def flatfield(self):
        """
        : divides bias-subtracted raw image by flat-field. Function exits with warning if the bias-subtract image has not yet been produced
        """
        if np.any(self.bias):
            self.red = self.biassub/self.flat
        else:
            print("\nWARNING: you must first subtract the bias\n")
            exit()
            
        
    def savefits(self,filen,name):
        """
        : saves input image or 2d-array to a fits file with name provided as an argument
        """
        hdu = fits.HDUList()
        hdu.append(fits.PrimaryHDU(header=self.headprim))
        hdu.append(fits.ImageHDU(data=filen, header=self.headprim))
    
        fileout = 'Data/%s.fits' % name
        if os.path.isfile(fileout):
            os.remove(fileout)
        hdu.writeto(fileout)
