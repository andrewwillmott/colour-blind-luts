#ifndef COLOUR_MAPS_H
#define COLOUR_MAPS_H

// Various "false colour" mono->rgba maps for visualisation, 
// assembled in gamma-corrected RGBA32 format by Andrew Willmott. 
// See particularly the viridis and cividis maps for colourblind-savvy
// ramps.

// kMagmaLUT, kInfernoLUT, kPlasmaLUT, kViridisLUT
//
// New matplotlib colormaps by Nathaniel J. Smith, Stefan van der Walt,
// and (in the case of viridis) Eric Firing.
//
// This file and the colormaps in it are released under the CC0 license /
// public domain dedication. We would appreciate credit if you use or
// redistribute these colormaps, but do not impose any legal restrictions.
//
// To the extent possible under law, the persons who associated CC0 with
// mpl-colormaps have waived all copyright and related or neighboring rights
// to mpl-colormaps.
//
// You should have received a copy of the CC0 legalcode along with this
// work.  If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
//
// See e.g. https://raw.githubusercontent.com/BIDS/colormap/master/colormaps.py

extern const unsigned char kMagmaLUT  [256][4];
extern const unsigned char kInfernoLUT[256][4];
extern const unsigned char kPlasmaLUT [256][4];
extern const unsigned char kViridisLUT[256][4];

// kCividisLUT
//
// This is a version of Viridis optimised for CVD. It works better for those
// with (red-green) CVD, but is not as aesthetically pleasing as Viridis, due
// to being a straight blue/yellow ramp. 
//
// See http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0199239

extern const unsigned char kCividisLUT[256][4];

#endif
