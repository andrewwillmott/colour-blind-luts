//
//  File:       CBLuts.h
//
//  Function:   Utilities for colour-blind modelling and LUT construction
//
//  Copyright:  Andrew Willmott 2018
//

#ifndef CB_LUTS_H
#define CB_LUTS_H

#include <stdint.h>

namespace CBLut
{
    struct Vec3f { float x; float y; float z; };
    struct Mat3f { Vec3f x; Vec3f y; Vec3f z; };

    // LMS colour space, models human eye response: https://en.wikipedia.org/wiki/LMS_color_space
    extern const Mat3f kLMSFromRGB;     ///< Convert to LMS colour system
    extern const Mat3f kRGBFromLMS;     ///< Convert back from LMS colour system

    extern const Mat3f kLMSDeuteranope; ///< Deuteranope: greens are greatly reduced (1% men)
    extern const Mat3f kLMSProtanope;   ///< Protanope: reds are greatly reduced (1% men)
    extern const Mat3f kLMSTritanope;   ///< Tritanope: blues are greatly reduced (0.003% population)

    enum tLMS
    {
        kL,
        kM,
        kS,
    };

    // Colour blindness simulation and correction
    Vec3f Simulate (Vec3f rgb, tLMS lmsType, float strength = 1.0f); ///< Simulate given form of colour blindness, with optional 0-1 strength for e.g. protanomaly (< 1) rather than protanopia (= 1). 
    
    Vec3f Daltonise(Vec3f rgb, tLMS lmsType, float strength = 1.0f); ///< "Daltonise" 'rgb' to enhance it for the given type of colour blindness, using Fidaner et al.
    Vec3f Correct  (Vec3f rgb, tLMS lmsType, float strength = 1.0f); ///< Correct image for given type of colour blindness using a mixture of amplification and hue shifting.


    // Simple 32-bit RGBA handling
    struct RGBA32
    {
        union
        {
            uint8_t  c[4];
            uint32_t u32;
        };
    };
    
    RGBA32 ToRGBA32   (Vec3f c, uint8_t alpha);
    RGBA32 ToRGBA32u  (Vec3f c);
    Vec3f  FromRGBA32 (RGBA32 rgb);
    Vec3f  FromRGBA32u(RGBA32 rgb);

    // RGB LUT support
    constexpr int kLUTBits = 5; // 32 x 32 x 32, compromise between accuracy and memory.
    constexpr int kLUTSize = 1 << kLUTBits;

    void CreateIdentityLUT(RGBA32 rgbLUT[kLUTSize][kLUTSize][kLUTSize]);    // Create identity
    void ApplyLUT      (RGBA32 rgbLUT[kLUTSize][kLUTSize][kLUTSize], int n, const RGBA32 dataIn[], RGBA32 dataOut[]); ///< Apply lut to the given image 
    void ApplyLUTNoLerp(RGBA32 rgbLUT[kLUTSize][kLUTSize][kLUTSize], int n, const RGBA32 dataIn[], RGBA32 dataOut[]); ///< Apply lut to the given image, using point sampling

    // Mono LUT support
    void ApplyMonoLUT(const RGBA32 monoLUT[256], int n, const RGBA32 dataIn[], RGBA32 dataOut[], int channel = -1);
    ///< Apply given mono->rgba ramp to either sRGB (D65) luminance, or the specified channel. 
}

#endif
