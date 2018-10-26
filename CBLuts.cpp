//
//  File:       CBLuts.cpp
//
//  Function:   Utilities for colour-blind modelling and LUT construction
//
//  Copyright:  Andrew Willmott 2018
//

#include "CBLuts.h"

#include "CividisLUT.h"
#include <math.h>
#include <assert.h>

using namespace CBLut;

// --- Colour-blind support ---------------------------------------------------

// LMS colour space, models human eye response: https://en.wikipedia.org/wiki/LMS_color_space

// https://ixora.io/projects/colorblindness/color-blindness-simulation-research/
// More recent version of original approach below, uses more up-to-date LMS
// conversion, different approach to colour constraints, and observes Tritanope
// conversion has an issue in that it appears to have been derived by ensuring
// blue remains constant rather than red or green...

const Mat3f CBLut::kLMSFromRGB =
{
    0.31399022,    0.63951294,    0.04649755,
    0.15537241,    0.75789446,    0.08670142,
    0.01775239,    0.10944209,    0.87256922,
};

const Mat3f CBLut::kRGBFromLMS =
{
    5.47221206,   -4.6419601,     0.16963708,
    -1.1252419,    2.29317094,   -0.1678952,
    0.02980165,   -0.19318073,    1.16364789,
};

const Mat3f CBLut::kLMSProtanope =      /// Protanope: red sensitivity is greatly reduced, reds/yellows appear darker (1% men).
{
    0, 1.05118294, -0.05116099,
    0, 1, 0,
    0, 0, 1,
};

const Mat3f CBLut::kLMSDeuteranope =    /// Deuteranope: green sensivitity is greatly reduced, no brightness issues (1% men)
{
    1, 0, 0,
    0.9513092, 0,  0.04866992,
    0, 0, 1,
};

const Mat3f CBLut::kLMSTritanope =      /// Tritanope: blue sensitivity greatly reduced (0.003% population)
{
    1, 0, 0,
    0, 1, 0,
    -0.86744736, 1.86727089, 0
};

namespace
{
    inline float& elt(      Vec3f& v, int i) { return (&v.x)[i]; } 
//  inline Vec3f& row(      Mat3f& m, int i) { return (&m.x)[i]; } 
    inline float  elt(const Vec3f& v, int i) { return (&v.x)[i]; } 
    inline Vec3f  row(const Mat3f& m, int i) { return (&m.x)[i]; } 
    inline Vec3f  col(const Mat3f& m, int i) { return { (&m.x.x)[i], (&m.y.x)[i], (&m.z.x)[i] }; } 

    inline Vec3f operator+(Vec3f a, Vec3f b) { return { a.x + b.x, a.y + b.y, a.z + b.z}; }
    inline Vec3f operator-(Vec3f a, Vec3f b) { return { a.x - b.x, a.y - b.y, a.z - b.z}; }
    inline Vec3f operator*(float s, Vec3f a) { return { s   * a.x, s   * a.y, s   * a.z}; }
    
    inline float dot      (Vec3f a, Vec3f b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
    inline Vec3f pow      (Vec3f v, float p) { return { powf(v.x, p), powf(v.y, p), powf(v.z, p) }; }

    inline Vec3f operator*(const Mat3f& m, const Vec3f& v) { return Vec3f { dot(m.x, v), dot(m.y, v), dot(m.z, v) }; }
}

// "Digital Video Colourmaps for Checking the Legibility of Displays by Dichromats", Viénot et al.
//
// Example implementations: http://www.daltonize.org,
// Unfortunately now dead (sublinks redirect to a new site), copy(?) here:
//   https://github.com/joergdietrich/daltonize/blob/master/daltonize.py

namespace
{
    // Note that unlike kLMSFromRGB, LMS are weighted, e.g., red -> (17.8, 3.4, 0.02), blue ->  
    const Mat3f kLMSFromRGBV =
    {
         { 17.8824,   43.5161,  4.11935 },
         { 3.45565,   27.1554,  3.86714 },
         { 0.0299566, 0.184309, 1.46709 },
    };

    const Mat3f kRGBFromLMSV =
    {
        {  0.0809444479f,   -0.130504409f,    0.116721066f },
        { -0.0102485335f,    0.0540193266f,  -0.113614708f },
        { -0.000365296938f, -0.00412161469f,  0.693511405f },
    };

    // These transforms to LMS colours simulate particular forms of colour blindness
    const Mat3f kLMSProtanopeV =      /// Protanope: red sensitivity is greatly reduced, reds/yellows appear darker (1% men).
    {
         { 0.0, 2.02344, -2.52581, },
         { 0.0, 1.0,      0.0,     },
         { 0.0, 0.0,      1.0      },
    };

    const Mat3f kLMSDeuteranopeV =    /// Deuteranope: green sensivitity is greatly reduced, no brightness issues (1% men)
    {
         { 1.0,      0.0, 0.0,      },
         { 0.494207, 0.0, 1.24827,  },
         { 0.0,      0.0, 1.0       },
    };

    const Mat3f kLMSTritanopeV =      /// Tritanope: blue sensitivity greatly reduced (0.003% population)
    {
         1.0,       0.0,      0.0,
         0.0,       1.0,      0.0,
         -0.395913, 0.801109, 0.0
    };
}

namespace
{
    // From Onur Fidaner, Poliang Lin, and Nevran Ozguven. http://scien.stanford.edu/class/psych221/ projects/05/ofidaner/project report.pdf.
    // Unfortunately SCIEN have seen fit to break the links on their website, and there seems to be no author copy.
    // Update: a copy can be found here for now: https://github.com/joergdietrich/daltonize/blob/master/doc/project_report.pdf
    // The matrix values can be found here: https://github.com/joergdietrich/daltonize/blob/master/doc/conv_img.m.
    // Their precise values aren't discussed or justified in the paper.
    const Mat3f kDaltonErrorToDeltaP =
    {
        { 0.0, 0.0, 0.0, },
        { 0.7, 1.0, 0.0, },
        { 0.7, 0.0, 1.0  },
    };
    const Mat3f kDaltonErrorToDeltaD =
    {
        { 1.0, 0.7, 0.0, },
        { 0.0, 0.0, 0.0, },
        { 0.0, 0.7, 1.0  },
    };
    const Mat3f kDaltonErrorToDeltaT =
    {
        { 1.0, 0.0, 0.7, },
        { 0.0, 1.0, 0.7, },
        { 0.0, 0.0, 0.0  },
    };

    Vec3f SimulateV(Vec3f rgb, const Mat3f& lmsTransform)
    {
        Vec3f lms = lmsTransform * (kLMSFromRGBV * rgb);
        return kRGBFromLMSV * lms;
    }
}

Vec3f CBLut::Daltonise(Vec3f rgb, tLMS lmsType, float strength)
{
    // Daltonisation: take delta from original RGB
    // + use to shift colors towards visible spectrum
    Vec3f rgbSim;
    Vec3f rgbDelta;
    
    switch (lmsType)
    {
    case kL:
        rgbSim = SimulateV(rgb, kLMSProtanopeV);
        rgbDelta = kDaltonErrorToDeltaP * (strength * (rgb - rgbSim));
        break;
    case kM:
        rgbSim = SimulateV(rgb, kLMSDeuteranopeV);
        rgbDelta = kDaltonErrorToDeltaD * (strength * (rgb - rgbSim));
        break;
    case kS:
        rgbSim = SimulateV(rgb, kLMSTritanopeV);
        rgbDelta = kDaltonErrorToDeltaT * (strength * (rgb - rgbSim));
        break;
    }
    
    return rgb + rgbDelta;
}

namespace
{
    const Mat3f kLMSSimulate =      // P/D/T simulation amalgamated into one matrix
    {
        0,           1.05118294, -0.05116099,
        0.9513092,   0,           0.04866992,
        -0.86744736, 1.86727089,  0
    };

    const Mat3f kNCDeltaBrighten =     // vanilla transfer error to remaining channels, like Daltonise approach but in lms space. amount=2.5 gets closish match
    {
        0, 1, 1,
        1, 0, 1,
        1, 1, 0,
    };

    const Mat3f kNCDeltaContrast =     // vanilla increase contrast between remaining channels. Not very controllable.
    {
         0, +1, +1,
        -1,  0, -1,
        +1, -1,  0,
    };

    const Mat3f kNCDeltaMix1 =          // mix of contrast/brighten
    {
        0,   0.2, 0.8,
        0.2, 0,   0.2,
        0.8, 0.8, 0,
    };

    // const Mat3f kLMSDaltonError = kLMSFromRGB * kDaltonErrorToDeltaP * kRGBFromLMS;
    const Mat3f kNCDalton =         // The Fidaner RGB delta transformed to LMS space.
    {
        { 1.90957534, -0.771573185, 0.0281965993f },
        { 2.38503432, -1.02317369,  0.0739354342f },
        { 3.66449642, -3.10851407,  1.11359835f   },
    };

    const Mat3f kNCDeltaRecip =         // trans(1/kLMSSimulate)
    {
        0,             1.05118299, -1.15280771,
        0.951309144,    0,          0.535540938,
        -19.5461426, 20.5465717,    0,
    };

    const Mat3f kNCDeltaRecipAbs =      // abs(trans(1/kLMSSimulate))
    {
        0,             1.05118299, 1.15280771,
        0.951309144,    0,         0.535540938,
        19.5461426, 20.5465717,    0,
    };

    const Mat3f kNCDeltaInv =     // inv(kSimulate)
    {
        0.672, 0.706, -0.378,
        0.312, 0.328, 0.360,
        -13.133, 6.741, 7.393
    };
}

Vec3f CBLut::Simulate(Vec3f rgb, tLMS lmsType, float strength)
{
    Vec3f lms = kLMSFromRGB * rgb;

    float&      eltx   = elt(lms, lmsType);                     // affected channel
    const float simElt = dot(row(kLMSSimulate, lmsType), lms);  // 'sim' weighted combo of the other two.

    eltx += strength * (simElt - eltx);

    rgb = kRGBFromLMS * lms;

    return rgb;
}

Vec3f CBLut::Correct(Vec3f rgb, tLMS lmsType, float strength)
{
    const Vec3f lms = kLMSFromRGB * rgb;

    const float orgElt = elt(lms, lmsType);                     // original value of affected channel
    const float simElt = dot(row(kLMSSimulate, lmsType), lms);  // simulated full-strength value
    const float error  = strength * (orgElt - simElt);          // error

    float mc = strength * strength; // How much to use strategy 1: redistributing error into other channels in a way that shifts hue
    float ms = 1.0f - strength;     // How much to use stragegy 2: simply brighten affected channel

    Vec3f amount3Recip = { -0.25f, -0.3f, -0.07f };    // tuning values for redistribution
    float amount = elt(amount3Recip, lmsType);

    Vec3f correct = mc * amount * col(kNCDeltaRecip, lmsType);
    elt(correct, lmsType) = ms * 2.0f;

    Vec3f lmsCorrect = lms + error * correct;
    
    rgb = kRGBFromLMS * lmsCorrect;

    return rgb;
}


// --- LUT support ------------------------------------------------------------

namespace
{
    constexpr float kGamma = 2.2f;

    inline uint8_t ToU8(float f)
    {
        if (f <= 0.0f)
            return 0;
        if (f >= 1.0f)
            return 255;

        return uint8_t(f * 255.0f + 0.5f);
    }

    inline uint8_t ToU8u(float f)   // 0-256 variant used for LUT construction
    {
        if (f <= 0.0f)
            return 0;
        if (f >= 1.0f)
            return 255;

        return uint8_t(f * 256.0f);
    }
}

RGBA32 CBLut::ToRGBA32(Vec3f c)
{
    c = pow(c, 1.0f / kGamma);
    RGBA32 result;

    result.c[0] = ToU8(c.x);
    result.c[1] = ToU8(c.y);
    result.c[2] = ToU8(c.z);
    result.c[3] = 255;

    return result;
}

RGBA32 CBLut::ToRGBA32u(Vec3f c)
{
    c = pow(c, 1.0f / kGamma);
    RGBA32 result;

    result.c[0] = ToU8u(c.x);
    result.c[1] = ToU8u(c.y);
    result.c[2] = ToU8u(c.z);
    result.c[3] = 255;

    return result;
}

Vec3f CBLut::FromRGBA32(RGBA32 rgb)
{
    Vec3f c = { rgb.c[0] / 255.0f, rgb.c[1] / 255.0f, rgb.c[2] / 255.0f };
    return pow(c, kGamma);
}

Vec3f CBLut::FromRGBA32u(RGBA32 rgb)
{
    Vec3f c = { rgb.c[0] / 256.0f, rgb.c[1] / 256.0f, rgb.c[2] / 256.0f };
    return pow(c, kGamma);
}
    

void CBLut::CreateCividisLUT(RGBA32 linearLUT[256])
{
    for (int i = 0; i < 256; i++)
    {
        Vec3f c0 = { kCividisLUT[i][2], kCividisLUT[i][1], kCividisLUT[i][0] };
        
        linearLUT[i] = ToRGBA32(c0);
    }
}

void CBLut::CreateIdentityLUT(RGBA32 rgbLUT[kLUTSize][kLUTSize][kLUTSize])
{
    constexpr int scale  = 256 / kLUTSize;
    constexpr int offset = scale / 2;

    for (int i = 0; i < kLUTSize; i++)
    for (int j = 0; j < kLUTSize; j++)
    for (int k = 0; k < kLUTSize; k++)
    {
        RGBA32& p = rgbLUT[i][j][k];

        p.c[0] = k * scale + offset;
        p.c[1] = j * scale + offset;
        p.c[2] = i * scale + offset;
        p.c[3] = 255;
    }
}

#define EXTRAPOLATE_LUT 1

void CBLut::ApplyLUT(RGBA32 rgbLUT[kLUTSize][kLUTSize][kLUTSize], int n, const RGBA32 dataIn[], RGBA32 dataOut[])
{
    constexpr int lutShift = kLUTBits;
    constexpr int lutSize  = 1 << lutShift;
    constexpr int fShift   = 8 - lutShift;
    constexpr int fHalf    = 1 << (fShift - 1);
    constexpr int fMask    = (1 << fShift) - 1;

    for (int i = 0; i < n; i++)
    {
        const uint8_t* ci = dataIn[i].c;

        int co[3] = { ci[0] + fHalf,   ci[1] + fHalf,   ci[2] + fHalf   };
        int i1[3] = { co[0] >> fShift, co[1] >> fShift, co[2] >> fShift };
        int i0[3] = { i1[0] - 1,       i1[1] - 1,       i1[2] - 1       };
        int s [3] = { co[0] & fMask,   co[1] & fMask,   co[2] & fMask   };

        for (int j = 0; j < 3; j++)
        {
            if (i0[j] < 0)
            {
                i0[j]++;
            #ifdef EXTRAPOLATE_LUT
                i1[j]++;
                s [j] -= 8;
            #endif
            }
            else
            if (i1[j] >= lutSize)
            {
                i1[j]--;
            #ifdef EXTRAPOLATE_LUT
                i0[j]--;
                s [j] += 8;
            #endif
            }

            assert(0 <= i0[j] && i0[j] < kLUTSize);
            assert(0 <= i1[j] && i1[j] < kLUTSize);
        }
        
        RGBA32 lutC0 = rgbLUT[i0[2]][i0[1]][i0[0]];
        RGBA32 lutC1 = rgbLUT[i1[2]][i1[1]][i1[0]];

        int ch0 = (((8 - s[0]) * lutC0.c[0] + s[0] * lutC1.c[0])) >> 3;
        int ch1 = (((8 - s[1]) * lutC0.c[1] + s[1] * lutC1.c[1])) >> 3;
        int ch2 = (((8 - s[2]) * lutC0.c[2] + s[2] * lutC1.c[2])) >> 3;

    #ifdef EXTRAPOLATE_LUT
        ch0 = ch0 < 0 ? 0 : ch0 > 255 ? 255 : ch0;
        ch1 = ch1 < 0 ? 0 : ch1 > 255 ? 255 : ch1;
        ch2 = ch2 < 0 ? 0 : ch2 > 255 ? 255 : ch2;
    #endif

        assert(0 <= ch0 && ch0 <= 255);
        assert(0 <= ch1 && ch1 <= 255);
        assert(0 <= ch2 && ch2 <= 255);

        dataOut[i].c[0] = ch0;
        dataOut[i].c[1] = ch1;
        dataOut[i].c[2] = ch2;
        dataOut[i].c[3] = 255;
    }
}

void CBLut::ApplyLUTNoLerp(RGBA32 rgbLUT[kLUTSize][kLUTSize][kLUTSize], int n, const RGBA32 dataIn[], RGBA32 dataOut[])
{
    constexpr int fShift = 8 - kLUTBits;

    for (int i = 0; i < n; i++)
    {
        const uint8_t* ci = dataIn[i].c;

        dataOut[i] = rgbLUT[ci[2] >> fShift][ci[1] >> fShift][ci[0] >> fShift];
    }
}

